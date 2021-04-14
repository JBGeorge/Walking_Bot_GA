#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <GL/glut.h>
#include <vector>
#include <stdarg.h>
#include <fstream>
#include "omp.h"
#include <algorithm>
#include <math.h>
#include <numeric>
#include <random>
#include <queue>
#include <ctime>
#include <string>


#define PI 3.1412
#define Cos(th) cos(PI/180*(th))
#define Sin(th) sin(PI/180*(th))
#define NUM_THREADS 16

struct MASS
{
    double m;      
    double pos[3];    
    double vel[3];    
    double acc[3];    
    int springFactor;
};

struct SPRING
{
    double k;       
    double L_0;     
    double L;
    int m1;         
    int m2;        
    int type;
};

struct GENE
{
    bool exist;
    int springFactor;
    std::vector<double> p;
};



#define location 3
#define depth 4
#define W 4
#define pMutation 0.9
#define GRAPHICS
#define UNIQUEPOSITION


int th = 0;            
int ph = 0;           
int axes = 0;         
int light = 0;
double asp = 1;     
int fov = 40;         
double dim = 10.0;  

int generationNumber = 1;
int cubotNumber = 1;
int simulationTime = 12;

double mass = 0.5;
double length = 0.8;
double gravity = 9.8;
double T = 0;

double timeStep = 0.001;
double restoreConstant = 100000;
double springConstant = 10000;
double dampingConstant = 0.99;
double frictionCoefficient = 0.8;

static GLint Frames = 0;
static GLfloat fps = -1;
static GLint T0 = 0;

double start_time;

GLfloat worldRotation[16] = { 1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1 };


double calcmag(double x[], std::size_t sz)
{
    return std::sqrt(std::inner_product(x, x + sz, x, 0.0));
}

std::vector<int> sort_indexes(const std::vector<double>& v) {

    std::vector<int> idx(v.size());
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(),
        [&v](int i1, int i2) {return v[i1] > v[i2]; });
    return idx;
}

std::vector<int> generateChildrenIndex(int crossOverPoint) {
    std::vector<int> childrenIndex;
    childrenIndex.push_back(crossOverPoint);
    std::queue<int> myQueue;
    myQueue.push(crossOverPoint);
    while (!myQueue.empty()) {
        int index = myQueue.front();
        if (3 * index + 3 <= (pow(3, depth) - 1) / 2) {
            int indexArray[3] = { 3 * index + 1, 3 * index + 2, 3 * index + 3 };
            for (int i = 0; i < 3; i++) {
                childrenIndex.push_back(indexArray[i]);
                myQueue.push(indexArray[i]);
            }
        }
        myQueue.pop();
    }
    return childrenIndex;
}

std::vector<int> generateCrossOverRange(int crossOverLevel) {
    std::vector<int> crossOverRange;
    int start = (pow(3, crossOverLevel + 1) - 1) / 2;
    int end = (pow(3, crossOverLevel + 2) - 1) / 2;
    for (int i = start; i < end; i++) {
        crossOverRange.push_back(i);
    }
    return crossOverRange;
}

int randomNumber(int min, int max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(min, max);
    return dis(gen);
}


int randomSpringFactor()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 3);
    return dis(gen);
}

class cubot
{
public:
    std::vector<GENE> cubotGene;
    std::vector<MASS> cubemasses;
    std::vector<SPRING> cubesprings;
    double initialLocation[3] = { 0,0,0 };

    cubot(double initialX, double initialY, double initialZ, std::vector<GENE> cubotGenes) {
        cubotGene = cubotGenes;
        generatecubot(initialX, initialY, initialZ);
        generateSpring();
        calculateInitialLocation();
    }

    void calculateInitialLocation() {
        double x = 0; double y = 0;
        for (int j = 0; j < cubemasses.size(); j++) {
            x = x + cubemasses[j].pos[0];
            y = y + cubemasses[j].pos[1];
        }
        x = x / cubemasses.size();
        y = y / cubemasses.size();
        initialLocation[0] = x; initialLocation[1] = y;
    }


    void generateSpring() {

        for (int i = 0; i < cubemasses.size() - 1; i++) {
            for (int j = i + 1; j < cubemasses.size(); j++) {
                double positionDiff[3] = { cubemasses[j].pos[0] - cubemasses[i].pos[0],cubemasses[j].pos[1] - cubemasses[i].pos[1],cubemasses[j].pos[2] - cubemasses[i].pos[2] };
                if (calcmag(positionDiff, 3) < 2 * length) {
                    cubesprings.push_back({ springConstant,calcmag(positionDiff,3),calcmag(positionDiff,3),i,j,(cubemasses[i].springFactor + cubemasses[j].springFactor) / 2 });
                }
            }
        }
    }

    void generatecubot(double iX, double iY, double iZ) {
        std::queue<int> myQueue;
        myQueue.push(0);

        while (!myQueue.empty()) {
            int index = myQueue.front();
            cubemasses.push_back(
                { 0.5, {iX + cubotGene[index].p[0], iY + cubotGene[index].p[1], iZ + cubotGene[index].p[2]}, {0, 0, 0}, {0, 0, 0}, cubotGene[index].springFactor });

            if (3 * index + 3 <= (pow(3, depth) - 1) / 2) {
                int indexArray2[3] = { 3 * index + 1, 3 * index + 2, 3 * index + 3 };
                for (int i = 0; i < 3; i++) {
                    if (cubotGene[indexArray2[i]].exist == true) {
                        myQueue.push(indexArray2[i]);
                    }
                }
            }
            myQueue.pop();
        }


    }


    void cubotDraw() {
        glColor3f(1, 0, 0);

        GLUquadric* quad;
        quad = gluNewQuadric();
        for (int i = 0; i < (int)cubemasses.size(); i++) {
            glPushMatrix();
            glMultMatrixf(worldRotation);
            glTranslated(cubemasses[i].pos[0], cubemasses[i].pos[1], cubemasses[i].pos[2]);
            gluSphere(quad, length / 20, 10, 10);
            glPopMatrix();
        }

        for (int i = 0; i < (int)cubesprings.size(); i++) {
            if (cubesprings[i].type == 0) { glColor3f(0.0, 0.7, 0.8); }
            if (cubesprings[i].type == 1) { glColor3f(0.3, 0.4, 0.1); }
            if (cubesprings[i].type == 2) { glColor3f(0.1, 0.9, 0.5); }
            if (cubesprings[i].type == 3) { glColor3f(0.6, 0.7, 1.0); }
            glPushMatrix();
            glMultMatrixf(worldRotation);
            glBegin(GL_LINES);
            glLineWidth(1000);
            glVertex3f(GLfloat(cubemasses[cubesprings[i].m1].pos[0]), GLfloat(cubemasses[cubesprings[i].m1].pos[1]),
                GLfloat(cubemasses[cubesprings[i].m1].pos[2]));
            glVertex3f(GLfloat(cubemasses[cubesprings[i].m2].pos[0]), GLfloat(cubemasses[cubesprings[i].m2].pos[1]),
                GLfloat(cubemasses[cubesprings[i].m2].pos[2]));
            glPopMatrix();
        }
        double x = 0; double y = 0;
        for (int j = 0; j < cubemasses.size(); j++) {
            x = x + cubemasses[j].pos[0];
            y = y + cubemasses[j].pos[1];
        }
        x = x / cubemasses.size();
        y = y / cubemasses.size();
        glPushMatrix();
        glMultMatrixf(worldRotation);
        glBegin(GL_LINES);
        glColor3f(0.0, 1.0, 0.0);
        glVertex3f(GLfloat(x), GLfloat(y), GLfloat(0.0));
        glVertex3f(GLfloat(initialLocation[0]), GLfloat(initialLocation[1]), GLfloat(0.0));
        glEnd();
        glPopMatrix();

    }

    void cubotUpdate()
    {
        std::vector<std::vector<double>> cubeforces((int)cubemasses.size(), std::vector<double>(3));
        for (int i = 0; i < (int)cubesprings.size(); i++) {
            if (cubesprings[i].type == 0) {
                cubesprings[i].k = 2000;
            }
            else if (cubesprings[i].type == 1) {
                cubesprings[i].k = 1000;
                cubesprings[i].L = cubesprings[i].L_0 - 0.3 * length * sin(W * T);
            }
            else if (cubesprings[i].type == 2) {
                cubesprings[i].k = 1000;
                cubesprings[i].L = cubesprings[i].L_0 + 0.3 * length * sin(W * T);
            }
            else if (cubesprings[i].type == 3) {
                cubesprings[i].k = 1000;
            }
            MASS mass1 = cubemasses[cubesprings[i].m1];
            MASS mass2 = cubemasses[cubesprings[i].m2];
            double positionDiff[3] = { mass2.pos[0] - mass1.pos[0], mass2.pos[1] - mass1.pos[1], mass2.pos[2] - mass1.pos[2] };
            double L = calcmag(positionDiff, 3);
            double force = cubesprings[i].k * fabs(cubesprings[i].L - L);
            double direction[3];
            if (L == 0) {
                direction[0] = 0; direction[1] = 0; direction[2] = 0;
            }
            else {
                direction[0] = positionDiff[0] / L; direction[1] = positionDiff[1] / L; direction[2] = positionDiff[2] / L;
            }
            if (L > cubesprings[i].L) {
                cubeforces[cubesprings[i].m1][0] = cubeforces[cubesprings[i].m1][0] + direction[0] * force;
                cubeforces[cubesprings[i].m1][1] = cubeforces[cubesprings[i].m1][1] + direction[1] * force;
                cubeforces[cubesprings[i].m1][2] = cubeforces[cubesprings[i].m1][2] + direction[2] * force;
                cubeforces[cubesprings[i].m2][0] = cubeforces[cubesprings[i].m2][0] - direction[0] * force;
                cubeforces[cubesprings[i].m2][1] = cubeforces[cubesprings[i].m2][1] - direction[1] * force;
                cubeforces[cubesprings[i].m2][2] = cubeforces[cubesprings[i].m2][2] - direction[2] * force;
            }
           
            else if (L < cubesprings[i].L) {
                cubeforces[cubesprings[i].m1][0] = cubeforces[cubesprings[i].m1][0] - direction[0] * force;
                cubeforces[cubesprings[i].m1][1] = cubeforces[cubesprings[i].m1][1] - direction[1] * force;
                cubeforces[cubesprings[i].m1][2] = cubeforces[cubesprings[i].m1][2] - direction[2] * force;
                cubeforces[cubesprings[i].m2][0] = cubeforces[cubesprings[i].m2][0] + direction[0] * force;
                cubeforces[cubesprings[i].m2][1] = cubeforces[cubesprings[i].m2][1] + direction[1] * force;
                cubeforces[cubesprings[i].m2][2] = cubeforces[cubesprings[i].m2][2] + direction[2] * force;
            }
        }

        for (int i = 0; i < (int)cubemasses.size(); i++) {
        
            cubeforces[i][2] = cubeforces[i][2] - cubemasses[i].m * gravity;
            
            if (cubemasses[i].pos[2] <= 0) {
                cubeforces[i][2] = cubeforces[i][2] + restoreConstant * fabs(cubemasses[i].pos[2]);
             
                double F_h = sqrt(pow(cubeforces[i][0], 2) + pow(cubeforces[i][1], 2));
                double F_v = cubeforces[i][2];
                if (F_h < F_v * frictionCoefficient) {
                    cubeforces[i][0] = 0;
                    cubeforces[i][1] = 0;
                    cubemasses[i].vel[0] = 0;
                    cubemasses[i].vel[1] = 0;
                }
            }
          
            cubemasses[i].acc[0] = cubeforces[i][0] / cubemasses[i].m;
            cubemasses[i].acc[1] = cubeforces[i][1] / cubemasses[i].m;
            cubemasses[i].acc[2] = cubeforces[i][2] / cubemasses[i].m;

            cubemasses[i].vel[0] = dampingConstant * (cubemasses[i].vel[0] + cubemasses[i].acc[0] * timeStep);
            cubemasses[i].vel[1] = dampingConstant * (cubemasses[i].vel[1] + cubemasses[i].acc[1] * timeStep);
            cubemasses[i].vel[2] = dampingConstant * (cubemasses[i].vel[2] + cubemasses[i].acc[2] * timeStep);

            cubemasses[i].pos[0] = cubemasses[i].pos[0] + cubemasses[i].vel[0] * timeStep;
            cubemasses[i].pos[1] = cubemasses[i].pos[1] + cubemasses[i].vel[1] * timeStep;
            cubemasses[i].pos[2] = cubemasses[i].pos[2] + cubemasses[i].vel[2] * timeStep;
        }
    }
};

class Simulation {
private:
    int populationSize;
    std::vector<double> populationDistance;
    std::vector<std::vector<GENE>> populationGene;
    std::vector<std::vector<GENE>> newPopulationGene;
    std::vector<cubot> cubots;
    std::vector<std::vector<double>> massPositions;
    std::ofstream bestGene;
    std::ofstream popDis;
public:
    double averageDistance;
    double maxDistance;

    Simulation(int popSize) {
        populationSize = popSize;
        //generateGenes();
        generateBestGene();
        generatecubots();
        std::time_t t = std::time(0);
        popDis.open("PopulationMaxDist.txt");
        bestGene.open("BestGenomePop");
    }
    void startSim(double time) {
        if (T < time || T > time) {
            simUpdate();
#ifdef  GRAPHICS
            simDraw();
#endif
        }
        else {
            printf("### Generation %d ###", generationNumber);
            double time = omp_get_wtime() - start_time;
            printf(" Time: %f ###\n", time);
            calculatePopulationDistance();
            generationNumber++;
            cubots.clear();
            cubots.shrink_to_fit();
            generatecubots();
            T = 0;
            start_time = omp_get_wtime();
        }
    }

    void generateBestGene() {
        int gene00[160] = { 1,3,2,0,0,2,3,0,0,2,3,0,0,3,1,3,0,1,2,0,1,1,3,3,2,3,1,1,0,3,1,0,1,2,2,1,2,1,3,0,0,1,2,0,1,1,3,3,2,3,2,3,3,1,3,1,2,1,0,3,1,0,0,0,1,3,2,1,0,1,2,1,0,1,3,0,1,1,3,2,3,2,3,2,1,3,3,3,3,1,3,1,0,0,2,0,1,3,3,0,1,3,2,1,0,1,3,1,0,0,2,0,2,3,2,2,3,2,3,2,1,1,1,1,3,1,3,1,2,0,3,2,1,0,0,0,1,2,2,2,0,1,2,1,0,0,2,0,1,1,3,2,3,2,3,2,1,3,3,3};
        int gene98[52] = { 0,1,1,0,1,1,1,1,0,2,1,0,2,2,2,0,3,1,0,2,1,1,0,0,0,1,2,0,1,1,0,0,3,0,1,2,3,1,1,2 };
        int gene99[52] = { 0,1,1,0,2,1,2,1,0,2,1,0,2,2,2,0,3,1,0,2,1,3,0,0,0,1,2,0,2,1,0,0,3,0,1,2,3,1,2,2 };
        int gene190[52] = { 0,1,1,0,1,1,1,1,1,1,0,0,1,1,0,0,3,1,0,2,1,2,2,2,3,1,0,2,3,1,0,2,1,1,0,0,1,1,0,0,3,1,0,2,1,1,0,0,3,1,0,2 };
        int gene299[52] = { 0,1,1,0,1,1,1,1,1,1,1,1,1,1,2,2,3,1,0,2,1,1,0,0,1,1,0,0,3,1,0,2,1,1,0,0,3,1,0,2,3,1,2,1,1,1,0,0,1,1,0,0 };
        int gene417[52] = { 0,1,1,0,1,1,1,1,1,1,2,2,1,1,2,2,3,1,0,2,1,1,0,0,3,1,0,2,3,1,2,1,1,1,1,1,1,1,0,1,3,1,2,1,1,1,0,0,1,1,0,0 };
        int gene500[52] = { 0,1,1,0,1,1,0,1,1,1,2,2,1,1,2,2,3,1,0,2,1,1,0,0,1,1,2,2,1,1,2,2,1,1,2,2,1,1,0,0,3,1,2,1,1,1,2,2,1,1,1,1 };
        int gene1708[160] = { 1,2,2,1,2,3,2,0,2,3,2,0,3,0,3,0,0,2,0,1,0,3,1,1,1,1,3,1,0,2,0,1,1,2,0,0,3,0,3,0,2,0,3,1,3,2,1,0,3,1,3,1,0,2,1,2,2,1,2,2,3,1,3,1,0,1,2,1,2,1,0,0,3,0,3,0,0,3,1,1,2,1,0,0,0,3,1,1,0,2,0,1,0,2,0,1,3,0,3,0,0,1,3,2,3,0,3,0,2,1,0,1,3,3,1,1,2,1,3,0,3,1,3,1,0,1,3,0,2,1,2,2,1,1,1,1,0,0,3,0,0,1,2,1,1,3,2,1,0,0,0,1 };
        int gene1600[484] = { 0,4,4,1,2,4,3,3,3,4,0,1,2,2,3,3,3,4,0,1,3,1,0,1,1,1,4,3,3,0,2,4,3,2,0,1,3,3,4,2,3,0,1,4,1,1,4,3,1,2,1,1,0,3,1,0,2,3,4,2,2,4,2,0,3,3,4,2,3,4,2,1,2,0,3,3,2,1,3,2,2,0,2,0,0,4,0,0,0,1,0,3,2,3,2,0,2,1,0,2,2,3,1,2,3,4,0,1,1,4,1,4,0,1,0,3,2,3,3,0,2,1,1,4,0,2,0,1,2,3,2,0,2,2,0,2,2,4,1,1,0,4,0,0,0,0,1,3,1,2,2,2,0,1,1,1,0,2,1,3,0,0,4,3,0,0,3,4,3,1,0,3,2,3,3,0,1,0,3,2,2,3,1,0,1,1,0,1,1,3,3,4,0,2,0,2,3,1,4,3,1,3,3,4,3,0,4,2,1,0,0,4,1,1,4,4,1,2,4,2,2,3,3,0,1,4,2,1,2,1,1,2,2,0,3,1,2,3,1,0,2,4,2,1,2,3,3,3,1,3,2,1,0,0,0,4,2,2,3,2,0,1,3,4,2,0,1,1,0,0,0,0,0,4,3,4,3,3,3,0,2,4,1,2,1,4,3,0,2,2,2,3,2,4,3,0,1,3,2,1,0,1,2,3,2,3,3,0,1,4,2,1,2,1,0,2,1,0,4,3,1,0,0,4,1,2,3,3,2,0,3,1,3,3,3,0,2,3,1,0,0,0,4,3,0,0,2,1,3,3,0,1,2,3,3,0,2,0,4,4,3,4,3,3,0,4,0,1,1,1,1,4,0,3,4,3,0,4,2,2,0,0,1,3,0,1,3,2,2,3,3,0,1,0,0,4,2,3,1,1,2,4,2,2,1,2,3,4,0,1,3,0,2,3,3,0,1,4,1,1,2,0,3,4,1,4,1,3,2,1,2,4,3,3,0,0,2,2,3,2,0,1,3,4,2,1,0,2,2,1,0,1,0,3,3,2,0,1,1,1,1,1,4,2,2,0,1,1,3,4,2,2,2,0,2,0,0,4,2,4,0,1,2,3 };
        std::vector<GENE> temp;
        std::vector<GENE> temp1;
        std::vector<GENE> temp2;
        std::vector<GENE> temp3;
        std::vector<GENE> temp4;
        std::vector<GENE> temp5;
        std::vector<GENE> temp6;
        std::vector<GENE> temp7;
        std::vector<GENE> temp8;
        for (int i = 0; i < 484; i = i + 4) {
            std::vector<double> p{ (double)gene1708[i + 1], (double)gene1708[i + 2], (double)gene1708[i + 3] };
            temp.push_back({ true, gene1708[i], p });
        }
        populationGene.push_back(temp);
        /*for (int i = 0; i < 484; i = i + 4) {
            std::vector<double> p{ (double)gene99[i + 1], (double)gene99[i + 2], (double)gene99[i + 3] };
            temp2.push_back({ true, gene99[i], p });
        }
        populationGene.push_back(temp2);
        for (int i = 0; i < 484; i = i + 4) {
            std::vector<double> p{ (double)gene299[i + 1], (double)gene299[i + 2], (double)gene299[i + 3] };
            temp3.push_back({ true, gene299[i], p });
        }
        populationGene.push_back(temp3);
        for (int i = 0; i < 484; i = i + 4) {
            std::vector<double> p{ (double)gene1600[i + 1], (double)gene1600[i + 2], (double)gene1600[i + 3] };
            temp4.push_back({ true, gene1600[i], p });
        }
        populationGene.push_back(temp4);*/

        //for (int i = 0; i < 484; i = i + 4) {
        //    std::vector<double> p{ (double)gene1708[i + 1], (double)gene1708[i + 2], (double)gene1708[i + 3] };
        //    temp5.push_back({ true, gene1708[i], p });
        //}
        //populationGene.push_back(temp5);

        //for (int i = 0; i < 484; i = i + 4) {
        //    std::vector<double> p{ (double)gene1708[i + 1], (double)gene1708[i + 2], (double)gene1708[i + 3] };
        //    temp6.push_back({ true, gene1708[i], p });
        //}
        //populationGene.push_back(temp6);

        //for (int i = 0; i < 484; i = i + 4) {
        //    std::vector<double> p{ (double)gene1708[i + 1], (double)gene1708[i + 2], (double)gene1708[i + 3] };
        //    temp7.push_back({ true, gene1708[i], p });
        //}
        //populationGene.push_back(temp7);

        //for (int i = 0; i < 484; i = i + 4) {
        //    std::vector<double> p{ (double)gene1708[i + 1], (double)gene1708[i + 2], (double)gene1708[i + 3] };
        //    temp8.push_back({ true, gene1708[i], p });
        //}
        //populationGene.push_back(temp8);



    }

    int x = populationGene.size();
   

    void selection() {
        std::vector<int> index = sort_indexes(populationDistance);
        newPopulationGene.clear();
        newPopulationGene.shrink_to_fit();
        for (int i = 0; i < index.size() / 2; i++) {
            newPopulationGene.push_back(populationGene[index[i]]);
        }
        for (int i = 0; i < newPopulationGene[0].size(); i++) {
            bestGene << newPopulationGene[0][i].springFactor << " " << newPopulationGene[0][i].p[0] << " " << newPopulationGene[0][i].p[1] << " " << newPopulationGene[0][i].p[2] << " ";
        }
        bestGene << newPopulationGene[0].size();
        bestGene << "\n";
        bestGene.flush();
    }

    void crossOver() {
        for (int n = 0; n < populationSize / 4; n++) {
            int parentIndex1 = randomNumber(0, newPopulationGene.size() - 1);
            int parentIndex2 = randomNumber(0, newPopulationGene.size() - 1);
            std::vector<GENE> parent1 = newPopulationGene[parentIndex1];
            std::vector<GENE> parent2 = newPopulationGene[parentIndex2];
            int maxLevel = std::max(log(parent1.size() * 2 - 1) / log(3), log(parent2.size() * 2 - 1) / log(3));
            int crossOverLevel = randomNumber(1, maxLevel);
            std::vector<int> crossOverRange = generateCrossOverRange(crossOverLevel - 1);
            int crossOverPoint1 = crossOverRange[randomNumber(0, crossOverRange.size() - 1)];
            int crossOverPoint2 = crossOverRange[randomNumber(0, crossOverRange.size() - 1)];
            std::vector<int> crossOverLocations1 = generateChildrenIndex(crossOverPoint1);
            std::vector<int> crossOverLocations2 = generateChildrenIndex(crossOverPoint2);
            std::vector<GENE> subHeap1; std::vector<GENE> subHeap2;
            for (int i = 0; i < crossOverLocations1.size(); i++) {
                subHeap1.push_back(parent1[crossOverLocations1[i]]);
            }
            for (int i = 0; i < crossOverLocations2.size(); i++) {
                subHeap2.push_back(parent2[crossOverLocations2[i]]);
            }
            std::vector<GENE> offSpring1(parent1);
            std::vector<GENE> offSpring2(parent2);
            for (int i = 0; i < crossOverLocations1.size(); i++) {
                offSpring1[crossOverLocations1[i]] = subHeap2[i];
                offSpring2[crossOverLocations2[i]] = subHeap1[i];
            }
            offSpring1 = mutation(offSpring1);
            offSpring2 = mutation(offSpring2);
            newPopulationGene.push_back(offSpring1);
            newPopulationGene.push_back(offSpring2);
        }
    }

    std::vector<GENE> mutation(std::vector<GENE> offSpring) {
        for (int i = 0; i < offSpring.size(); i++) {
            float r = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 1.0));
            if (r > 0.9) {
                offSpring[i].p[0] = randomNumber(0, location);
                offSpring[i].p[1] = randomNumber(0, location);
                offSpring[i].p[2] = randomNumber(0, location);
            }
        }
        return offSpring;
    }

    std::vector<double> generatePosition(int MAX) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, MAX);
        std::vector<double> position;
#ifndef UNIQUEPOSITION
        for (int i = 0; i < 3; i++)
        {
            int randN = dis(gen);
            position.push_back(randN);
        }
#endif
#ifdef UNIQUEPOSITION
        while (1) {
            position.clear();
            for (int i = 0; i < 3; i++)
            {
                int randN = dis(gen);
                position.push_back(randN);
            }
            if (std::find(massPositions.begin(), massPositions.end(), position) == massPositions.end()) {
                break;
            }
        }
        massPositions.push_back(position);
#endif
        return position;
    }

    void generateGenes() {
        for (int i = 0; i < populationSize; i++) {
            massPositions.clear();
            std::vector<GENE> cubotGene;
            for (int i = 0; i < (pow(3, depth) - 1) / 2; i++) {
                int level = (int)(log(i * 2 - 1) / log(3));
                bool existence = false;
                std::random_device rd;
                std::mt19937 gen(rd());
                std::uniform_int_distribution<> dis(1, 20);
                if (level < dis(gen)) {
                    existence = true;
                }
                cubotGene.push_back({ existence, randomSpringFactor(), generatePosition(location) });
            }
            populationGene.push_back(cubotGene);
        }
    }

    void generatecubots() {
        std::cout << x;
        for (int i = 0; i < populationSize; i++) {
            double X = 0;
            double Y = 0;
            cubots.push_back(cubot(X, Y, 10.0, populationGene[i]));
        }
    }

    void simUpdate() {
#ifndef GRAPHICS
#pragma omp parallel for num_threads(NUM_THREADS)
#endif
        for (int i = 0; i < populationSize; i++) {
            cubots[i].cubotUpdate();
        }
    }

    void simDraw() {
        for (int i = 0; i < populationSize; i++) {
            cubots[i].cubotDraw();
        }
    }

    void calculatePopulationDistance() {
        populationDistance.clear();
        populationDistance.shrink_to_fit();
        for (int i = 0; i < populationSize; i++) {
            double x = 0; double y = 0;
            for (int j = 0; j < cubots[i].cubemasses.size(); j++) {
                x = x + cubots[i].cubemasses[j].pos[0];
                y = y + cubots[i].cubemasses[j].pos[1];
            }
            x = x / cubots[i].cubemasses.size();
            y = y / cubots[i].cubemasses.size();
            double distance[2] = { fabs(x - cubots[i].initialLocation[0]), fabs(y - cubots[i].initialLocation[1]) };
            double distanceNorm = calcmag(distance, 2);
            populationDistance.push_back(distanceNorm);
        }
        averageDistance = 0;
        maxDistance = 0;
        int counter = 0;
        for (int i = 0; i < populationSize; i++) {
            averageDistance = averageDistance + populationDistance[i];
            maxDistance = std::max(maxDistance, populationDistance[i]);
            popDis << populationDistance[i] << " ";
        }
        popDis << "\n"; popDis.flush();
        averageDistance = averageDistance / populationSize;
        std::cout << "Maximum Distance: " << maxDistance << std::endl;
        std::cout << "Average Distance: " << averageDistance << std::endl;
    }

};

#ifdef GRAPHICS
Simulation sim1(cubotNumber);
#endif

void drawGrid() {

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glRotatef(0, 1.0, 0.0, 2.0);
    glRotatef(0, 0.0, 1.0, 2.0);
    glColor3f(0.2, 0.2, 0.2);
    glBegin(GL_LINES);
    for (GLfloat i = -30; i < 30; i += 1)
    {
        glVertex3f(i, 0, 30);
        glVertex3f(i, 0, -30);
        glVertex3f(30, 0, i);
        glVertex3f(-30, 0, i);
    }
    glEnd();

}

void display()
{
    const double len = 2.0; 
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glLoadIdentity();
   
    double Ex = -2 * dim * Sin(th) * Cos(ph);
    double Ey = +2 * dim * Sin(ph);
    double Ez = +2 * dim * Cos(th) * Cos(ph);
    gluLookAt(Ex, Ey, Ez, 0, 0, 0, 0, Cos(ph), 0);
   
   

    drawGrid();



#ifdef GRAPHICS
    glColor3f(1, 1, 1);
    sim1.startSim(simulationTime);
    T = T + timeStep;
#endif



    Frames++;
    GLint t = glutGet(GLUT_ELAPSED_TIME);
    if (t - T0 >= 1000) {
        GLfloat seconds = (t - T0) / 1000.0;
        fps = Frames / seconds;
        T0 = t;
        Frames = 0;
    }

    glRasterPos3d(0.0, 2, 0.0);
  
    glColor3f(1, 1, 1);
    
    glFlush();
   
    glutSwapBuffers();
}
void special(int key, int x, int y)
{
    if (key == GLUT_KEY_RIGHT)
        th += 10;
    else if (key == GLUT_KEY_LEFT)
        th -= 10;
    else if (key == GLUT_KEY_UP)
    {
        if (ph + 5 < 90)
        {
            ph += 5;
        }
    }
    else if (key == GLUT_KEY_DOWN)
    {
        if (ph - 5 > 0)
        {
            ph -= 5;
        }
    }
    th %= 360;
    ph %= 360;
    glutPostRedisplay();
}

void Project(double fov, double asp, double dim){

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if (fov)
        gluPerspective(fov, asp, dim / 16, 16 * dim);
    else
        glOrtho(-asp * dim, asp * dim, -dim, +dim, -dim, +dim);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void reshape(int width, int height){
    asp = (height > 0) ? (double)width / height : 1;
    glViewport(0, 0, width, height);
    Project(fov, asp, dim);
}

void idle()
{
    glutPostRedisplay();
}

int main(int argc, char* argv[]) {
#ifdef GRAPHICS
    glutInit(&argc, argv);
    glutInitWindowSize(1000, 800);
    glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
    glutCreateWindow("jb4512");
    glutIdleFunc(idle);
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutSpecialFunc(special);
    glutMainLoop();
#endif

#ifndef GRAPHICS
    Simulation sim1(cubotNumber);
    while (1) {
        double start_time = omp_get_wtime();
        sim1.startSim(simulationTime);
        T = T + timeStep;
        //printf("%f\n", T);
        double time = omp_get_wtime() - start_time;
    }
#endif

    return 0;
}
