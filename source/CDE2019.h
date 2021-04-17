#pragma once

#include "Mesh.h"

#include <cmath>
#include <iostream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>

#define DIFF_SCALING

#define w0 4.0 / 9
#define wc 1.0 / 9
#define wd 1.0 / 36

#define delta0 0.0001
#define C0 0.01

const double w[9] = {w0, wc, wc, wc, wc, wd, wd, wd, wd};

const int dirx[] = {0,1,0,-1, 0,1,-1,-1, 1};
const int diry[] = {0,0,1, 0,-1,1, 1,-1,-1};

struct BADE
{
    double Bx;
    double By;
};

class CDE2019
{
private:
    const int Q9 = 9;
    int NX;
    int NY;

    double dx;

#ifdef DIFF_SCALING
    double dt;
    double mu;
    double c;
    double Cs;
#elif
    double dt, c, Cs, Cs2;
    double D;
    double tau;

#endif

    double *s;

    double *C;
    double *C_old;
    double *F;
    double *f1; // 1 is post collision
    double *f0; // 0 is post streaming

    size_t scalar_index(int i, int j)
    {
        return (NX * j + i);
    }
    size_t field_index(unsigned int x, unsigned int y, unsigned int d)
    {
        return (NX * (NY * (d) + y) + x);
        // return (Q9 * (NX * y + x) + d);
    }

    void getMR(double *m, double *R);
    void getFi(double *Fi, double F);
    double meqADE2019(int k, double phi, double D, BADE B);
    double feqADE2019(int k, double phi, double D, BADE B);

public:
    CDE2019(int NX_, int NY_);
    ~CDE2019();

#ifdef DIFF_SCALING
    void initCDE2019_Diff(double mu, double D_, double s_add);
#elif
    void initCDE2019(double c_, double D_, double s_add);
#endif

    void getCOld();

    void initScalar(size_t n, double *x, double *y);

    void initFeqAndMeq();

    void AnsScalar(size_t n, double *x, double *y);

    void outputTec(double t, double *x, double *y);

    void getF(size_t t_, double *x, double *y);

    void collision();
    void streaming();

    void streamingBGK();
    void collisionBGK();

    void m2f();
    double errC();

    double time(size_t n)
    {
        return n * dt;
    }
};
