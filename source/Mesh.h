#pragma once
#include <iostream>

#define PI 3.14159265359


enum
{
    OBSTACLE = -2,
    SILID = -1,
    COMP,
    RIGHT,
    TOP,
    LEFT,
    BOTTOM,
    RIGHT_TOP,
    LEFT_TOP,
    LEFT_BOTTOM,
    RIGHT_BOTTOM
};

class Mesh
{
protected:
    int NX;
    int NY;
    int *solid;
    double *x, *y;

    size_t scalar_index(int i, int j)
    {
        return (NX * j + i);
    }

public:
    double *xNorm, *yNorm;

    Mesh();
    Mesh(int NX_, int NY_);
    ~Mesh();
    void InitMesh();
};

template<typename T>
inline  T p2(const T x) {
	return(x * x);
}