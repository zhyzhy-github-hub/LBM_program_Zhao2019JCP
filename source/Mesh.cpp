#include "Mesh.h"
Mesh::~Mesh()
{
    delete[] x;
    delete[] y;
    delete[] solid;
    delete[] xNorm;
    delete[] yNorm;
}

Mesh::Mesh(int NX_, int NY_)
{
    NX = NX_;
    NY = NY_;
    const unsigned int NN = NX * NY;
    x = new double[NN];
    y = new double[NN];
    xNorm = new double [NN];
    yNorm = new double [NN];
    solid = new int[NN];
}

void Mesh::InitMesh()
{
    const double LC = NY - 1;
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t index = scalar_index(i,j);
            x[index] = i;
            y[index] = j;
            solid[index] = COMP;
            xNorm[index] = double(i) / LC;
            yNorm[index] = double(j) / LC;
        }
    }
}