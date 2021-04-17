#include "source/CDE2019.h"
#include "source/Mesh.h"

int main()
{
    const int NX = 41;
    const int NY = 41;
    Mesh mesh0(NX, NY);
    mesh0.InitMesh();

    CDE2019 CDE0(NX, NY);

    double snu = 1;
    double nu = 0.1;
    double mu = (1.0 / snu - 0.5) / (3 * nu);

    CDE0.initCDE2019_Diff(mu, nu, 0);
    CDE0.initScalar(0, mesh0.xNorm, mesh0.yNorm);
    CDE0.getCOld();
    CDE0.initFeqAndMeq();

    CDE0.outputTec(0, mesh0.xNorm, mesh0.yNorm);

    for (int n = 1; n < 6000; ++n)
    {
        CDE0.AnsScalar(n, mesh0.xNorm, mesh0.yNorm);
        CDE0.getF(n, mesh0.xNorm, mesh0.yNorm);
        CDE0.m2f();
        CDE0.streaming();
        CDE0.collision();


        // CDE0.collisionBGK();
        // CDE0.streamingBGK();


        if (fabs(CDE0.time(n) - 0.25) < 1e-7)
        {
            CDE0.outputTec(n, mesh0.xNorm, mesh0.yNorm);
            break;
        }

        // if (n % 500 == 0)
        // {
        //     double time = CDE0.time(n);
        //     std::cout << " n = " << n << ", time = " << time << std::endl;
        //     CDE0.outputTec(n, mesh0.xNorm, mesh0.yNorm);
        // }
    }
    double Er = CDE0.errC();
    std::cout << "E_r = " << Er << std::endl;
}