#include "CDE2019.h"

CDE2019::CDE2019(int NX_, int NY_)
{

    NX = NX_;
    NY = NY_;
    const unsigned int NN = NX * NY;
    C = new double[NN];
    C_old = new double[NN];
    F = new double[NN];

    f1 = new double[NN * Q9];
    f0 = new double[NN * Q9];

    s = new double[Q9];
}

CDE2019::~CDE2019()
{
    delete[] C;
    delete[] C_old;
    delete[] F;

    delete[] f0;
    delete[] f1;

    delete[] s;
}

#ifdef DIFF_SCALING
void CDE2019::initCDE2019_Diff(double mu_, double D_, double s_add)
{
    mu = mu_;
    dx = 1.0 / (NY - 1);
    dt = mu * dx * dx;
    // dt = mu * dx ;
    c = dx / dt;
    Cs = c / sqrt(3);

    double s3 = 2.0 / (6 * mu * D_ + 1);
    for (int i = 0; i < Q9; ++i)
    {
        if (i == 3 || i == 5)
        {
            s[i] = s3;
        }
        else
        {
            s[i] = s3; //1.0 + s_add;
        }
    }

    using std::cout;
    using std::endl;
    cout << "dx = " << dx << endl;
    cout << "dt = " << dt << endl;
    cout << "mu = " << mu << endl;
    cout << " c = " << c << endl;
    for (int i = 0; i < Q9; ++i)
    {
        cout << "s[" << i << "] " << s[i] << " | ";
    }
    cout << endl;
}
#elif
void CDE2019::initCDE2019(double c_, double D_, double s_add)
{
    c = c_;
    dt = dx / c;
    Cs = c / sqrt(3);
    Cs2 = c * c / 3;

    D = D_;
    tau = 0.5 + D / (Cs2 * dt);
    double tau_inv = 1.0 / tau;

    for (int i = 0; i < Q9; ++i)
    {
        if (i == 3 || i == 5)
        {
            s[i] = tau_inv;
        }
        else
        {
            s[i] = 1.0 + s_add;
        }
    }
}
#endif

void CDE2019::initScalar(size_t n, double *x, double *y)
{
    double t = dt * n;
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t index = scalar_index(i, j);
            double y_ = y[index];
            double x_ = x[index];

            // C[index] = C0 * exp(-(p2(x_ - 0.5) + p2(y_ - 0.5)) / (2 * p2(delta0)));

            C[index] = (t + 1) * sin(2 * PI * x_) * cos(2 * PI * y_);
        }
    }
}

void CDE2019::getCOld()
{
    size_t NN = NX * NY;
    for (int i = 0; i < NN; ++i)
    {
        C_old[i] = C[i];
    }
}

void CDE2019::getFi(double *Fi, double F)
{
    for (int i = 0; i < Q9; ++i)
    {
        if (i == 0)
        {
            Fi[i] = dt * w0 * F;
        }
        // else if (i < 5 && i > 0)
        else if (i == 1 || i == 2 || i == 3 || i == 4)
        {
            Fi[i] = dt * wc * F;
        }
        else
        {
            Fi[i] = dt * wd * F;
        }
    }
}

double CDE2019::meqADE2019(int k, double phi, double D, BADE B)
{
    switch (k)
    {
    case 0:
        return phi;
    case 1:
        return 2 * D - 4 * phi;
    case 2:
        return 3 * phi - 2 * D;
    case 3:
        return mu * dx * B.Bx;
    case 4:
        return -mu * dx * B.Bx;
    case 5:
        return mu * dx * B.By;
    case 6:
        return -mu * dx * B.By;
    case 7:
        return 0;
    case 8:
        return 0;
    }
}

// double CDE2019::meqADE2019(int k, double phi, double D, BADE B)
// {
//     double feq[9];
//     double meq[9];

//     feq[0] = w0 * (2 * phi - D);
//     feq[1] = wc * (2 * phi - D + 3 * B.Bx / c + (3 * (D - phi) * 1) / 2);
//     feq[2] = wc * (2 * phi - D + 3 * B.By / c + (3 * (D - phi) * 1) / 2);
//     feq[3] = wc * (2 * phi - D - 3 * B.Bx / c + (3 * (D - phi) * 1) / 2);
//     feq[4] = wc * (2 * phi - D - 3 * B.By / c + (3 * (D - phi) * 1) / 2);
//     feq[5] = wd * (2 * phi - D + (3 * (B.Bx + B.By)) / c + (3 * (D - phi) * 2) / 2);
//     feq[6] = wd * (2 * phi - D + (3 * (-B.Bx + B.By)) / c + (3 * (D - phi) * 2) / 2);
//     feq[7] = wd * (2 * phi - D - (3 * (B.Bx + B.By )) / c + (3 * (D - phi) * 2) / 2);
//     feq[8] = wd * (2 * phi - D - (3 * (-B.Bx + B.By)) / c + (3 * (D - phi) * 2) / 2);

//     getMR(meq, feq);
//     return meq[k];
// }

double CDE2019::feqADE2019(int k, double phi, double D, BADE B)
{
    switch (k)
    {
    case 0:
        return (w0 * (-D + 2 * phi));
    case 1:
        return (wc * (6 * B.Bx + c * (D + phi)) / (2 * c));
    case 2:
        return (wc * (6 * B.By + c * (D + phi)) / (2 * c));
    case 3:
        return (wc * (-6 * B.Bx + c * (D + phi)) / (2 * c));
    case 4:
        return (wc * (-6 * B.By + c * (D + phi)) / (2 * c));
    case 5:
        return (wd * (3 * B.Bx + 3 * B.By + c * (2 * D - phi)) / c);
    case 6:
        return (wd * (-3 * B.Bx + 3 * B.By + c * (2 * D - phi)) / c);
    case 7:
        return (-wd * (3 * B.Bx + 3 * B.By - c * (2 * D - phi)) / c);
    case 8:
        return (wd * (3 * B.Bx - 3 * B.By + c * (2 * D - phi)) / c);
        // case 1:
        //     return wc * (2 * phi - D + 3 * B.Bx / c + (3 * (D - phi) * 1) / 2);
        // case 2:
        //     return wc * (2 * phi - D + 3 * B.By / c + (3 * (D - phi) * 1) / 2);
        // case 3:
        //     return wc * (2 * phi - D - 3 * B.Bx / c + (3 * (D - phi) * 1) / 2);
        // case 4:
        //     return wc * (2 * phi - D - 3 * B.By / c + (3 * (D - phi) * 1) / 2);
        // case 5:
        //     return wd * (2 * phi - D + (3 * (B.Bx + B.By)) / c + (3 * (D - phi) * 2) / 2);
        // case 6:
        //     return wd * (2 * phi - D + (3 * (-B.Bx + B.By)) / c + (3 * (D - phi) * 2) / 2);
        // case 7:
        //     return wd * (2 * phi - D - (3 * (B.Bx + B.By)) / c + (3 * (D - phi) * 2) / 2);
        // case 8:
        //     return wd * (2 * phi - D - (3 * (-B.Bx + B.By)) / c + (3 * (D - phi) * 2) / 2);
    }
}

void CDE2019::AnsScalar(size_t n, double *x, double *y)
{
    // initScalar(t, x, y);
    double t = dt * n;
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t index = scalar_index(i, j);
            double y_ = y[index];
            double x_ = x[index];
            C_old[index] = (t + 1) * sin(2 * PI * x_) * cos(2 * PI * y_);

            // double deltaD = sqrt(0.1 * t * 2);
            // C_old[index] = p2(delta0) / (p2(delta0) + p2(deltaD)) * C0 * exp(-(p2(x_ - 0.5) + p2(y_ - 0.5)) / (2 * (p2(delta0) + p2(deltaD))));
        }
    }
}

void CDE2019::initFeqAndMeq()
{
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t s_idx = scalar_index(i, j);
            double phi = C[s_idx];
            double D_ = sin(phi);
            BADE B_ = {phi, phi};

            for (int k = 0; k < Q9; ++k)
            {
                size_t f_idx = field_index(i, j, k);
                f0[f_idx] = feqADE2019(k, phi, D_, B_);
                f1[f_idx] = meqADE2019(k, phi, D_, B_);
                // f1[f_idx] = feqADE2019(k, phi, D_, B_);
            }
        }
    }
}

void CDE2019::getF(size_t t_, double *x, double *y)
{
    double t = t_ * dt;
    // std::cout << t << std::endl;
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t index = scalar_index(i, j);
            double y_ = y[index];
            double x_ = x[index];
            F[index] = sin(2 * PI * x_) * cos(2 * PI * y_)                                                                                           //
                       + 2 * PI * (t + 1) * cos(2 * PI * x_ + 2 * PI * y_)                                                                           //
                       + 0.4 * p2(PI) * p2(t + 1) * sin((t + 1) * sin(2 * PI * x_) * cos(2 * PI * y_)) * p2(cos(2 * PI * x_)) * p2(cos(2 * PI * y_)) //
                       + 0.4 * p2(PI) * p2(t + 1) * sin((t + 1) * sin(2 * PI * x_) * cos(2 * PI * y_)) * p2(sin(2 * PI * x_)) * p2(sin(2 * PI * y_)) //
                       + 0.8 * p2(PI) * (t + 1) * cos((t + 1) * sin(2 * PI * x_) * cos(2 * PI * y_)) * (sin(2 * PI * x_)) * (cos(2 * PI * y_));      //
        }
    }
}

void CDE2019::getMR(double *m, double *R)
{
    const double f0 = R[0];
    const double f1 = R[1];
    const double f2 = R[2];
    const double f3 = R[3];
    const double f4 = R[4];
    const double f5 = R[5];
    const double f6 = R[6];
    const double f7 = R[7];
    const double f8 = R[8];
    m[0] = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
    m[1] = (-4 * f0 - f1 - f2 - f3 - f4 + 2 * f5 + 2 * f6 + 2 * f7 + 2 * f8);
    m[2] = (4 * f0 - 2 * f1 - 2 * f2 - 2 * f3 - 2 * f4 + f5 + f6 + f7 + f8);
    m[3] = (f1 - f3 + f5 - f6 - f7 + f8);
    m[4] = (-2 * f1 + 2 * f3 + f5 - f6 - f7 + f8);
    m[5] = (f2 - f4 + f5 + f6 - f7 - f8);
    m[6] = (-2 * f2 + 2 * f4 + f5 + f6 - f7 - f8);
    m[7] = (f1 - f2 + f3 - f4);
    m[8] = (f5 - f6 + f7 - f8);
}

void CDE2019::m2f()
{
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            double F_ = F[scalar_index(i, j)];

            double m0 = f1[field_index(i, j, 0)];
            double m1 = f1[field_index(i, j, 1)];
            double m2 = f1[field_index(i, j, 2)];
            double m3 = f1[field_index(i, j, 3)];
            double m4 = f1[field_index(i, j, 4)];
            double m5 = f1[field_index(i, j, 5)];
            double m6 = f1[field_index(i, j, 6)];
            double m7 = f1[field_index(i, j, 7)];
            double m8 = f1[field_index(i, j, 8)];

            f1[field_index(i, j, 0)] = m0 / 9 - m1 / 9 + m2 / 9;
            f1[field_index(i, j, 1)] = m0 / 9 - m1 / 36 - m2 / 18 + m3 / 6 - m4 / 6 + m7 / 4;
            f1[field_index(i, j, 2)] = m0 / 9 - m1 / 36 - m2 / 18 + m5 / 6 - m6 / 6 - m7 / 4;
            f1[field_index(i, j, 3)] = m0 / 9 - m1 / 36 - m2 / 18 - m3 / 6 + m4 / 6 + m7 / 4;
            f1[field_index(i, j, 4)] = m0 / 9 - m1 / 36 - m2 / 18 - m5 / 6 + m6 / 6 - m7 / 4;
            f1[field_index(i, j, 5)] = m0 / 9 + m1 / 18 + m2 / 36 + m3 / 6 + m4 / 12 + m5 / 6 + m6 / 12 + m8 / 4;
            f1[field_index(i, j, 6)] = m0 / 9 + m1 / 18 + m2 / 36 - m3 / 6 - m4 / 12 + m5 / 6 + m6 / 12 - m8 / 4;
            f1[field_index(i, j, 7)] = m0 / 9 + m1 / 18 + m2 / 36 - m3 / 6 - m4 / 12 - m5 / 6 - m6 / 12 + m8 / 4;
            f1[field_index(i, j, 8)] = m0 / 9 + m1 / 18 + m2 / 36 + m3 / 6 + m4 / 12 - m5 / 6 - m6 / 12 - m8 / 4;
            // f1[field_index(i, j, 0)] += +dt * w0 * F_;
            // f1[field_index(i, j, 1)] += +dt * wc * F_;
            // f1[field_index(i, j, 2)] += +dt * wc * F_;
            // f1[field_index(i, j, 3)] += +dt * wc * F_;
            // f1[field_index(i, j, 4)] += +dt * wc * F_;
            // f1[field_index(i, j, 5)] += +dt * wd * F_;
            // f1[field_index(i, j, 6)] += +dt * wd * F_;
            // f1[field_index(i, j, 7)] += +dt * wd * F_;
            // f1[field_index(i, j, 8)] += +dt * wd * F_;
        }
    }
}

void CDE2019::streaming()
{
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            unsigned int ip1 = (i + 1) % (NX-1);
            unsigned int jp1 = (j + 1) % (NY-1);
            unsigned int im1 = (NX -1 + i - 1) % (NX-1);
            unsigned int jm1 = (NY -1 + j - 1) % (NY-1);

            double ft0 = f1[field_index(i, j, 0)]    ;
            double ft1 = f1[field_index(im1, j, 1)]  ;
            double ft2 = f1[field_index(i, jm1, 2)]  ;
            double ft3 = f1[field_index(ip1, j, 3)]  ;
            double ft4 = f1[field_index(i, jp1, 4)]  ;
            double ft5 = f1[field_index(im1, jm1, 5)];
            double ft6 = f1[field_index(ip1, jm1, 6)];
            double ft7 = f1[field_index(ip1, jp1, 7)];
            double ft8 = f1[field_index(im1, jp1, 8)];

            C[scalar_index(i, j)] = ft0                     //
                                    + ft1 + ft2 + ft3 + ft4 //
                                    + ft5 + ft6 + ft7 + ft8;

            f0[field_index(i, j, 0)] = ft0;
            f0[field_index(i, j, 1)] = ft1;
            f0[field_index(i, j, 2)] = ft2;
            f0[field_index(i, j, 3)] = ft3;
            f0[field_index(i, j, 4)] = ft4;
            f0[field_index(i, j, 5)] = ft5;
            f0[field_index(i, j, 6)] = ft6;
            f0[field_index(i, j, 7)] = ft7;
            f0[field_index(i, j, 8)] = ft8;
        }
    }
}

void CDE2019::collision()
{
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t index = scalar_index(i, j);

            double Fi[9];
            getFi(Fi, F[index]);

            double MF[9];
            getMR(MF, Fi);

            // double phi = C_old[index];
            double phi = C[index];
            // double D_ = sin(phi);
            double D_ = sin(phi);
            BADE B_; // = {phi, phi};
            B_.By = phi;
            B_.Bx = phi;

            double f0_[9];
            double m0[9];
            for (int k = 0; k < Q9; ++k)
            {
                f0_[k] = f0[field_index(i, j, k)];
            }
            getMR(m0, f0_);

            for (int k = 0; k < Q9; ++k)
            {
                size_t fidx = field_index(i, j, k);
                f1[fidx] = m0[k] - s[k] * (m0[k] - meqADE2019(k, phi, D_, B_)) + MF[k];
            }
        }
    }
}

void CDE2019::collisionBGK()
{
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t index = scalar_index(i, j);
            double phi = C[index];
            // double phi = C_old[index];
            double D_ = sin(phi);
            // double D_ = sin(phi);
            BADE B_ = {phi, phi};

            for (int k = 0; k < Q9; ++k)
            {
                size_t f_idx = field_index(i, j, k);
                // f1[f_idx] = f0[f_idx] - s[3] * (f0[f_idx] - feqADE2019(k, phi, D_, B_)) + w[k] * dt * F[index];
                f1[f_idx] = (1 - s[3]) * f0[f_idx] + s[3] * feqADE2019(k, phi, D_, B_) + w[k] * dt * F[index];
                // f1[f_idx] = feqADE2019(k, phi, D_, B_)  + w[k] * dt * F[index];

            }
        }
    }
}

void CDE2019::streamingBGK()
{
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {

            // double C_ = 0;
            // for (int k = 0; k < Q9; ++k)
            // {
            //     unsigned int xmd = (NX + i - dirx[k]) % NX;
            //     unsigned int ymd = (NY + j - diry[k]) % NY;

            //     f0[field_index(i, j, k)] = f1[field_index(xmd, ymd, k)];

            //     C_ += f0[field_index(i, j, k)];
            // }
            // C[scalar_index(i, j)] = C_;

            unsigned int ip1 = (i + 1) % (NX-1);
            unsigned int jp1 = (j + 1) % (NY-1);
            unsigned int im1 = (NX - 1 + i - 1) % (NX-1);
            unsigned int jm1 = (NY - 1 + j - 1) % (NY-1);

            double ft0 = f1[field_index(i, j, 0)];
            double ft1 = f1[field_index(im1, j, 1)];
            double ft2 = f1[field_index(i, jm1, 2)];
            double ft3 = f1[field_index(ip1, j, 3)];
            double ft4 = f1[field_index(i, jp1, 4)];
            double ft5 = f1[field_index(im1, jm1, 5)];
            double ft6 = f1[field_index(ip1, jm1, 6)];
            double ft7 = f1[field_index(ip1, jp1, 7)];
            double ft8 = f1[field_index(im1, jp1, 8)];

            C[scalar_index(i, j)] = ft0 + ft1 + ft2 + ft3 + ft4 + ft5 + ft6 + ft7 + ft8;
            // C[scalar_index(i, j)] = C_old[scalar_index(i, j)];

            f0[field_index(i, j, 0)] = ft0;
            f0[field_index(i, j, 1)] = ft1;
            f0[field_index(i, j, 2)] = ft2;
            f0[field_index(i, j, 3)] = ft3;
            f0[field_index(i, j, 4)] = ft4;
            f0[field_index(i, j, 5)] = ft5;
            f0[field_index(i, j, 6)] = ft6;
            f0[field_index(i, j, 7)] = ft7;
            f0[field_index(i, j, 8)] = ft8;
        }
    }
}

double CDE2019::errC()
{
    double temp1 = 0;
    double temp2 = 0;
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t idx = scalar_index(i, j);
            temp1 += p2(C[idx] - C_old[idx]);
            temp2 += p2(C_old[idx]);
        }
    }
    double err = sqrt(temp1) / sqrt(temp2);
    return err;
}

void CDE2019::outputTec(double t, double *x, double *y)
{
    std::ostringstream name;
    name << "t_" << t << "_"
         << ".dat";
    std::ofstream out(name.str().c_str());
    out << "Title= \"Lid Driven Flow\"\n";
    out << "VARIABLES = \"X\", \"Y\", \"C\",  \"Co\", \n";
    out << "ZONE T= \"BOX\",I=" << NX << ",J=" << NY << ",F=	POINT" << std::endl;

    for (size_t j = 0; j < NY; ++j)
    {
        for (size_t i = 0; i < NX; ++i)
        {
            const size_t index = scalar_index(i, j);
            const size_t findex = field_index(i, j, 0);
            out << std::fixed << std::setprecision(10)
                << " " << x[index] << " " << y[index]
                // << " " << f0[findex]
                // << " " << f1[findex]
                << " " << C[index]
                << " " << C_old[index]
                << std::endl;
        }
    }
}
