#define _USE_MATH_DEFINES
#include<iostream>
#include<vector>
#include<random>
#include<cmath>
#include<ctime>
#include "fun_final.h"
#include<algorithm>

using namespace  std;

int main()
{
    clock_t t_begin = clock();
    // Input
    // nmove moves per cycle, nycle cycles, total of nmove*ncycle MC steps
    int ncycl = 1000, nmove = 10000;
    int npav = 250, nexc = 1;
    int npavexc = npav + nexc;
    int icycl, imove;
    
    
    int npart = 250;
    double L = 100;                                       // Angstrom
    double T = 300;                                             // Kelvin
    double P = 1e6;                                               // Pascal

    // Constant
    double h_c = 6.626e-34;                                      // (m2.kg/s)
    double k_B = 1.381e-23;                                      // J/K
    double N_A = 6.022e23;
    double m = 0.03994 / N_A;                                   // kg/molecule

    // Parameters
    double lambda = h_c / pow(2 * M_PI * m * k_B * T, 0.5);     // m
    double eps = k_B * 119.8;                               // J
    double beta = 1 / (k_B * T);                               // 1/J
    double lambda3 = lambda * lambda * lambda;
    double mu0 = (k_B * T) * log(lambda3);                 // J
    double mu = ((beta * mu0) + log(beta * P)) / beta;     // J
    double V = pow(L, 3) * pow(10, -30);      // m3
    double sig = 3.404;                                    // Angstrom
    double sig6 = pow(sig, 6.0);
    double sig12 = pow(sig, 12.0);
    double rc2 = 900;                                     // Angstrom2
    double A = 4 * eps * sig12;
    double B = 4 * eps * sig6;
    double E_trunc_corr = (A / pow(rc2, 6.0)) - (B / pow(rc2, 3.0));
    double zz = exp(beta * mu) / lambda3;
    double zzV = zz * V;



    // Position of Particles
    vector<vector<double>> x(npart, vector<double>(3));
    vector<vector<double>> r2(npart, vector<double>(npart, 0));
    vector<double> en(npart);

    // Variable declarations
    int i, j;
    vector<vector<double>> E_pair(npart, vector<double>(npart, 0));
    double U = 0, enn = 0, delta = 25.0, dU = 0;
    int ipart;
    float ndisp = 0, ndispaccpt = 0;
    vector<double> ipart_x(3);
    double arg;

    srand(time(nullptr));
    for (i = 0; i < npart; i++)
    {
        for (j = 0; j < 3; j++)
        {
            x[i][j] = MyRand(0, L);
        }
    }

    // Minimum Image Convention

    for (i = 0; i < npart; i++)
    {
        r2[i][i] = 0;
        for (j = i + 1; j < npart; j++)
        {
            r2[i][j] = Findr2(x[i], x[j], L);
            r2[j][i] = r2[i][j];
        }

    }

    // Energy of system

    for (i = 0; i < npart; i++)
    {
        for (j = i + 1; j < npart; j++)
        {
            E_pair[i][j] = FindEPair(r2[i][j], rc2, A, B, E_trunc_corr);
            E_pair[j][i] = E_pair[i][j];
        }
        E_pair[i][i] = 0;
        en[i] = 0;
        for (j = 0; j < npart; j++)
            en[i] += E_pair[i][j];

    }

    for (i = 0; i < npart; i++)
    {
        U += en[i];
        U = U / 2.0;
    }

    // GCMC Simulation

    for (icycl = 0; icycl < ncycl; icycl++)
    {
        for (imove = 0; imove < nmove; imove++)
        {
            // displace particle
            if (MyRand(0, npavexc) < npav)
            {
                ndisp += 1;
                ipart = rand() % npart;

                // Periodic Condition
                for (j = 0; j < 3; j++)
                    ipart_x[j] = ApplyPB(x[ipart][j] + delta * MyRand(-0.5, 0.5), L);

                enn = 0;

                for (i = 0; i < ipart; i++)
                    enn += FindEPair(Findr2(ipart_x, x[i], L), rc2, A, B, E_trunc_corr);

                for (i = ipart + 1; i < npart; i++)
                    enn += FindEPair(Findr2(ipart_x, x[i], L), rc2, A, B, E_trunc_corr);

                dU = enn - en[ipart];

                if (MyRand(0, 1) < exp(-beta * dU))
                {
                    x[ipart][0] = ipart_x[0];
                    x[ipart][1] = ipart_x[1];
                    x[ipart][2] = ipart_x[2];
                    en[ipart] = enn;
                    U += dU;
                    ndispaccpt += 1;
                }
            }

            // exchange particle
            else
            {
                // add particle
                if (MyRand(0, 1) < 0.5)
                {
                    for (j = 0; j < 3; j++)
                        ipart_x[j] = MyRand(0, L);

                    enn = 0;
                    for (i = 0; i < npart; i++)
                        enn += FindEPair(Findr2(ipart_x, x[i], L), rc2, A, B, E_trunc_corr);

                    arg = zzV * exp(-beta * enn) / (npart + 1);

                    if (MyRand(0, 1) < arg)
                    {
                        npart += 1;
                        x.push_back({ ipart_x[0], ipart_x[1] , ipart_x[2] });
                        en.push_back(enn);
                        U += enn;
                    }
                }
                // delete particle
                else
                {
                    ipart = rand() % npart;
                    arg = npart / zzV * exp(beta * en[ipart]);
                    if (MyRand(0, 1) < arg)
                    {
                        npart -= 1;
                        U -= en[ipart];
                        x.erase(x.begin() + ipart);
                        en.erase(en.begin() + ipart);
                    }
                }

            }
        }

        delta = ReviseDelta(ndispaccpt / ndisp, delta);
        ndispaccpt = 0;
        ndisp = 0;
        //Log();
        Display(icycl, npart, delta, U, (clock() - t_begin) / (float) CLOCKS_PER_SEC);

    }
}
