#include<iostream>
#include<vector>
#include<random>
#include<cmath>
#include<ctime>
#include<chrono>
using namespace  std;


/*double my_pow(double x, size_t n) {
    double r = 1.0;

    while(n > 0){
        r *= x;
        --n;
    }

    return r;
}*/

double MyRand(double a,double b)
{
   return ((double) rand() / (RAND_MAX + 1.0) * (b - a) + a);
  
    
}
double Findr2(vector <double> pos1, vector <double> pos2, double L)
{
    double dx = fabs(pos1[0] - pos2[0]);
    double dy = fabs(pos1[1] - pos2[1]);
    double dz = fabs(pos1[2] - pos2[2]);
    dx = dx - L * (floor(dx / L));
    dy = dy - L * (floor(dy / L));
    dz = dz - L * (floor(dz / L));
    double r2 = dx*dx + dy*dy + dz*dz;
    return r2;
}

double FindEPair(double r2, double r_c2, double A, double B, double E_trunc_corr)
{
    double E_pair;
    if (r2 <= r_c2)
    {
     //value = (A/my_pow(r,6)) - (B/my_pow(r,3)) - E_trunc_corr;
       E_pair = (A / pow(r2, 6.0)) - (B / pow(r2, 3.0));
    }
    else{
        E_pair = 0;
    }
    return E_pair;
}

double ApplyPB(double coordinate, double L)
{
    if(coordinate > L)
    {
        coordinate -= L;
    }
    else if(coordinate < 0)
    {
        coordinate += L;
    }
    return  coordinate;
}
double ReviseDelta(float ratio, double delta)
{
    if (ratio > 0.55)
    delta = delta * 1.05;
    else if (ratio < 0.45)
    delta = delta * 0.95;
    
    return delta;
}

void Display(int icycl, int npart, double delta, double U, float t_elapsed)
{
    cout << "icycl =" << icycl << ", npart =" << npart << ", delta =" << delta << ", U =" << U << ", t =" << t_elapsed;
    cout << '\n';
}
