#ifndef EWALD_H
#define EWALD_H
#include<vector>

using namespace std;

double norm(double x,double y,double z);

double dp(double x,double y,double z,double xx,double yy,double zz);

double RealandReciprocalSpace(vector<double> const &r, double Lx, double Ly, double Lz, double kappa, int ds, int size, int rank);

double PointEnergy(vector<double> const &r, int size, int rank);

double Dipole(vector<double> const &r);

#endif
