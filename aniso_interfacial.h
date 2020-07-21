#ifndef ANISO_INTERFACIAL_H
#define ANISO_INTERFACIAL_H

#include <math.h>
double calc_f2x(double phix, double phiy, double phiz, double dphi[3][3], double ddphi[3][3][3]);
double calc_f2y(double phix, double phiy, double phiz, double dphi[3][3], double ddphi[3][3][3]);
double calc_f4x(double phix, double phiy, double phiz, double dphi[3][3], double ddphi[3][3][3]);
double calc_f4y(double phix, double phiy, double phiz, double dphi[3][3], double ddphi[3][3][3]);
double calc_f2z(double phix, double phiy, double phiz, double dphi[3][3], double ddphi[3][3][3]);
double calc_f4z(double phix, double phiy, double phiz, double dphi[3][3], double ddphi[3][3][3]);
#endif
