#include "aniso_interfacial.h"

double calc_f2x(double phix, double phiy, double phiz, double dphi[3][3], double ddphi[3][3][3])
{
	double f2x = (2*(-(phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*dphi[0][2]*
          (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
            pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
          pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
            pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
            pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
          (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
            pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))) - 
       phix*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*dphi[2][2]*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))) + 
       2*phix*phiz*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))) - 
       phiy*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*dphi[0][1]*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))) - 
       phix*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*dphi[1][1]*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))) + 
       2*phix*phiy*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))) - 
       2*phix*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*dphi[0][0]*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))) + 
       2*pow(phix,2)*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))) + 
       (pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[1][2] - dphi[2][1],2) + 
          pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + pow(dphi[0][1] - dphi[1][0],2) + 
          pow(dphi[0][2] - dphi[2][0],2))*(2*phix*dphi[0][0] + phiy*(dphi[0][1] + dphi[1][0]) + 
          phiz*(dphi[0][2] + dphi[2][0]))*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))) + 
       (pow(phix,2) + pow(phiy,2) + pow(phiz,2))*dphi[0][2]*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) - 
       2*phix*pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2)*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) + 
       phix*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[0][2],2) + pow(dphi[1][2],2) + pow(dphi[2][2],2) + 
          phix*ddphi[0][2][2] + phiy*ddphi[1][2][2] + phiz*ddphi[2][2][2])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) + 
       (pow(phix,2) + pow(phiy,2) + pow(phiz,2))*dphi[0][1]*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) - 
       2*phix*pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2)*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) + 
       phix*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[0][1],2) + pow(dphi[1][1],2) + pow(dphi[2][1],2) + 
          phix*ddphi[0][1][1] + phiy*ddphi[1][1][1] + phiz*ddphi[2][1][1])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) + 
       (pow(phix,2) + pow(phiy,2) + pow(phiz,2))*dphi[0][0]*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) - 
       2*phix*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2)*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) - 
       (pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[1][2] - dphi[2][1],2) + 
          pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + pow(dphi[0][1] - dphi[1][0],2) + 
          pow(dphi[0][2] - dphi[2][0],2))*(dphi[0][2]*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2]) + 
          dphi[0][1]*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1]) + 
          dphi[0][0]*(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0]))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) + 
       2*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*(dphi[0][1] - dphi[1][0])*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) - 
       phix*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) + 
       2*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*(dphi[0][2] - dphi[2][0])*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) + 
       2*(dphi[2][2] + dphi[1][1] + dphi[0][0])*(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) - 
       2*phix*phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        ((dphi[1][2] - dphi[2][1])*(ddphi[1][2][2] - ddphi[2][1][2]) + 
          (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][2][2] + ddphi[1][1][2] + ddphi[0][0][2]) + 
          (dphi[0][1] - dphi[1][0])*(ddphi[0][1][2] - ddphi[1][0][2]) + 
          (dphi[0][2] - dphi[2][0])*(ddphi[0][2][2] - ddphi[2][0][2])) + 
       2*phix*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((dphi[1][2] - dphi[2][1])*(ddphi[1][2][2] - ddphi[2][1][2]) + 
          (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][2][2] + ddphi[1][1][2] + ddphi[0][0][2]) + 
          (dphi[0][1] - dphi[1][0])*(ddphi[0][1][2] - ddphi[1][0][2]) + 
          (dphi[0][2] - dphi[2][0])*(ddphi[0][2][2] - ddphi[2][0][2])) - 
       (pow(phix,2) + pow(phiy,2) + pow(phiz,2))*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        (ddphi[0][2][2] - ddphi[2][0][2]) + 2*phix*phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (pow(dphi[0][2],2) + pow(dphi[1][2],2) + pow(dphi[2][2],2) + phix*ddphi[0][2][2] + phiy*ddphi[1][2][2] + 
             phiz*ddphi[2][2][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
             phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + 
          (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
             phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2])) - 
       4*phix*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (pow(dphi[0][2],2) + pow(dphi[1][2],2) + pow(dphi[2][2],2) + phix*ddphi[0][2][2] + phiy*ddphi[1][2][2] + 
             phiz*ddphi[2][2][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
             phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + 
          (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
             phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2])) + 
       2*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(dphi[0][2] - dphi[2][0])*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (pow(dphi[0][2],2) + pow(dphi[1][2],2) + pow(dphi[2][2],2) + phix*ddphi[0][2][2] + phiy*ddphi[1][2][2] + 
             phiz*ddphi[2][2][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
             phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + 
          (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
             phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2])) - 
       phix*phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[1][2] - dphi[2][1],2) + 
          pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + pow(dphi[0][1] - dphi[1][0],2) + 
          pow(dphi[0][2] - dphi[2][0],2))*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (2*phiz*pow(dphi[2][2],2) + pow(phiz,2)*ddphi[2][2][2] + 2*phiy*dphi[1][2]*dphi[1][1] + 
          phiz*dphi[1][2]*(dphi[1][2] + dphi[2][1]) + phiy*dphi[2][2]*(dphi[1][2] + dphi[2][1]) + 
          pow(phiy,2)*ddphi[1][1][2] + phiy*phiz*(ddphi[1][2][2] + ddphi[2][1][2]) + 2*phix*dphi[0][2]*dphi[0][0] + 
          dphi[0][2]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][2] + 
          phix*(dphi[1][2]*(dphi[0][1] + dphi[1][0]) + dphi[2][2]*(dphi[0][2] + dphi[2][0]) + 
             phiy*(ddphi[0][1][2] + ddphi[1][0][2]) + phiz*(ddphi[0][2][2] + ddphi[2][0][2]))) + 
       2*phix*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        (2*phiz*pow(dphi[2][2],2) + pow(phiz,2)*ddphi[2][2][2] + 2*phiy*dphi[1][2]*dphi[1][1] + 
          phiz*dphi[1][2]*(dphi[1][2] + dphi[2][1]) + phiy*dphi[2][2]*(dphi[1][2] + dphi[2][1]) + 
          pow(phiy,2)*ddphi[1][1][2] + phiy*phiz*(ddphi[1][2][2] + ddphi[2][1][2]) + 2*phix*dphi[0][2]*dphi[0][0] + 
          dphi[0][2]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][2] + 
          phix*(dphi[1][2]*(dphi[0][1] + dphi[1][0]) + dphi[2][2]*(dphi[0][2] + dphi[2][0]) + 
             phiy*(ddphi[0][1][2] + ddphi[1][0][2]) + phiz*(ddphi[0][2][2] + ddphi[2][0][2]))) - 
       2*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(dphi[0][2] - dphi[2][0])*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        (2*phiz*pow(dphi[2][2],2) + pow(phiz,2)*ddphi[2][2][2] + 2*phiy*dphi[1][2]*dphi[1][1] + 
          phiz*dphi[1][2]*(dphi[1][2] + dphi[2][1]) + phiy*dphi[2][2]*(dphi[1][2] + dphi[2][1]) + 
          pow(phiy,2)*ddphi[1][1][2] + phiy*phiz*(ddphi[1][2][2] + ddphi[2][1][2]) + 2*phix*dphi[0][2]*dphi[0][0] + 
          dphi[0][2]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][2] + 
          phix*(dphi[1][2]*(dphi[0][1] + dphi[1][0]) + dphi[2][2]*(dphi[0][2] + dphi[2][0]) + 
             phiy*(ddphi[0][1][2] + ddphi[1][0][2]) + phiz*(ddphi[0][2][2] + ddphi[2][0][2]))) - 
       (pow(phix,2) + pow(phiy,2) + pow(phiz,2))*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        (ddphi[0][1][1] - ddphi[1][0][1]) - 2*phix*phiy*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        ((dphi[1][2] - dphi[2][1])*(ddphi[1][1][2] - ddphi[2][1][1]) + 
          (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][1][2] + ddphi[1][1][1] + ddphi[0][0][1]) + 
          (dphi[0][1] - dphi[1][0])*(ddphi[0][1][1] - ddphi[1][0][1]) + 
          (dphi[0][2] - dphi[2][0])*(ddphi[0][1][2] - ddphi[2][0][1])) + 
       2*phix*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((dphi[1][2] - dphi[2][1])*(ddphi[1][1][2] - ddphi[2][1][1]) + 
          (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][1][2] + ddphi[1][1][1] + ddphi[0][0][1]) + 
          (dphi[0][1] - dphi[1][0])*(ddphi[0][1][1] - ddphi[1][0][1]) + 
          (dphi[0][2] - dphi[2][0])*(ddphi[0][1][2] - ddphi[2][0][1])) + 
       2*phix*phiy*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[1][2] - dphi[2][1],2) + 
          pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + pow(dphi[0][1] - dphi[1][0],2) + 
          pow(dphi[0][2] - dphi[2][0],2))*(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
             phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + 
          (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (pow(dphi[0][1],2) + pow(dphi[1][1],2) + pow(dphi[2][1],2) + phix*ddphi[0][1][1] + phiy*ddphi[1][1][1] + 
             phiz*ddphi[2][1][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
             phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1])) - 
       4*phix*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
             phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + 
          (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (pow(dphi[0][1],2) + pow(dphi[1][1],2) + pow(dphi[2][1],2) + phix*ddphi[0][1][1] + phiy*ddphi[1][1][1] + 
             phiz*ddphi[2][1][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
             phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1])) + 
       2*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(dphi[0][1] - dphi[1][0])*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
             phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + 
          (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (pow(dphi[0][1],2) + pow(dphi[1][1],2) + pow(dphi[2][1],2) + phix*ddphi[0][1][1] + phiy*ddphi[1][1][1] + 
             phiz*ddphi[2][1][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
             phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1])) - 
       phix*phiy*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[1][2] - dphi[2][1],2) + 
          pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + pow(dphi[0][1] - dphi[1][0],2) + 
          pow(dphi[0][2] - dphi[2][0],2))*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (2*phiy*pow(dphi[1][1],2) + 2*phiz*dphi[2][2]*dphi[2][1] + phiz*dphi[1][1]*(dphi[1][2] + dphi[2][1]) + 
          phiy*dphi[2][1]*(dphi[1][2] + dphi[2][1]) + pow(phiz,2)*ddphi[2][1][2] + pow(phiy,2)*ddphi[1][1][1] + 
          phiy*phiz*(ddphi[1][1][2] + ddphi[2][1][1]) + 2*phix*dphi[0][1]*dphi[0][0] + 
          dphi[0][1]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][1] + 
          phix*(dphi[1][1]*(dphi[0][1] + dphi[1][0]) + dphi[2][1]*(dphi[0][2] + dphi[2][0]) + 
             phiy*(ddphi[0][1][1] + ddphi[1][0][1]) + phiz*(ddphi[0][1][2] + ddphi[2][0][1]))) + 
       2*phix*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        (2*phiy*pow(dphi[1][1],2) + 2*phiz*dphi[2][2]*dphi[2][1] + phiz*dphi[1][1]*(dphi[1][2] + dphi[2][1]) + 
          phiy*dphi[2][1]*(dphi[1][2] + dphi[2][1]) + pow(phiz,2)*ddphi[2][1][2] + pow(phiy,2)*ddphi[1][1][1] + 
          phiy*phiz*(ddphi[1][1][2] + ddphi[2][1][1]) + 2*phix*dphi[0][1]*dphi[0][0] + 
          dphi[0][1]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][1] + 
          phix*(dphi[1][1]*(dphi[0][1] + dphi[1][0]) + dphi[2][1]*(dphi[0][2] + dphi[2][0]) + 
             phiy*(ddphi[0][1][1] + ddphi[1][0][1]) + phiz*(ddphi[0][1][2] + ddphi[2][0][1]))) - 
       2*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(dphi[0][1] - dphi[1][0])*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        (2*phiy*pow(dphi[1][1],2) + 2*phiz*dphi[2][2]*dphi[2][1] + phiz*dphi[1][1]*(dphi[1][2] + dphi[2][1]) + 
          phiy*dphi[2][1]*(dphi[1][2] + dphi[2][1]) + pow(phiz,2)*ddphi[2][1][2] + pow(phiy,2)*ddphi[1][1][1] + 
          phiy*phiz*(ddphi[1][1][2] + ddphi[2][1][1]) + 2*phix*dphi[0][1]*dphi[0][0] + 
          dphi[0][1]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][1] + 
          phix*(dphi[1][1]*(dphi[0][1] + dphi[1][0]) + dphi[2][1]*(dphi[0][2] + dphi[2][0]) + 
             phiy*(ddphi[0][1][1] + ddphi[1][0][1]) + phiz*(ddphi[0][1][2] + ddphi[2][0][1]))) - 
       (pow(phix,2) + pow(phiy,2) + pow(phiz,2))*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        (ddphi[2][0][2] + ddphi[1][0][1] + ddphi[0][0][0]) - 
       2*pow(phix,2)*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        ((dphi[1][2] - dphi[2][1])*(ddphi[1][0][2] - ddphi[2][0][1]) + 
          (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][0][2] + ddphi[1][0][1] + ddphi[0][0][0]) + 
          (dphi[0][1] - dphi[1][0])*(ddphi[0][0][1] - ddphi[1][0][0]) + 
          (dphi[0][2] - dphi[2][0])*(ddphi[0][0][2] - ddphi[2][0][0])) + 
       2*phix*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((dphi[1][2] - dphi[2][1])*(ddphi[1][0][2] - ddphi[2][0][1]) + 
          (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][0][2] + ddphi[1][0][1] + ddphi[0][0][0]) + 
          (dphi[0][1] - dphi[1][0])*(ddphi[0][0][1] - ddphi[1][0][0]) + 
          (dphi[0][2] - dphi[2][0])*(ddphi[0][0][2] - ddphi[2][0][0])) + 
       phix*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[1][2] - dphi[2][1],2) + 
          pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + pow(dphi[0][1] - dphi[1][0],2) + 
          pow(dphi[0][2] - dphi[2][0],2))*(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        (pow(dphi[0][0],2) + pow(dphi[1][0],2) + pow(dphi[2][0],2) + phix*ddphi[0][0][0] + phiy*ddphi[1][0][0] + 
          phiz*ddphi[2][0][0]) + 2*pow(phix,2)*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
             phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2]) + 
          (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
             phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1]) + 
          (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (pow(dphi[0][0],2) + pow(dphi[1][0],2) + pow(dphi[2][0],2) + phix*ddphi[0][0][0] + phiy*ddphi[1][0][0] + 
             phiz*ddphi[2][0][0])) - 4*phix*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
             phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2]) + 
          (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
             phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1]) + 
          (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (pow(dphi[0][0],2) + pow(dphi[1][0],2) + pow(dphi[2][0],2) + phix*ddphi[0][0][0] + phiy*ddphi[1][0][0] + 
             phiz*ddphi[2][0][0])) + 2*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(dphi[2][2] + dphi[1][1] + dphi[0][0])*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
             phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2]) + 
          (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
             phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1]) + 
          (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (pow(dphi[0][0],2) + pow(dphi[1][0],2) + pow(dphi[2][0],2) + phix*ddphi[0][0][0] + phiy*ddphi[1][0][0] + 
             phiz*ddphi[2][0][0])) - pow(phix,2)*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (2*phix*pow(dphi[0][0],2) + 2*phiy*dphi[1][1]*dphi[1][0] + phiz*(dphi[1][2] + dphi[2][1])*dphi[1][0] + 
          2*phiz*dphi[2][2]*dphi[2][0] + phiy*(dphi[1][2] + dphi[2][1])*dphi[2][0] + 
          dphi[0][0]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phiz,2)*ddphi[2][0][2] + 
          pow(phiy,2)*ddphi[1][0][1] + phiy*phiz*(ddphi[1][0][2] + ddphi[2][0][1]) + pow(phix,2)*ddphi[0][0][0] + 
          phix*(dphi[1][0]*(dphi[0][1] + dphi[1][0]) + dphi[2][0]*(dphi[0][2] + dphi[2][0]) + 
             phiy*(ddphi[0][0][1] + ddphi[1][0][0]) + phiz*(ddphi[0][0][2] + ddphi[2][0][0]))) + 
       2*phix*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[1][2] - dphi[2][1],2) + 
          pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + pow(dphi[0][1] - dphi[1][0],2) + 
          pow(dphi[0][2] - dphi[2][0],2))*(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        (2*phix*pow(dphi[0][0],2) + 2*phiy*dphi[1][1]*dphi[1][0] + phiz*(dphi[1][2] + dphi[2][1])*dphi[1][0] + 
          2*phiz*dphi[2][2]*dphi[2][0] + phiy*(dphi[1][2] + dphi[2][1])*dphi[2][0] + 
          dphi[0][0]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phiz,2)*ddphi[2][0][2] + 
          pow(phiy,2)*ddphi[1][0][1] + phiy*phiz*(ddphi[1][0][2] + ddphi[2][0][1]) + pow(phix,2)*ddphi[0][0][0] + 
          phix*(dphi[1][0]*(dphi[0][1] + dphi[1][0]) + dphi[2][0]*(dphi[0][2] + dphi[2][0]) + 
             phiy*(ddphi[0][0][1] + ddphi[1][0][0]) + phiz*(ddphi[0][0][2] + ddphi[2][0][0]))) - 
       2*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(dphi[2][2] + dphi[1][1] + dphi[0][0])*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        (2*phix*pow(dphi[0][0],2) + 2*phiy*dphi[1][1]*dphi[1][0] + phiz*(dphi[1][2] + dphi[2][1])*dphi[1][0] + 
          2*phiz*dphi[2][2]*dphi[2][0] + phiy*(dphi[1][2] + dphi[2][1])*dphi[2][0] + 
          dphi[0][0]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phiz,2)*ddphi[2][0][2] + 
          pow(phiy,2)*ddphi[1][0][1] + phiy*phiz*(ddphi[1][0][2] + ddphi[2][0][1]) + pow(phix,2)*ddphi[0][0][0] + 
          phix*(dphi[1][0]*(dphi[0][1] + dphi[1][0]) + dphi[2][0]*(dphi[0][2] + dphi[2][0]) + 
             phiy*(ddphi[0][0][1] + ddphi[1][0][0]) + phiz*(ddphi[0][0][2] + ddphi[2][0][0])))))/
   (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
       pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
       pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3));

   return f2x;
}

double calc_f2y(double phix, double phiy, double phiz, double dphi[3][3], double ddphi[3][3][3])
{
	double f2y = (2*(-(phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*dphi[1][2]*
          (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
            pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
          pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
            pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
            pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
          (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
            pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))) - 
       phiy*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*dphi[2][2]*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))) + 
       2*phiy*phiz*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))) - 
       2*phiy*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*dphi[1][1]*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))) + 
       2*pow(phiy,2)*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))) - 
       phiy*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*dphi[0][0]*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))) - 
       phix*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*dphi[1][0]*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))) + 
       (pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(2*phiy*dphi[1][1] + phiz*(dphi[1][2] + dphi[2][1]) + 
          phix*(dphi[0][1] + dphi[1][0]))*(pow(dphi[1][2] - dphi[2][1],2) + 
          pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + pow(dphi[0][1] - dphi[1][0],2) + 
          pow(dphi[0][2] - dphi[2][0],2))*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))) + 
       2*phix*phiy*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))) + 
       (pow(phix,2) + pow(phiy,2) + pow(phiz,2))*dphi[1][2]*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) - 
       2*phiy*pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2)*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) + 
       phiy*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[0][2],2) + pow(dphi[1][2],2) + pow(dphi[2][2],2) + 
          phix*ddphi[0][2][2] + phiy*ddphi[1][2][2] + phiz*ddphi[2][2][2])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) + 
       (pow(phix,2) + pow(phiy,2) + pow(phiz,2))*dphi[1][1]*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) - 
       2*phiy*pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2)*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) + 
       phiy*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[0][1],2) + pow(dphi[1][1],2) + pow(dphi[2][1],2) + 
          phix*ddphi[0][1][1] + phiy*ddphi[1][1][1] + phiz*ddphi[2][1][1])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) + 
       (pow(phix,2) + pow(phiy,2) + pow(phiz,2))*dphi[1][0]*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) - 
       2*phiy*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2)*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) - 
       (pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[1][2] - dphi[2][1],2) + 
          pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + pow(dphi[0][1] - dphi[1][0],2) + 
          pow(dphi[0][2] - dphi[2][0],2))*(dphi[1][2]*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2]) + 
          dphi[1][1]*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1]) + 
          dphi[1][0]*(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0]))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) + 
       2*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*(dphi[1][2] - dphi[2][1])*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) - 
       (pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(ddphi[1][2][2] - ddphi[2][1][2])*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) + 
       2*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*(dphi[2][2] + dphi[1][1] + dphi[0][0])*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) - 
       phiy*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) - 
       2*(dphi[0][1] - dphi[1][0])*(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) - 
       2*phiy*phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        ((dphi[1][2] - dphi[2][1])*(ddphi[1][2][2] - ddphi[2][1][2]) + 
          (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][2][2] + ddphi[1][1][2] + ddphi[0][0][2]) + 
          (dphi[0][1] - dphi[1][0])*(ddphi[0][1][2] - ddphi[1][0][2]) + 
          (dphi[0][2] - dphi[2][0])*(ddphi[0][2][2] - ddphi[2][0][2])) + 
       2*phiy*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((dphi[1][2] - dphi[2][1])*(ddphi[1][2][2] - ddphi[2][1][2]) + 
          (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][2][2] + ddphi[1][1][2] + ddphi[0][0][2]) + 
          (dphi[0][1] - dphi[1][0])*(ddphi[0][1][2] - ddphi[1][0][2]) + 
          (dphi[0][2] - dphi[2][0])*(ddphi[0][2][2] - ddphi[2][0][2])) + 
       2*phiy*phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[1][2] - dphi[2][1],2) + 
          pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + pow(dphi[0][1] - dphi[1][0],2) + 
          pow(dphi[0][2] - dphi[2][0],2))*(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (pow(dphi[0][2],2) + pow(dphi[1][2],2) + pow(dphi[2][2],2) + phix*ddphi[0][2][2] + phiy*ddphi[1][2][2] + 
             phiz*ddphi[2][2][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
             phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + 
          (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
             phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2])) - 
       4*phiy*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (pow(dphi[0][2],2) + pow(dphi[1][2],2) + pow(dphi[2][2],2) + phix*ddphi[0][2][2] + phiy*ddphi[1][2][2] + 
             phiz*ddphi[2][2][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
             phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + 
          (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
             phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2])) + 
       2*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(dphi[1][2] - dphi[2][1])*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (pow(dphi[0][2],2) + pow(dphi[1][2],2) + pow(dphi[2][2],2) + phix*ddphi[0][2][2] + phiy*ddphi[1][2][2] + 
             phiz*ddphi[2][2][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
             phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + 
          (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
             phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2])) - 
       phiy*phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[1][2] - dphi[2][1],2) + 
          pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + pow(dphi[0][1] - dphi[1][0],2) + 
          pow(dphi[0][2] - dphi[2][0],2))*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (2*phiz*pow(dphi[2][2],2) + pow(phiz,2)*ddphi[2][2][2] + 2*phiy*dphi[1][2]*dphi[1][1] + 
          phiz*dphi[1][2]*(dphi[1][2] + dphi[2][1]) + phiy*dphi[2][2]*(dphi[1][2] + dphi[2][1]) + 
          pow(phiy,2)*ddphi[1][1][2] + phiy*phiz*(ddphi[1][2][2] + ddphi[2][1][2]) + 2*phix*dphi[0][2]*dphi[0][0] + 
          dphi[0][2]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + 
          pow(phix,2)*ddphi[0][0][2] + phix*(dphi[1][2]*(dphi[0][1] + dphi[1][0]) + 
             dphi[2][2]*(dphi[0][2] + dphi[2][0]) + phiy*(ddphi[0][1][2] + ddphi[1][0][2]) + 
             phiz*(ddphi[0][2][2] + ddphi[2][0][2]))) + 2*phiy*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*
        (phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        (2*phiz*pow(dphi[2][2],2) + pow(phiz,2)*ddphi[2][2][2] + 2*phiy*dphi[1][2]*dphi[1][1] + 
          phiz*dphi[1][2]*(dphi[1][2] + dphi[2][1]) + phiy*dphi[2][2]*(dphi[1][2] + dphi[2][1]) + 
          pow(phiy,2)*ddphi[1][1][2] + phiy*phiz*(ddphi[1][2][2] + ddphi[2][1][2]) + 2*phix*dphi[0][2]*dphi[0][0] + 
          dphi[0][2]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + 
          pow(phix,2)*ddphi[0][0][2] + phix*(dphi[1][2]*(dphi[0][1] + dphi[1][0]) + 
             dphi[2][2]*(dphi[0][2] + dphi[2][0]) + phiy*(ddphi[0][1][2] + ddphi[1][0][2]) + 
             phiz*(ddphi[0][2][2] + ddphi[2][0][2]))) - 2*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(dphi[1][2] - dphi[2][1])*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        (2*phiz*pow(dphi[2][2],2) + pow(phiz,2)*ddphi[2][2][2] + 2*phiy*dphi[1][2]*dphi[1][1] + 
          phiz*dphi[1][2]*(dphi[1][2] + dphi[2][1]) + phiy*dphi[2][2]*(dphi[1][2] + dphi[2][1]) + 
          pow(phiy,2)*ddphi[1][1][2] + phiy*phiz*(ddphi[1][2][2] + ddphi[2][1][2]) + 2*phix*dphi[0][2]*dphi[0][0] + 
          dphi[0][2]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + 
          pow(phix,2)*ddphi[0][0][2] + phix*(dphi[1][2]*(dphi[0][1] + dphi[1][0]) + 
             dphi[2][2]*(dphi[0][2] + dphi[2][0]) + phiy*(ddphi[0][1][2] + ddphi[1][0][2]) + 
             phiz*(ddphi[0][2][2] + ddphi[2][0][2]))) - (pow(phix,2) + pow(phiy,2) + pow(phiz,2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        (ddphi[2][1][2] + ddphi[1][1][1] + ddphi[0][0][1]) - 
       2*pow(phiy,2)*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        ((dphi[1][2] - dphi[2][1])*(ddphi[1][1][2] - ddphi[2][1][1]) + 
          (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][1][2] + ddphi[1][1][1] + ddphi[0][0][1]) + 
          (dphi[0][1] - dphi[1][0])*(ddphi[0][1][1] - ddphi[1][0][1]) + 
          (dphi[0][2] - dphi[2][0])*(ddphi[0][1][2] - ddphi[2][0][1])) + 
       2*phiy*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((dphi[1][2] - dphi[2][1])*(ddphi[1][1][2] - ddphi[2][1][1]) + 
          (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][1][2] + ddphi[1][1][1] + ddphi[0][0][1]) + 
          (dphi[0][1] - dphi[1][0])*(ddphi[0][1][1] - ddphi[1][0][1]) + 
          (dphi[0][2] - dphi[2][0])*(ddphi[0][1][2] - ddphi[2][0][1])) + 
       2*pow(phiy,2)*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[1][2] - dphi[2][1],2) + 
          pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + pow(dphi[0][1] - dphi[1][0],2) + 
          pow(dphi[0][2] - dphi[2][0],2))*(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
             phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + 
          (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (pow(dphi[0][1],2) + pow(dphi[1][1],2) + pow(dphi[2][1],2) + phix*ddphi[0][1][1] + phiy*ddphi[1][1][1] + 
             phiz*ddphi[2][1][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
             phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1])) - 
       4*phiy*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
             phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + 
          (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (pow(dphi[0][1],2) + pow(dphi[1][1],2) + pow(dphi[2][1],2) + phix*ddphi[0][1][1] + phiy*ddphi[1][1][1] + 
             phiz*ddphi[2][1][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
             phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1])) + 
       2*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(dphi[2][2] + dphi[1][1] + dphi[0][0])*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
             phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + 
          (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (pow(dphi[0][1],2) + pow(dphi[1][1],2) + pow(dphi[2][1],2) + phix*ddphi[0][1][1] + phiy*ddphi[1][1][1] + 
             phiz*ddphi[2][1][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
             phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1])) - 
       pow(phiy,2)*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[1][2] - dphi[2][1],2) + 
          pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + pow(dphi[0][1] - dphi[1][0],2) + 
          pow(dphi[0][2] - dphi[2][0],2))*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (2*phiy*pow(dphi[1][1],2) + 2*phiz*dphi[2][2]*dphi[2][1] + phiz*dphi[1][1]*(dphi[1][2] + dphi[2][1]) + 
          phiy*dphi[2][1]*(dphi[1][2] + dphi[2][1]) + pow(phiz,2)*ddphi[2][1][2] + pow(phiy,2)*ddphi[1][1][1] + 
          phiy*phiz*(ddphi[1][1][2] + ddphi[2][1][1]) + 2*phix*dphi[0][1]*dphi[0][0] + 
          dphi[0][1]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + 
          pow(phix,2)*ddphi[0][0][1] + phix*(dphi[1][1]*(dphi[0][1] + dphi[1][0]) + 
             dphi[2][1]*(dphi[0][2] + dphi[2][0]) + phiy*(ddphi[0][1][1] + ddphi[1][0][1]) + 
             phiz*(ddphi[0][1][2] + ddphi[2][0][1]))) + 2*phiy*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*
        (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        (2*phiy*pow(dphi[1][1],2) + 2*phiz*dphi[2][2]*dphi[2][1] + phiz*dphi[1][1]*(dphi[1][2] + dphi[2][1]) + 
          phiy*dphi[2][1]*(dphi[1][2] + dphi[2][1]) + pow(phiz,2)*ddphi[2][1][2] + pow(phiy,2)*ddphi[1][1][1] + 
          phiy*phiz*(ddphi[1][1][2] + ddphi[2][1][1]) + 2*phix*dphi[0][1]*dphi[0][0] + 
          dphi[0][1]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + 
          pow(phix,2)*ddphi[0][0][1] + phix*(dphi[1][1]*(dphi[0][1] + dphi[1][0]) + 
             dphi[2][1]*(dphi[0][2] + dphi[2][0]) + phiy*(ddphi[0][1][1] + ddphi[1][0][1]) + 
             phiz*(ddphi[0][1][2] + ddphi[2][0][1]))) - 2*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*
        (dphi[2][2] + dphi[1][1] + dphi[0][0])*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        (2*phiy*pow(dphi[1][1],2) + 2*phiz*dphi[2][2]*dphi[2][1] + phiz*dphi[1][1]*(dphi[1][2] + dphi[2][1]) + 
          phiy*dphi[2][1]*(dphi[1][2] + dphi[2][1]) + pow(phiz,2)*ddphi[2][1][2] + pow(phiy,2)*ddphi[1][1][1] + 
          phiy*phiz*(ddphi[1][1][2] + ddphi[2][1][1]) + 2*phix*dphi[0][1]*dphi[0][0] + 
          dphi[0][1]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + 
          pow(phix,2)*ddphi[0][0][1] + phix*(dphi[1][1]*(dphi[0][1] + dphi[1][0]) + 
             dphi[2][1]*(dphi[0][2] + dphi[2][0]) + phiy*(ddphi[0][1][1] + ddphi[1][0][1]) + 
             phiz*(ddphi[0][1][2] + ddphi[2][0][1]))) + (pow(phix,2) + pow(phiy,2) + pow(phiz,2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        (ddphi[0][0][1] - ddphi[1][0][0]) - 2*phix*phiy*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        ((dphi[1][2] - dphi[2][1])*(ddphi[1][0][2] - ddphi[2][0][1]) + 
          (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][0][2] + ddphi[1][0][1] + ddphi[0][0][0]) + 
          (dphi[0][1] - dphi[1][0])*(ddphi[0][0][1] - ddphi[1][0][0]) + 
          (dphi[0][2] - dphi[2][0])*(ddphi[0][0][2] - ddphi[2][0][0])) + 
       2*phiy*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((dphi[1][2] - dphi[2][1])*(ddphi[1][0][2] - ddphi[2][0][1]) + 
          (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][0][2] + ddphi[1][0][1] + ddphi[0][0][0]) + 
          (dphi[0][1] - dphi[1][0])*(ddphi[0][0][1] - ddphi[1][0][0]) + 
          (dphi[0][2] - dphi[2][0])*(ddphi[0][0][2] - ddphi[2][0][0])) + 
       phiy*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[1][2] - dphi[2][1],2) + 
          pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + pow(dphi[0][1] - dphi[1][0],2) + 
          pow(dphi[0][2] - dphi[2][0],2))*(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        (pow(dphi[0][0],2) + pow(dphi[1][0],2) + pow(dphi[2][0],2) + phix*ddphi[0][0][0] + phiy*ddphi[1][0][0] + 
          phiz*ddphi[2][0][0]) + 2*phix*phiy*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
             phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2]) + 
          (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
             phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1]) + 
          (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (pow(dphi[0][0],2) + pow(dphi[1][0],2) + pow(dphi[2][0],2) + phix*ddphi[0][0][0] + phiy*ddphi[1][0][0] + 
             phiz*ddphi[2][0][0])) - 4*phiy*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
             phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2]) + 
          (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
             phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1]) + 
          (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (pow(dphi[0][0],2) + pow(dphi[1][0],2) + pow(dphi[2][0],2) + phix*ddphi[0][0][0] + phiy*ddphi[1][0][0] + 
             phiz*ddphi[2][0][0])) - 2*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(dphi[0][1] - dphi[1][0])*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
             phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2]) + 
          (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
             phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1]) + 
          (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (pow(dphi[0][0],2) + pow(dphi[1][0],2) + pow(dphi[2][0],2) + phix*ddphi[0][0][0] + phiy*ddphi[1][0][0] + 
             phiz*ddphi[2][0][0])) - phix*phiy*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (2*phix*pow(dphi[0][0],2) + 2*phiy*dphi[1][1]*dphi[1][0] + phiz*(dphi[1][2] + dphi[2][1])*dphi[1][0] + 
          2*phiz*dphi[2][2]*dphi[2][0] + phiy*(dphi[1][2] + dphi[2][1])*dphi[2][0] + 
          dphi[0][0]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + 
          pow(phiz,2)*ddphi[2][0][2] + pow(phiy,2)*ddphi[1][0][1] + phiy*phiz*(ddphi[1][0][2] + ddphi[2][0][1]) + 
          pow(phix,2)*ddphi[0][0][0] + phix*(dphi[1][0]*(dphi[0][1] + dphi[1][0]) + 
             dphi[2][0]*(dphi[0][2] + dphi[2][0]) + phiy*(ddphi[0][0][1] + ddphi[1][0][0]) + 
             phiz*(ddphi[0][0][2] + ddphi[2][0][0]))) + 2*phiy*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        (2*phix*pow(dphi[0][0],2) + 2*phiy*dphi[1][1]*dphi[1][0] + phiz*(dphi[1][2] + dphi[2][1])*dphi[1][0] + 
          2*phiz*dphi[2][2]*dphi[2][0] + phiy*(dphi[1][2] + dphi[2][1])*dphi[2][0] + 
          dphi[0][0]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + 
          pow(phiz,2)*ddphi[2][0][2] + pow(phiy,2)*ddphi[1][0][1] + phiy*phiz*(ddphi[1][0][2] + ddphi[2][0][1]) + 
          pow(phix,2)*ddphi[0][0][0] + phix*(dphi[1][0]*(dphi[0][1] + dphi[1][0]) + 
             dphi[2][0]*(dphi[0][2] + dphi[2][0]) + phiy*(ddphi[0][0][1] + ddphi[1][0][0]) + 
             phiz*(ddphi[0][0][2] + ddphi[2][0][0]))) + 2*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(dphi[0][1] - dphi[1][0])*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        (2*phix*pow(dphi[0][0],2) + 2*phiy*dphi[1][1]*dphi[1][0] + phiz*(dphi[1][2] + dphi[2][1])*dphi[1][0] + 
          2*phiz*dphi[2][2]*dphi[2][0] + phiy*(dphi[1][2] + dphi[2][1])*dphi[2][0] + 
          dphi[0][0]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + 
          pow(phiz,2)*ddphi[2][0][2] + pow(phiy,2)*ddphi[1][0][1] + phiy*phiz*(ddphi[1][0][2] + ddphi[2][0][1]) + 
          pow(phix,2)*ddphi[0][0][0] + phix*(dphi[1][0]*(dphi[0][1] + dphi[1][0]) + 
             dphi[2][0]*(dphi[0][2] + dphi[2][0]) + phiy*(ddphi[0][0][1] + ddphi[1][0][0]) + 
             phiz*(ddphi[0][0][2] + ddphi[2][0][0])))))/
   (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
       pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
       pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3));

    return f2y;
}

double calc_f2z(double phix, double phiy, double phiz, double dphi[3][3], double ddphi[3][3][3])
{
	double f2z = (2*(-2*phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*dphi[2][2]*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))) + 
       2*pow(phiz,2)*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))) - 
       phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*dphi[1][1]*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))) - 
       phiy*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*dphi[2][1]*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))) + 
       2*phiy*phiz*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))) - 
       phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*dphi[0][0]*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))) - 
       phix*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[1][2] - dphi[2][1],2) + 
          pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + pow(dphi[0][1] - dphi[1][0],2) + 
          pow(dphi[0][2] - dphi[2][0],2))*dphi[2][0]*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))) + 
       2*phix*phiz*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))) + 
       (pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[1][2] - dphi[2][1],2) + 
          pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + pow(dphi[0][1] - dphi[1][0],2) + 
          pow(dphi[0][2] - dphi[2][0],2))*(2*phiz*dphi[2][2] + phiy*(dphi[1][2] + dphi[2][1]) + 
          phix*(dphi[0][2] + dphi[2][0]))*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))) + 
       (pow(phix,2) + pow(phiy,2) + pow(phiz,2))*dphi[2][2]*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) - 
       2*phiz*pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2)*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) + 
       phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[0][2],2) + pow(dphi[1][2],2) + pow(dphi[2][2],2) + 
          phix*ddphi[0][2][2] + phiy*ddphi[1][2][2] + phiz*ddphi[2][2][2])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) + 
       (pow(phix,2) + pow(phiy,2) + pow(phiz,2))*dphi[2][1]*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) - 
       2*phiz*pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2)*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) + 
       phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[0][1],2) + pow(dphi[1][1],2) + pow(dphi[2][1],2) + 
          phix*ddphi[0][1][1] + phiy*ddphi[1][1][1] + phiz*ddphi[2][1][1])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) + 
       (pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[1][2] - dphi[2][1],2) + 
          pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + pow(dphi[0][1] - dphi[1][0],2) + 
          pow(dphi[0][2] - dphi[2][0],2))*dphi[2][0]*(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) - 
       2*phiz*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2)*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) - 
       (pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[1][2] - dphi[2][1],2) + 
          pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + pow(dphi[0][1] - dphi[1][0],2) + 
          pow(dphi[0][2] - dphi[2][0],2))*(dphi[2][2]*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2]) + 
          dphi[2][1]*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1]) + 
          dphi[2][0]*(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0]))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) - 
       2*(dphi[1][2] - dphi[2][1])*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) + 
       (pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(ddphi[1][1][2] - ddphi[2][1][1])*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) + 
       2*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*(dphi[2][2] + dphi[1][1] + dphi[0][0])*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) - 
       phiz*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) - 
       2*(dphi[0][2] - dphi[2][0])*(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2) - 
       (pow(phix,2) + pow(phiy,2) + pow(phiz,2))*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        (ddphi[2][2][2] + ddphi[1][1][2] + ddphi[0][0][2]) - 
       2*pow(phiz,2)*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        ((dphi[1][2] - dphi[2][1])*(ddphi[1][2][2] - ddphi[2][1][2]) + 
          (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][2][2] + ddphi[1][1][2] + ddphi[0][0][2]) + 
          (dphi[0][1] - dphi[1][0])*(ddphi[0][1][2] - ddphi[1][0][2]) + 
          (dphi[0][2] - dphi[2][0])*(ddphi[0][2][2] - ddphi[2][0][2])) + 
       2*phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((dphi[1][2] - dphi[2][1])*(ddphi[1][2][2] - ddphi[2][1][2]) + 
          (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][2][2] + ddphi[1][1][2] + ddphi[0][0][2]) + 
          (dphi[0][1] - dphi[1][0])*(ddphi[0][1][2] - ddphi[1][0][2]) + 
          (dphi[0][2] - dphi[2][0])*(ddphi[0][2][2] - ddphi[2][0][2])) + 
       2*pow(phiz,2)*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[1][2] - dphi[2][1],2) + 
          pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + pow(dphi[0][1] - dphi[1][0],2) + 
          pow(dphi[0][2] - dphi[2][0],2))*(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (pow(dphi[0][2],2) + pow(dphi[1][2],2) + pow(dphi[2][2],2) + phix*ddphi[0][2][2] + phiy*ddphi[1][2][2] + 
             phiz*ddphi[2][2][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
             phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + 
          (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
             phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2])) - 
       4*phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (pow(dphi[0][2],2) + pow(dphi[1][2],2) + pow(dphi[2][2],2) + phix*ddphi[0][2][2] + phiy*ddphi[1][2][2] + 
             phiz*ddphi[2][2][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
             phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + 
          (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
             phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2])) + 
       2*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(dphi[2][2] + dphi[1][1] + dphi[0][0])*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (pow(dphi[0][2],2) + pow(dphi[1][2],2) + pow(dphi[2][2],2) + phix*ddphi[0][2][2] + phiy*ddphi[1][2][2] + 
             phiz*ddphi[2][2][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
             phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + 
          (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
             phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2])) - 
       pow(phiz,2)*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[1][2] - dphi[2][1],2) + 
          pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + pow(dphi[0][1] - dphi[1][0],2) + 
          pow(dphi[0][2] - dphi[2][0],2))*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (2*phiz*pow(dphi[2][2],2) + pow(phiz,2)*ddphi[2][2][2] + 2*phiy*dphi[1][2]*dphi[1][1] + 
          phiz*dphi[1][2]*(dphi[1][2] + dphi[2][1]) + phiy*dphi[2][2]*(dphi[1][2] + dphi[2][1]) + 
          pow(phiy,2)*ddphi[1][1][2] + phiy*phiz*(ddphi[1][2][2] + ddphi[2][1][2]) + 2*phix*dphi[0][2]*dphi[0][0] + 
          dphi[0][2]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][2] + 
          phix*(dphi[1][2]*(dphi[0][1] + dphi[1][0]) + dphi[2][2]*(dphi[0][2] + dphi[2][0]) + 
             phiy*(ddphi[0][1][2] + ddphi[1][0][2]) + phiz*(ddphi[0][2][2] + ddphi[2][0][2]))) + 
       2*phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        (2*phiz*pow(dphi[2][2],2) + pow(phiz,2)*ddphi[2][2][2] + 2*phiy*dphi[1][2]*dphi[1][1] + 
          phiz*dphi[1][2]*(dphi[1][2] + dphi[2][1]) + phiy*dphi[2][2]*(dphi[1][2] + dphi[2][1]) + 
          pow(phiy,2)*ddphi[1][1][2] + phiy*phiz*(ddphi[1][2][2] + ddphi[2][1][2]) + 2*phix*dphi[0][2]*dphi[0][0] + 
          dphi[0][2]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][2] + 
          phix*(dphi[1][2]*(dphi[0][1] + dphi[1][0]) + dphi[2][2]*(dphi[0][2] + dphi[2][0]) + 
             phiy*(ddphi[0][1][2] + ddphi[1][0][2]) + phiz*(ddphi[0][2][2] + ddphi[2][0][2]))) - 
       2*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(dphi[2][2] + dphi[1][1] + dphi[0][0])*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        (2*phiz*pow(dphi[2][2],2) + pow(phiz,2)*ddphi[2][2][2] + 2*phiy*dphi[1][2]*dphi[1][1] + 
          phiz*dphi[1][2]*(dphi[1][2] + dphi[2][1]) + phiy*dphi[2][2]*(dphi[1][2] + dphi[2][1]) + 
          pow(phiy,2)*ddphi[1][1][2] + phiy*phiz*(ddphi[1][2][2] + ddphi[2][1][2]) + 2*phix*dphi[0][2]*dphi[0][0] + 
          dphi[0][2]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][2] + 
          phix*(dphi[1][2]*(dphi[0][1] + dphi[1][0]) + dphi[2][2]*(dphi[0][2] + dphi[2][0]) + 
             phiy*(ddphi[0][1][2] + ddphi[1][0][2]) + phiz*(ddphi[0][2][2] + ddphi[2][0][2]))) - 
       2*phiy*phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        ((dphi[1][2] - dphi[2][1])*(ddphi[1][1][2] - ddphi[2][1][1]) + 
          (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][1][2] + ddphi[1][1][1] + ddphi[0][0][1]) + 
          (dphi[0][1] - dphi[1][0])*(ddphi[0][1][1] - ddphi[1][0][1]) + 
          (dphi[0][2] - dphi[2][0])*(ddphi[0][1][2] - ddphi[2][0][1])) + 
       2*phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((dphi[1][2] - dphi[2][1])*(ddphi[1][1][2] - ddphi[2][1][1]) + 
          (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][1][2] + ddphi[1][1][1] + ddphi[0][0][1]) + 
          (dphi[0][1] - dphi[1][0])*(ddphi[0][1][1] - ddphi[1][0][1]) + 
          (dphi[0][2] - dphi[2][0])*(ddphi[0][1][2] - ddphi[2][0][1])) + 
       2*phiy*phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[1][2] - dphi[2][1],2) + 
          pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + pow(dphi[0][1] - dphi[1][0],2) + 
          pow(dphi[0][2] - dphi[2][0],2))*(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
             phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + 
          (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (pow(dphi[0][1],2) + pow(dphi[1][1],2) + pow(dphi[2][1],2) + phix*ddphi[0][1][1] + phiy*ddphi[1][1][1] + 
             phiz*ddphi[2][1][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
             phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1])) - 
       4*phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
             phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + 
          (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (pow(dphi[0][1],2) + pow(dphi[1][1],2) + pow(dphi[2][1],2) + phix*ddphi[0][1][1] + phiy*ddphi[1][1][1] + 
             phiz*ddphi[2][1][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
             phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1])) - 
       2*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(dphi[1][2] - dphi[2][1])*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
             phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + 
          (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (pow(dphi[0][1],2) + pow(dphi[1][1],2) + pow(dphi[2][1],2) + phix*ddphi[0][1][1] + phiy*ddphi[1][1][1] + 
             phiz*ddphi[2][1][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
             phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1])) - 
       phiy*phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[1][2] - dphi[2][1],2) + 
          pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + pow(dphi[0][1] - dphi[1][0],2) + 
          pow(dphi[0][2] - dphi[2][0],2))*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (2*phiy*pow(dphi[1][1],2) + 2*phiz*dphi[2][2]*dphi[2][1] + phiz*dphi[1][1]*(dphi[1][2] + dphi[2][1]) + 
          phiy*dphi[2][1]*(dphi[1][2] + dphi[2][1]) + pow(phiz,2)*ddphi[2][1][2] + pow(phiy,2)*ddphi[1][1][1] + 
          phiy*phiz*(ddphi[1][1][2] + ddphi[2][1][1]) + 2*phix*dphi[0][1]*dphi[0][0] + 
          dphi[0][1]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][1] + 
          phix*(dphi[1][1]*(dphi[0][1] + dphi[1][0]) + dphi[2][1]*(dphi[0][2] + dphi[2][0]) + 
             phiy*(ddphi[0][1][1] + ddphi[1][0][1]) + phiz*(ddphi[0][1][2] + ddphi[2][0][1]))) + 
       2*phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        (2*phiy*pow(dphi[1][1],2) + 2*phiz*dphi[2][2]*dphi[2][1] + phiz*dphi[1][1]*(dphi[1][2] + dphi[2][1]) + 
          phiy*dphi[2][1]*(dphi[1][2] + dphi[2][1]) + pow(phiz,2)*ddphi[2][1][2] + pow(phiy,2)*ddphi[1][1][1] + 
          phiy*phiz*(ddphi[1][1][2] + ddphi[2][1][1]) + 2*phix*dphi[0][1]*dphi[0][0] + 
          dphi[0][1]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][1] + 
          phix*(dphi[1][1]*(dphi[0][1] + dphi[1][0]) + dphi[2][1]*(dphi[0][2] + dphi[2][0]) + 
             phiy*(ddphi[0][1][1] + ddphi[1][0][1]) + phiz*(ddphi[0][1][2] + ddphi[2][0][1]))) + 
       2*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(dphi[1][2] - dphi[2][1])*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        (2*phiy*pow(dphi[1][1],2) + 2*phiz*dphi[2][2]*dphi[2][1] + phiz*dphi[1][1]*(dphi[1][2] + dphi[2][1]) + 
          phiy*dphi[2][1]*(dphi[1][2] + dphi[2][1]) + pow(phiz,2)*ddphi[2][1][2] + pow(phiy,2)*ddphi[1][1][1] + 
          phiy*phiz*(ddphi[1][1][2] + ddphi[2][1][1]) + 2*phix*dphi[0][1]*dphi[0][0] + 
          dphi[0][1]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][1] + 
          phix*(dphi[1][1]*(dphi[0][1] + dphi[1][0]) + dphi[2][1]*(dphi[0][2] + dphi[2][0]) + 
             phiy*(ddphi[0][1][1] + ddphi[1][0][1]) + phiz*(ddphi[0][1][2] + ddphi[2][0][1]))) - 
       2*phix*phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        ((dphi[1][2] - dphi[2][1])*(ddphi[1][0][2] - ddphi[2][0][1]) + 
          (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][0][2] + ddphi[1][0][1] + ddphi[0][0][0]) + 
          (dphi[0][1] - dphi[1][0])*(ddphi[0][0][1] - ddphi[1][0][0]) + 
          (dphi[0][2] - dphi[2][0])*(ddphi[0][0][2] - ddphi[2][0][0])) + 
       2*phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((dphi[1][2] - dphi[2][1])*(ddphi[1][0][2] - ddphi[2][0][1]) + 
          (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][0][2] + ddphi[1][0][1] + ddphi[0][0][0]) + 
          (dphi[0][1] - dphi[1][0])*(ddphi[0][0][1] - ddphi[1][0][0]) + 
          (dphi[0][2] - dphi[2][0])*(ddphi[0][0][2] - ddphi[2][0][0])) + 
       (pow(phix,2) + pow(phiy,2) + pow(phiz,2))*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        (ddphi[0][0][2] - ddphi[2][0][0]) + phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        (pow(dphi[0][0],2) + pow(dphi[1][0],2) + pow(dphi[2][0],2) + phix*ddphi[0][0][0] + phiy*ddphi[1][0][0] + 
          phiz*ddphi[2][0][0]) + 2*phix*phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
             phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2]) + 
          (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
             phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1]) + 
          (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (pow(dphi[0][0],2) + pow(dphi[1][0],2) + pow(dphi[2][0],2) + phix*ddphi[0][0][0] + phiy*ddphi[1][0][0] + 
             phiz*ddphi[2][0][0])) - 4*phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
             phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2]) + 
          (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
             phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1]) + 
          (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (pow(dphi[0][0],2) + pow(dphi[1][0],2) + pow(dphi[2][0],2) + phix*ddphi[0][0][0] + phiy*ddphi[1][0][0] + 
             phiz*ddphi[2][0][0])) - 2*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(dphi[0][2] - dphi[2][0])*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
        ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
           (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
             phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2]) + 
          (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
           (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
             phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1]) + 
          (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
           (pow(dphi[0][0],2) + pow(dphi[1][0],2) + pow(dphi[2][0],2) + phix*ddphi[0][0][0] + phiy*ddphi[1][0][0] + 
             phiz*ddphi[2][0][0])) - phix*phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*
        (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
          pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (2*phix*pow(dphi[0][0],2) + 2*phiy*dphi[1][1]*dphi[1][0] + phiz*(dphi[1][2] + dphi[2][1])*dphi[1][0] + 
          2*phiz*dphi[2][2]*dphi[2][0] + phiy*(dphi[1][2] + dphi[2][1])*dphi[2][0] + 
          dphi[0][0]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phiz,2)*ddphi[2][0][2] + 
          pow(phiy,2)*ddphi[1][0][1] + phiy*phiz*(ddphi[1][0][2] + ddphi[2][0][1]) + pow(phix,2)*ddphi[0][0][0] + 
          phix*(dphi[1][0]*(dphi[0][1] + dphi[1][0]) + dphi[2][0]*(dphi[0][2] + dphi[2][0]) + 
             phiy*(ddphi[0][0][1] + ddphi[1][0][0]) + phiz*(ddphi[0][0][2] + ddphi[2][0][0]))) + 
       2*phiz*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(pow(dphi[1][2] - dphi[2][1],2) + 
          pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + pow(dphi[0][1] - dphi[1][0],2) + 
          pow(dphi[0][2] - dphi[2][0],2))*(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
        (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2))*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        (2*phix*pow(dphi[0][0],2) + 2*phiy*dphi[1][1]*dphi[1][0] + phiz*(dphi[1][2] + dphi[2][1])*dphi[1][0] + 
          2*phiz*dphi[2][2]*dphi[2][0] + phiy*(dphi[1][2] + dphi[2][1])*dphi[2][0] + 
          dphi[0][0]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phiz,2)*ddphi[2][0][2] + 
          pow(phiy,2)*ddphi[1][0][1] + phiy*phiz*(ddphi[1][0][2] + ddphi[2][0][1]) + pow(phix,2)*ddphi[0][0][0] + 
          phix*(dphi[1][0]*(dphi[0][1] + dphi[1][0]) + dphi[2][0]*(dphi[0][2] + dphi[2][0]) + 
             phiy*(ddphi[0][0][1] + ddphi[1][0][0]) + phiz*(ddphi[0][0][2] + ddphi[2][0][0]))) + 
       2*(pow(phix,2) + pow(phiy,2) + pow(phiz,2))*(dphi[0][2] - dphi[2][0])*
        pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
          pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
          pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)*
        (pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
          pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])))*
        (2*phix*pow(dphi[0][0],2) + 2*phiy*dphi[1][1]*dphi[1][0] + phiz*(dphi[1][2] + dphi[2][1])*dphi[1][0] + 
          2*phiz*dphi[2][2]*dphi[2][0] + phiy*(dphi[1][2] + dphi[2][1])*dphi[2][0] + 
          dphi[0][0]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phiz,2)*ddphi[2][0][2] + 
          pow(phiy,2)*ddphi[1][0][1] + phiy*phiz*(ddphi[1][0][2] + ddphi[2][0][1]) + pow(phix,2)*ddphi[0][0][0] + 
          phix*(dphi[1][0]*(dphi[0][1] + dphi[1][0]) + dphi[2][0]*(dphi[0][2] + dphi[2][0]) + 
             phiy*(ddphi[0][0][1] + ddphi[1][0][0]) + phiz*(ddphi[0][0][2] + ddphi[2][0][0])))))/
   (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
       pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
       pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3));

   return f2z;
}

double calc_f4x(double phix, double phiy, double phiz, double dphi[3][3], double ddphi[3][3][3])
{
	double f4x = (-4*phiz*dphi[0][2]*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (4*phix*dphi[2][2]*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (16*phix*phiz*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (4*phiy*dphi[0][1]*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (4*phix*dphi[1][1]*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (16*phix*phiy*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (8*phix*dphi[0][0]*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (16*pow(phix,2)*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (4*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      (2*phix*dphi[0][0] + phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0]))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (4*dphi[0][2]*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (16*phix*pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2)*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (4*phix*(pow(dphi[0][2],2) + pow(dphi[1][2],2) + pow(dphi[2][2],2) + phix*ddphi[0][2][2] + phiy*ddphi[1][2][2] + 
        phiz*ddphi[2][2][2])*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (4*dphi[0][1]*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (16*phix*pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2)*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (4*phix*(pow(dphi[0][1],2) + pow(dphi[1][1],2) + pow(dphi[2][1],2) + phix*ddphi[0][1][1] + phiy*ddphi[1][1][1] + 
        phiz*ddphi[2][1][1])*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (4*dphi[0][0]*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (16*phix*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2)*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (4*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      (dphi[0][2]*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2]) + 
        dphi[0][1]*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1]) + 
        dphi[0][0]*(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0]))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (8*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*(dphi[0][1] - dphi[1][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (4*phix*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (8*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*(dphi[0][2] - dphi[2][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (8*(dphi[2][2] + dphi[1][1] + dphi[0][0])*(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (8*phix*phiz*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      ((dphi[1][2] - dphi[2][1])*(ddphi[1][2][2] - ddphi[2][1][2]) + 
        (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][2][2] + ddphi[1][1][2] + ddphi[0][0][2]) + 
        (dphi[0][1] - dphi[1][0])*(ddphi[0][1][2] - ddphi[1][0][2]) + 
        (dphi[0][2] - dphi[2][0])*(ddphi[0][2][2] - ddphi[2][0][2])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (8*phix*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((dphi[1][2] - dphi[2][1])*(ddphi[1][2][2] - ddphi[2][1][2]) + 
        (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][2][2] + ddphi[1][1][2] + ddphi[0][0][2]) + 
        (dphi[0][1] - dphi[1][0])*(ddphi[0][1][2] - ddphi[1][0][2]) + 
        (dphi[0][2] - dphi[2][0])*(ddphi[0][2][2] - ddphi[2][0][2])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (2*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      (ddphi[0][2][2] - ddphi[2][0][2]))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (16*phix*phiz*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (pow(dphi[0][2],2) + pow(dphi[1][2],2) + pow(dphi[2][2],2) + phix*ddphi[0][2][2] + phiy*ddphi[1][2][2] + 
           phiz*ddphi[2][2][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
           phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
           phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (24*phix*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (pow(dphi[0][2],2) + pow(dphi[1][2],2) + pow(dphi[2][2],2) + phix*ddphi[0][2][2] + phiy*ddphi[1][2][2] + 
           phiz*ddphi[2][2][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
           phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
           phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),4)) + 
   (8*(dphi[0][2] - dphi[2][0])*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + 
        phiy*phiz*(dphi[1][2] + dphi[2][1]) + pow(phix,2)*dphi[0][0] + 
        phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (pow(dphi[0][2],2) + pow(dphi[1][2],2) + pow(dphi[2][2],2) + phix*ddphi[0][2][2] + phiy*ddphi[1][2][2] + 
           phiz*ddphi[2][2][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
           phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
           phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (12*phix*phiz*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
      (2*phiz*pow(dphi[2][2],2) + pow(phiz,2)*ddphi[2][2][2] + 2*phiy*dphi[1][2]*dphi[1][1] + 
        phiz*dphi[1][2]*(dphi[1][2] + dphi[2][1]) + phiy*dphi[2][2]*(dphi[1][2] + dphi[2][1]) + 
        pow(phiy,2)*ddphi[1][1][2] + phiy*phiz*(ddphi[1][2][2] + ddphi[2][1][2]) + 2*phix*dphi[0][2]*dphi[0][0] + 
        dphi[0][2]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][2] + 
        phix*(dphi[1][2]*(dphi[0][1] + dphi[1][0]) + dphi[2][2]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][1][2] + ddphi[1][0][2]) + phiz*(ddphi[0][2][2] + ddphi[2][0][2]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (16*phix*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      (2*phiz*pow(dphi[2][2],2) + pow(phiz,2)*ddphi[2][2][2] + 2*phiy*dphi[1][2]*dphi[1][1] + 
        phiz*dphi[1][2]*(dphi[1][2] + dphi[2][1]) + phiy*dphi[2][2]*(dphi[1][2] + dphi[2][1]) + 
        pow(phiy,2)*ddphi[1][1][2] + phiy*phiz*(ddphi[1][2][2] + ddphi[2][1][2]) + 2*phix*dphi[0][2]*dphi[0][0] + 
        dphi[0][2]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][2] + 
        phix*(dphi[1][2]*(dphi[0][1] + dphi[1][0]) + dphi[2][2]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][1][2] + ddphi[1][0][2]) + phiz*(ddphi[0][2][2] + ddphi[2][0][2]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (8*(dphi[0][2] - dphi[2][0])*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + 
        phiy*phiz*(dphi[1][2] + dphi[2][1]) + pow(phix,2)*dphi[0][0] + 
        phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      (2*phiz*pow(dphi[2][2],2) + pow(phiz,2)*ddphi[2][2][2] + 2*phiy*dphi[1][2]*dphi[1][1] + 
        phiz*dphi[1][2]*(dphi[1][2] + dphi[2][1]) + phiy*dphi[2][2]*(dphi[1][2] + dphi[2][1]) + 
        pow(phiy,2)*ddphi[1][1][2] + phiy*phiz*(ddphi[1][2][2] + ddphi[2][1][2]) + 2*phix*dphi[0][2]*dphi[0][0] + 
        dphi[0][2]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][2] + 
        phix*(dphi[1][2]*(dphi[0][1] + dphi[1][0]) + dphi[2][2]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][1][2] + ddphi[1][0][2]) + phiz*(ddphi[0][2][2] + ddphi[2][0][2]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (2*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      (ddphi[0][1][1] - ddphi[1][0][1]))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (8*phix*phiy*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      ((dphi[1][2] - dphi[2][1])*(ddphi[1][1][2] - ddphi[2][1][1]) + 
        (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][1][2] + ddphi[1][1][1] + ddphi[0][0][1]) + 
        (dphi[0][1] - dphi[1][0])*(ddphi[0][1][1] - ddphi[1][0][1]) + 
        (dphi[0][2] - dphi[2][0])*(ddphi[0][1][2] - ddphi[2][0][1])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (8*phix*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((dphi[1][2] - dphi[2][1])*(ddphi[1][1][2] - ddphi[2][1][1]) + 
        (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][1][2] + ddphi[1][1][1] + ddphi[0][0][1]) + 
        (dphi[0][1] - dphi[1][0])*(ddphi[0][1][1] - ddphi[1][0][1]) + 
        (dphi[0][2] - dphi[2][0])*(ddphi[0][1][2] - ddphi[2][0][1])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (16*phix*phiy*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
           phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (pow(dphi[0][1],2) + pow(dphi[1][1],2) + pow(dphi[2][1],2) + phix*ddphi[0][1][1] + phiy*ddphi[1][1][1] + 
           phiz*ddphi[2][1][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
           phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (24*phix*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
           phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (pow(dphi[0][1],2) + pow(dphi[1][1],2) + pow(dphi[2][1],2) + phix*ddphi[0][1][1] + phiy*ddphi[1][1][1] + 
           phiz*ddphi[2][1][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
           phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),4)) + 
   (8*(dphi[0][1] - dphi[1][0])*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + 
        phiy*phiz*(dphi[1][2] + dphi[2][1]) + pow(phix,2)*dphi[0][0] + 
        phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
           phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (pow(dphi[0][1],2) + pow(dphi[1][1],2) + pow(dphi[2][1],2) + phix*ddphi[0][1][1] + phiy*ddphi[1][1][1] + 
           phiz*ddphi[2][1][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
           phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (12*phix*phiy*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
      (2*phiy*pow(dphi[1][1],2) + 2*phiz*dphi[2][2]*dphi[2][1] + phiz*dphi[1][1]*(dphi[1][2] + dphi[2][1]) + 
        phiy*dphi[2][1]*(dphi[1][2] + dphi[2][1]) + pow(phiz,2)*ddphi[2][1][2] + pow(phiy,2)*ddphi[1][1][1] + 
        phiy*phiz*(ddphi[1][1][2] + ddphi[2][1][1]) + 2*phix*dphi[0][1]*dphi[0][0] + 
        dphi[0][1]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][1] + 
        phix*(dphi[1][1]*(dphi[0][1] + dphi[1][0]) + dphi[2][1]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][1][1] + ddphi[1][0][1]) + phiz*(ddphi[0][1][2] + ddphi[2][0][1]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (16*phix*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      (2*phiy*pow(dphi[1][1],2) + 2*phiz*dphi[2][2]*dphi[2][1] + phiz*dphi[1][1]*(dphi[1][2] + dphi[2][1]) + 
        phiy*dphi[2][1]*(dphi[1][2] + dphi[2][1]) + pow(phiz,2)*ddphi[2][1][2] + pow(phiy,2)*ddphi[1][1][1] + 
        phiy*phiz*(ddphi[1][1][2] + ddphi[2][1][1]) + 2*phix*dphi[0][1]*dphi[0][0] + 
        dphi[0][1]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][1] + 
        phix*(dphi[1][1]*(dphi[0][1] + dphi[1][0]) + dphi[2][1]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][1][1] + ddphi[1][0][1]) + phiz*(ddphi[0][1][2] + ddphi[2][0][1]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (8*(dphi[0][1] - dphi[1][0])*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + 
        phiy*phiz*(dphi[1][2] + dphi[2][1]) + pow(phix,2)*dphi[0][0] + 
        phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      (2*phiy*pow(dphi[1][1],2) + 2*phiz*dphi[2][2]*dphi[2][1] + phiz*dphi[1][1]*(dphi[1][2] + dphi[2][1]) + 
        phiy*dphi[2][1]*(dphi[1][2] + dphi[2][1]) + pow(phiz,2)*ddphi[2][1][2] + pow(phiy,2)*ddphi[1][1][1] + 
        phiy*phiz*(ddphi[1][1][2] + ddphi[2][1][1]) + 2*phix*dphi[0][1]*dphi[0][0] + 
        dphi[0][1]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][1] + 
        phix*(dphi[1][1]*(dphi[0][1] + dphi[1][0]) + dphi[2][1]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][1][1] + ddphi[1][0][1]) + phiz*(ddphi[0][1][2] + ddphi[2][0][1]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (2*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      (ddphi[2][0][2] + ddphi[1][0][1] + ddphi[0][0][0]))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (8*pow(phix,2)*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      ((dphi[1][2] - dphi[2][1])*(ddphi[1][0][2] - ddphi[2][0][1]) + 
        (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][0][2] + ddphi[1][0][1] + ddphi[0][0][0]) + 
        (dphi[0][1] - dphi[1][0])*(ddphi[0][0][1] - ddphi[1][0][0]) + 
        (dphi[0][2] - dphi[2][0])*(ddphi[0][0][2] - ddphi[2][0][0])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (8*phix*(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((dphi[1][2] - dphi[2][1])*(ddphi[1][0][2] - ddphi[2][0][1]) + 
        (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][0][2] + ddphi[1][0][1] + ddphi[0][0][0]) + 
        (dphi[0][1] - dphi[1][0])*(ddphi[0][0][1] - ddphi[1][0][0]) + 
        (dphi[0][2] - dphi[2][0])*(ddphi[0][0][2] - ddphi[2][0][0])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (4*phix*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      (pow(dphi[0][0],2) + pow(dphi[1][0],2) + pow(dphi[2][0],2) + phix*ddphi[0][0][0] + phiy*ddphi[1][0][0] + 
        phiz*ddphi[2][0][0]))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (16*pow(phix,2)*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
           phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
           phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (pow(dphi[0][0],2) + pow(dphi[1][0],2) + pow(dphi[2][0],2) + phix*ddphi[0][0][0] + phiy*ddphi[1][0][0] + 
           phiz*ddphi[2][0][0])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (24*phix*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
           phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
           phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (pow(dphi[0][0],2) + pow(dphi[1][0],2) + pow(dphi[2][0],2) + phix*ddphi[0][0][0] + phiy*ddphi[1][0][0] + 
           phiz*ddphi[2][0][0])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),4)) + 
   (8*(dphi[2][2] + dphi[1][1] + dphi[0][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
           phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
           phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (pow(dphi[0][0],2) + pow(dphi[1][0],2) + pow(dphi[2][0],2) + phix*ddphi[0][0][0] + phiy*ddphi[1][0][0] + 
           phiz*ddphi[2][0][0])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (12*pow(phix,2)*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
      (2*phix*pow(dphi[0][0],2) + 2*phiy*dphi[1][1]*dphi[1][0] + phiz*(dphi[1][2] + dphi[2][1])*dphi[1][0] + 
        2*phiz*dphi[2][2]*dphi[2][0] + phiy*(dphi[1][2] + dphi[2][1])*dphi[2][0] + 
        dphi[0][0]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phiz,2)*ddphi[2][0][2] + 
        pow(phiy,2)*ddphi[1][0][1] + phiy*phiz*(ddphi[1][0][2] + ddphi[2][0][1]) + pow(phix,2)*ddphi[0][0][0] + 
        phix*(dphi[1][0]*(dphi[0][1] + dphi[1][0]) + dphi[2][0]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][0][1] + ddphi[1][0][0]) + phiz*(ddphi[0][0][2] + ddphi[2][0][0]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (16*phix*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      (2*phix*pow(dphi[0][0],2) + 2*phiy*dphi[1][1]*dphi[1][0] + phiz*(dphi[1][2] + dphi[2][1])*dphi[1][0] + 
        2*phiz*dphi[2][2]*dphi[2][0] + phiy*(dphi[1][2] + dphi[2][1])*dphi[2][0] + 
        dphi[0][0]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phiz,2)*ddphi[2][0][2] + 
        pow(phiy,2)*ddphi[1][0][1] + phiy*phiz*(ddphi[1][0][2] + ddphi[2][0][1]) + pow(phix,2)*ddphi[0][0][0] + 
        phix*(dphi[1][0]*(dphi[0][1] + dphi[1][0]) + dphi[2][0]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][0][1] + ddphi[1][0][0]) + phiz*(ddphi[0][0][2] + ddphi[2][0][0]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (8*(dphi[2][2] + dphi[1][1] + dphi[0][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      (2*phix*pow(dphi[0][0],2) + 2*phiy*dphi[1][1]*dphi[1][0] + phiz*(dphi[1][2] + dphi[2][1])*dphi[1][0] + 
        2*phiz*dphi[2][2]*dphi[2][0] + phiy*(dphi[1][2] + dphi[2][1])*dphi[2][0] + 
        dphi[0][0]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phiz,2)*ddphi[2][0][2] + 
        pow(phiy,2)*ddphi[1][0][1] + phiy*phiz*(ddphi[1][0][2] + ddphi[2][0][1]) + pow(phix,2)*ddphi[0][0][0] + 
        phix*(dphi[1][0]*(dphi[0][1] + dphi[1][0]) + dphi[2][0]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][0][1] + ddphi[1][0][0]) + phiz*(ddphi[0][0][2] + ddphi[2][0][0]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2));

    return f4x;
}

double calc_f4y(double phix, double phiy, double phiz, double dphi[3][3], double ddphi[3][3][3])
{
	double f4y = (-4*phiz*dphi[1][2]*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (4*phiy*dphi[2][2]*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (16*phiy*phiz*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (8*phiy*dphi[1][1]*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (16*pow(phiy,2)*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (4*phiy*dphi[0][0]*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (4*phix*dphi[1][0]*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (4*(2*phiy*dphi[1][1] + phiz*(dphi[1][2] + dphi[2][1]) + phix*(dphi[0][1] + dphi[1][0]))*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (16*phix*phiy*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (4*dphi[1][2]*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (16*phiy*pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2)*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (4*phiy*(pow(dphi[0][2],2) + pow(dphi[1][2],2) + pow(dphi[2][2],2) + phix*ddphi[0][2][2] + phiy*ddphi[1][2][2] + 
        phiz*ddphi[2][2][2])*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (4*dphi[1][1]*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (16*phiy*pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2)*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (4*phiy*(pow(dphi[0][1],2) + pow(dphi[1][1],2) + pow(dphi[2][1],2) + phix*ddphi[0][1][1] + phiy*ddphi[1][1][1] + 
        phiz*ddphi[2][1][1])*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (4*dphi[1][0]*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (16*phiy*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2)*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (4*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      (dphi[1][2]*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2]) + 
        dphi[1][1]*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1]) + 
        dphi[1][0]*(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0]))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (8*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*(dphi[1][2] - dphi[2][1])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (2*(ddphi[1][2][2] - ddphi[2][1][2])*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + 
        phiy*phiz*(dphi[1][2] + dphi[2][1]) + pow(phix,2)*dphi[0][0] + 
        phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (8*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*(dphi[2][2] + dphi[1][1] + dphi[0][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (4*phiy*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (8*(dphi[0][1] - dphi[1][0])*(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (8*phiy*phiz*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      ((dphi[1][2] - dphi[2][1])*(ddphi[1][2][2] - ddphi[2][1][2]) + 
        (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][2][2] + ddphi[1][1][2] + ddphi[0][0][2]) + 
        (dphi[0][1] - dphi[1][0])*(ddphi[0][1][2] - ddphi[1][0][2]) + 
        (dphi[0][2] - dphi[2][0])*(ddphi[0][2][2] - ddphi[2][0][2])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (8*phiy*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((dphi[1][2] - dphi[2][1])*(ddphi[1][2][2] - ddphi[2][1][2]) + 
        (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][2][2] + ddphi[1][1][2] + ddphi[0][0][2]) + 
        (dphi[0][1] - dphi[1][0])*(ddphi[0][1][2] - ddphi[1][0][2]) + 
        (dphi[0][2] - dphi[2][0])*(ddphi[0][2][2] - ddphi[2][0][2])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (16*phiy*phiz*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (pow(dphi[0][2],2) + pow(dphi[1][2],2) + pow(dphi[2][2],2) + phix*ddphi[0][2][2] + phiy*ddphi[1][2][2] + 
           phiz*ddphi[2][2][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
           phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
           phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (24*phiy*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (pow(dphi[0][2],2) + pow(dphi[1][2],2) + pow(dphi[2][2],2) + phix*ddphi[0][2][2] + phiy*ddphi[1][2][2] + 
           phiz*ddphi[2][2][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
           phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
           phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),4)) + 
   (8*(dphi[1][2] - dphi[2][1])*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + 
        phiy*phiz*(dphi[1][2] + dphi[2][1]) + pow(phix,2)*dphi[0][0] + 
        phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (pow(dphi[0][2],2) + pow(dphi[1][2],2) + pow(dphi[2][2],2) + phix*ddphi[0][2][2] + phiy*ddphi[1][2][2] + 
           phiz*ddphi[2][2][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
           phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
           phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (12*phiy*phiz*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
      (2*phiz*pow(dphi[2][2],2) + pow(phiz,2)*ddphi[2][2][2] + 2*phiy*dphi[1][2]*dphi[1][1] + 
        phiz*dphi[1][2]*(dphi[1][2] + dphi[2][1]) + phiy*dphi[2][2]*(dphi[1][2] + dphi[2][1]) + 
        pow(phiy,2)*ddphi[1][1][2] + phiy*phiz*(ddphi[1][2][2] + ddphi[2][1][2]) + 2*phix*dphi[0][2]*dphi[0][0] + 
        dphi[0][2]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][2] + 
        phix*(dphi[1][2]*(dphi[0][1] + dphi[1][0]) + dphi[2][2]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][1][2] + ddphi[1][0][2]) + phiz*(ddphi[0][2][2] + ddphi[2][0][2]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (16*phiy*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      (2*phiz*pow(dphi[2][2],2) + pow(phiz,2)*ddphi[2][2][2] + 2*phiy*dphi[1][2]*dphi[1][1] + 
        phiz*dphi[1][2]*(dphi[1][2] + dphi[2][1]) + phiy*dphi[2][2]*(dphi[1][2] + dphi[2][1]) + 
        pow(phiy,2)*ddphi[1][1][2] + phiy*phiz*(ddphi[1][2][2] + ddphi[2][1][2]) + 2*phix*dphi[0][2]*dphi[0][0] + 
        dphi[0][2]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][2] + 
        phix*(dphi[1][2]*(dphi[0][1] + dphi[1][0]) + dphi[2][2]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][1][2] + ddphi[1][0][2]) + phiz*(ddphi[0][2][2] + ddphi[2][0][2]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (8*(dphi[1][2] - dphi[2][1])*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + 
        phiy*phiz*(dphi[1][2] + dphi[2][1]) + pow(phix,2)*dphi[0][0] + 
        phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      (2*phiz*pow(dphi[2][2],2) + pow(phiz,2)*ddphi[2][2][2] + 2*phiy*dphi[1][2]*dphi[1][1] + 
        phiz*dphi[1][2]*(dphi[1][2] + dphi[2][1]) + phiy*dphi[2][2]*(dphi[1][2] + dphi[2][1]) + 
        pow(phiy,2)*ddphi[1][1][2] + phiy*phiz*(ddphi[1][2][2] + ddphi[2][1][2]) + 2*phix*dphi[0][2]*dphi[0][0] + 
        dphi[0][2]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][2] + 
        phix*(dphi[1][2]*(dphi[0][1] + dphi[1][0]) + dphi[2][2]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][1][2] + ddphi[1][0][2]) + phiz*(ddphi[0][2][2] + ddphi[2][0][2]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (2*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      (ddphi[2][1][2] + ddphi[1][1][1] + ddphi[0][0][1]))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (8*pow(phiy,2)*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      ((dphi[1][2] - dphi[2][1])*(ddphi[1][1][2] - ddphi[2][1][1]) + 
        (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][1][2] + ddphi[1][1][1] + ddphi[0][0][1]) + 
        (dphi[0][1] - dphi[1][0])*(ddphi[0][1][1] - ddphi[1][0][1]) + 
        (dphi[0][2] - dphi[2][0])*(ddphi[0][1][2] - ddphi[2][0][1])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (8*phiy*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((dphi[1][2] - dphi[2][1])*(ddphi[1][1][2] - ddphi[2][1][1]) + 
        (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][1][2] + ddphi[1][1][1] + ddphi[0][0][1]) + 
        (dphi[0][1] - dphi[1][0])*(ddphi[0][1][1] - ddphi[1][0][1]) + 
        (dphi[0][2] - dphi[2][0])*(ddphi[0][1][2] - ddphi[2][0][1])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (16*pow(phiy,2)*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
           phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (pow(dphi[0][1],2) + pow(dphi[1][1],2) + pow(dphi[2][1],2) + phix*ddphi[0][1][1] + phiy*ddphi[1][1][1] + 
           phiz*ddphi[2][1][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
           phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (24*phiy*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
           phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (pow(dphi[0][1],2) + pow(dphi[1][1],2) + pow(dphi[2][1],2) + phix*ddphi[0][1][1] + phiy*ddphi[1][1][1] + 
           phiz*ddphi[2][1][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
           phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),4)) + 
   (8*(dphi[2][2] + dphi[1][1] + dphi[0][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
           phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (pow(dphi[0][1],2) + pow(dphi[1][1],2) + pow(dphi[2][1],2) + phix*ddphi[0][1][1] + phiy*ddphi[1][1][1] + 
           phiz*ddphi[2][1][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
           phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (12*pow(phiy,2)*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
      (2*phiy*pow(dphi[1][1],2) + 2*phiz*dphi[2][2]*dphi[2][1] + phiz*dphi[1][1]*(dphi[1][2] + dphi[2][1]) + 
        phiy*dphi[2][1]*(dphi[1][2] + dphi[2][1]) + pow(phiz,2)*ddphi[2][1][2] + pow(phiy,2)*ddphi[1][1][1] + 
        phiy*phiz*(ddphi[1][1][2] + ddphi[2][1][1]) + 2*phix*dphi[0][1]*dphi[0][0] + 
        dphi[0][1]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][1] + 
        phix*(dphi[1][1]*(dphi[0][1] + dphi[1][0]) + dphi[2][1]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][1][1] + ddphi[1][0][1]) + phiz*(ddphi[0][1][2] + ddphi[2][0][1]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (16*phiy*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      (2*phiy*pow(dphi[1][1],2) + 2*phiz*dphi[2][2]*dphi[2][1] + phiz*dphi[1][1]*(dphi[1][2] + dphi[2][1]) + 
        phiy*dphi[2][1]*(dphi[1][2] + dphi[2][1]) + pow(phiz,2)*ddphi[2][1][2] + pow(phiy,2)*ddphi[1][1][1] + 
        phiy*phiz*(ddphi[1][1][2] + ddphi[2][1][1]) + 2*phix*dphi[0][1]*dphi[0][0] + 
        dphi[0][1]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][1] + 
        phix*(dphi[1][1]*(dphi[0][1] + dphi[1][0]) + dphi[2][1]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][1][1] + ddphi[1][0][1]) + phiz*(ddphi[0][1][2] + ddphi[2][0][1]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (8*(dphi[2][2] + dphi[1][1] + dphi[0][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      (2*phiy*pow(dphi[1][1],2) + 2*phiz*dphi[2][2]*dphi[2][1] + phiz*dphi[1][1]*(dphi[1][2] + dphi[2][1]) + 
        phiy*dphi[2][1]*(dphi[1][2] + dphi[2][1]) + pow(phiz,2)*ddphi[2][1][2] + pow(phiy,2)*ddphi[1][1][1] + 
        phiy*phiz*(ddphi[1][1][2] + ddphi[2][1][1]) + 2*phix*dphi[0][1]*dphi[0][0] + 
        dphi[0][1]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][1] + 
        phix*(dphi[1][1]*(dphi[0][1] + dphi[1][0]) + dphi[2][1]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][1][1] + ddphi[1][0][1]) + phiz*(ddphi[0][1][2] + ddphi[2][0][1]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (2*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      (ddphi[0][0][1] - ddphi[1][0][0]))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (8*phix*phiy*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      ((dphi[1][2] - dphi[2][1])*(ddphi[1][0][2] - ddphi[2][0][1]) + 
        (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][0][2] + ddphi[1][0][1] + ddphi[0][0][0]) + 
        (dphi[0][1] - dphi[1][0])*(ddphi[0][0][1] - ddphi[1][0][0]) + 
        (dphi[0][2] - dphi[2][0])*(ddphi[0][0][2] - ddphi[2][0][0])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (8*phiy*(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((dphi[1][2] - dphi[2][1])*(ddphi[1][0][2] - ddphi[2][0][1]) + 
        (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][0][2] + ddphi[1][0][1] + ddphi[0][0][0]) + 
        (dphi[0][1] - dphi[1][0])*(ddphi[0][0][1] - ddphi[1][0][0]) + 
        (dphi[0][2] - dphi[2][0])*(ddphi[0][0][2] - ddphi[2][0][0])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (4*phiy*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      (pow(dphi[0][0],2) + pow(dphi[1][0],2) + pow(dphi[2][0],2) + phix*ddphi[0][0][0] + phiy*ddphi[1][0][0] + 
        phiz*ddphi[2][0][0]))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (16*phix*phiy*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
           phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
           phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (pow(dphi[0][0],2) + pow(dphi[1][0],2) + pow(dphi[2][0],2) + phix*ddphi[0][0][0] + phiy*ddphi[1][0][0] + 
           phiz*ddphi[2][0][0])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (24*phiy*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
           phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
           phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (pow(dphi[0][0],2) + pow(dphi[1][0],2) + pow(dphi[2][0],2) + phix*ddphi[0][0][0] + phiy*ddphi[1][0][0] + 
           phiz*ddphi[2][0][0])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),4)) - 
   (8*(dphi[0][1] - dphi[1][0])*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + 
        phiy*phiz*(dphi[1][2] + dphi[2][1]) + pow(phix,2)*dphi[0][0] + 
        phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
           phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
           phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (pow(dphi[0][0],2) + pow(dphi[1][0],2) + pow(dphi[2][0],2) + phix*ddphi[0][0][0] + phiy*ddphi[1][0][0] + 
           phiz*ddphi[2][0][0])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (12*phix*phiy*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
      (2*phix*pow(dphi[0][0],2) + 2*phiy*dphi[1][1]*dphi[1][0] + phiz*(dphi[1][2] + dphi[2][1])*dphi[1][0] + 
        2*phiz*dphi[2][2]*dphi[2][0] + phiy*(dphi[1][2] + dphi[2][1])*dphi[2][0] + 
        dphi[0][0]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phiz,2)*ddphi[2][0][2] + 
        pow(phiy,2)*ddphi[1][0][1] + phiy*phiz*(ddphi[1][0][2] + ddphi[2][0][1]) + pow(phix,2)*ddphi[0][0][0] + 
        phix*(dphi[1][0]*(dphi[0][1] + dphi[1][0]) + dphi[2][0]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][0][1] + ddphi[1][0][0]) + phiz*(ddphi[0][0][2] + ddphi[2][0][0]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (16*phiy*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      (2*phix*pow(dphi[0][0],2) + 2*phiy*dphi[1][1]*dphi[1][0] + phiz*(dphi[1][2] + dphi[2][1])*dphi[1][0] + 
        2*phiz*dphi[2][2]*dphi[2][0] + phiy*(dphi[1][2] + dphi[2][1])*dphi[2][0] + 
        dphi[0][0]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phiz,2)*ddphi[2][0][2] + 
        pow(phiy,2)*ddphi[1][0][1] + phiy*phiz*(ddphi[1][0][2] + ddphi[2][0][1]) + pow(phix,2)*ddphi[0][0][0] + 
        phix*(dphi[1][0]*(dphi[0][1] + dphi[1][0]) + dphi[2][0]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][0][1] + ddphi[1][0][0]) + phiz*(ddphi[0][0][2] + ddphi[2][0][0]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (8*(dphi[0][1] - dphi[1][0])*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + 
        phiy*phiz*(dphi[1][2] + dphi[2][1]) + pow(phix,2)*dphi[0][0] + 
        phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      (2*phix*pow(dphi[0][0],2) + 2*phiy*dphi[1][1]*dphi[1][0] + phiz*(dphi[1][2] + dphi[2][1])*dphi[1][0] + 
        2*phiz*dphi[2][2]*dphi[2][0] + phiy*(dphi[1][2] + dphi[2][1])*dphi[2][0] + 
        dphi[0][0]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phiz,2)*ddphi[2][0][2] + 
        pow(phiy,2)*ddphi[1][0][1] + phiy*phiz*(ddphi[1][0][2] + ddphi[2][0][1]) + pow(phix,2)*ddphi[0][0][0] + 
        phix*(dphi[1][0]*(dphi[0][1] + dphi[1][0]) + dphi[2][0]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][0][1] + ddphi[1][0][0]) + phiz*(ddphi[0][0][2] + ddphi[2][0][0]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2));

    return f4y;
}

double calc_f4z(double phix, double phiy, double phiz, double dphi[3][3], double ddphi[3][3][3])
{
	double f4z = (-8*phiz*dphi[2][2]*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (16*pow(phiz,2)*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (4*phiz*dphi[1][1]*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (4*phiy*dphi[2][1]*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (16*phiy*phiz*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (4*phiz*dphi[0][0]*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (4*phix*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*dphi[2][0]*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (16*phix*phiz*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (4*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      (2*phiz*dphi[2][2] + phiy*(dphi[1][2] + dphi[2][1]) + phix*(dphi[0][2] + dphi[2][0]))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (4*dphi[2][2]*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (16*phiz*pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2)*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (4*phiz*(pow(dphi[0][2],2) + pow(dphi[1][2],2) + pow(dphi[2][2],2) + phix*ddphi[0][2][2] + phiy*ddphi[1][2][2] + 
        phiz*ddphi[2][2][2])*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (4*dphi[2][1]*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (16*phiz*pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2)*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (4*phiz*(pow(dphi[0][1],2) + pow(dphi[1][1],2) + pow(dphi[2][1],2) + phix*ddphi[0][1][1] + phiy*ddphi[1][1][1] + 
        phiz*ddphi[2][1][1])*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (4*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*dphi[2][0]*
      (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (16*phiz*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2)*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (4*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      (dphi[2][2]*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2]) + 
        dphi[2][1]*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1]) + 
        dphi[2][0]*(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0]))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (8*(dphi[1][2] - dphi[2][1])*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (2*(ddphi[1][1][2] - ddphi[2][1][1])*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + 
        phiy*phiz*(dphi[1][2] + dphi[2][1]) + pow(phix,2)*dphi[0][0] + 
        phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (8*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*(dphi[2][2] + dphi[1][1] + dphi[0][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (4*phiz*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (8*(dphi[0][2] - dphi[2][0])*(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),3)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (2*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      (ddphi[2][2][2] + ddphi[1][1][2] + ddphi[0][0][2]))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (8*pow(phiz,2)*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      ((dphi[1][2] - dphi[2][1])*(ddphi[1][2][2] - ddphi[2][1][2]) + 
        (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][2][2] + ddphi[1][1][2] + ddphi[0][0][2]) + 
        (dphi[0][1] - dphi[1][0])*(ddphi[0][1][2] - ddphi[1][0][2]) + 
        (dphi[0][2] - dphi[2][0])*(ddphi[0][2][2] - ddphi[2][0][2])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (8*phiz*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((dphi[1][2] - dphi[2][1])*(ddphi[1][2][2] - ddphi[2][1][2]) + 
        (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][2][2] + ddphi[1][1][2] + ddphi[0][0][2]) + 
        (dphi[0][1] - dphi[1][0])*(ddphi[0][1][2] - ddphi[1][0][2]) + 
        (dphi[0][2] - dphi[2][0])*(ddphi[0][2][2] - ddphi[2][0][2])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (16*pow(phiz,2)*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (pow(dphi[0][2],2) + pow(dphi[1][2],2) + pow(dphi[2][2],2) + phix*ddphi[0][2][2] + phiy*ddphi[1][2][2] + 
           phiz*ddphi[2][2][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
           phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
           phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (24*phiz*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (pow(dphi[0][2],2) + pow(dphi[1][2],2) + pow(dphi[2][2],2) + phix*ddphi[0][2][2] + phiy*ddphi[1][2][2] + 
           phiz*ddphi[2][2][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
           phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
           phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),4)) + 
   (8*(dphi[2][2] + dphi[1][1] + dphi[0][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (pow(dphi[0][2],2) + pow(dphi[1][2],2) + pow(dphi[2][2],2) + phix*ddphi[0][2][2] + phiy*ddphi[1][2][2] + 
           phiz*ddphi[2][2][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
           phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
           phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (12*pow(phiz,2)*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
      (2*phiz*pow(dphi[2][2],2) + pow(phiz,2)*ddphi[2][2][2] + 2*phiy*dphi[1][2]*dphi[1][1] + 
        phiz*dphi[1][2]*(dphi[1][2] + dphi[2][1]) + phiy*dphi[2][2]*(dphi[1][2] + dphi[2][1]) + 
        pow(phiy,2)*ddphi[1][1][2] + phiy*phiz*(ddphi[1][2][2] + ddphi[2][1][2]) + 2*phix*dphi[0][2]*dphi[0][0] + 
        dphi[0][2]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][2] + 
        phix*(dphi[1][2]*(dphi[0][1] + dphi[1][0]) + dphi[2][2]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][1][2] + ddphi[1][0][2]) + phiz*(ddphi[0][2][2] + ddphi[2][0][2]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (16*phiz*(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      (2*phiz*pow(dphi[2][2],2) + pow(phiz,2)*ddphi[2][2][2] + 2*phiy*dphi[1][2]*dphi[1][1] + 
        phiz*dphi[1][2]*(dphi[1][2] + dphi[2][1]) + phiy*dphi[2][2]*(dphi[1][2] + dphi[2][1]) + 
        pow(phiy,2)*ddphi[1][1][2] + phiy*phiz*(ddphi[1][2][2] + ddphi[2][1][2]) + 2*phix*dphi[0][2]*dphi[0][0] + 
        dphi[0][2]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][2] + 
        phix*(dphi[1][2]*(dphi[0][1] + dphi[1][0]) + dphi[2][2]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][1][2] + ddphi[1][0][2]) + phiz*(ddphi[0][2][2] + ddphi[2][0][2]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (8*(dphi[2][2] + dphi[1][1] + dphi[0][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      (2*phiz*pow(dphi[2][2],2) + pow(phiz,2)*ddphi[2][2][2] + 2*phiy*dphi[1][2]*dphi[1][1] + 
        phiz*dphi[1][2]*(dphi[1][2] + dphi[2][1]) + phiy*dphi[2][2]*(dphi[1][2] + dphi[2][1]) + 
        pow(phiy,2)*ddphi[1][1][2] + phiy*phiz*(ddphi[1][2][2] + ddphi[2][1][2]) + 2*phix*dphi[0][2]*dphi[0][0] + 
        dphi[0][2]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][2] + 
        phix*(dphi[1][2]*(dphi[0][1] + dphi[1][0]) + dphi[2][2]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][1][2] + ddphi[1][0][2]) + phiz*(ddphi[0][2][2] + ddphi[2][0][2]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (8*phiy*phiz*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      ((dphi[1][2] - dphi[2][1])*(ddphi[1][1][2] - ddphi[2][1][1]) + 
        (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][1][2] + ddphi[1][1][1] + ddphi[0][0][1]) + 
        (dphi[0][1] - dphi[1][0])*(ddphi[0][1][1] - ddphi[1][0][1]) + 
        (dphi[0][2] - dphi[2][0])*(ddphi[0][1][2] - ddphi[2][0][1])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (8*phiz*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((dphi[1][2] - dphi[2][1])*(ddphi[1][1][2] - ddphi[2][1][1]) + 
        (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][1][2] + ddphi[1][1][1] + ddphi[0][0][1]) + 
        (dphi[0][1] - dphi[1][0])*(ddphi[0][1][1] - ddphi[1][0][1]) + 
        (dphi[0][2] - dphi[2][0])*(ddphi[0][1][2] - ddphi[2][0][1])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (16*phiy*phiz*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
           phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (pow(dphi[0][1],2) + pow(dphi[1][1],2) + pow(dphi[2][1],2) + phix*ddphi[0][1][1] + phiy*ddphi[1][1][1] + 
           phiz*ddphi[2][1][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
           phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (24*phiz*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
           phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (pow(dphi[0][1],2) + pow(dphi[1][1],2) + pow(dphi[2][1],2) + phix*ddphi[0][1][1] + phiy*ddphi[1][1][1] + 
           phiz*ddphi[2][1][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
           phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),4)) - 
   (8*(dphi[1][2] - dphi[2][1])*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + 
        phiy*phiz*(dphi[1][2] + dphi[2][1]) + pow(phix,2)*dphi[0][0] + 
        phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (dphi[0][2]*dphi[0][1] + dphi[1][2]*dphi[1][1] + dphi[2][2]*dphi[2][1] + phix*ddphi[0][1][2] + 
           phiy*ddphi[1][1][2] + phiz*ddphi[2][1][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (pow(dphi[0][1],2) + pow(dphi[1][1],2) + pow(dphi[2][1],2) + phix*ddphi[0][1][1] + phiy*ddphi[1][1][1] + 
           phiz*ddphi[2][1][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
           phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (12*phiy*phiz*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
      (2*phiy*pow(dphi[1][1],2) + 2*phiz*dphi[2][2]*dphi[2][1] + phiz*dphi[1][1]*(dphi[1][2] + dphi[2][1]) + 
        phiy*dphi[2][1]*(dphi[1][2] + dphi[2][1]) + pow(phiz,2)*ddphi[2][1][2] + pow(phiy,2)*ddphi[1][1][1] + 
        phiy*phiz*(ddphi[1][1][2] + ddphi[2][1][1]) + 2*phix*dphi[0][1]*dphi[0][0] + 
        dphi[0][1]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][1] + 
        phix*(dphi[1][1]*(dphi[0][1] + dphi[1][0]) + dphi[2][1]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][1][1] + ddphi[1][0][1]) + phiz*(ddphi[0][1][2] + ddphi[2][0][1]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (16*phiz*(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
      (pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      (2*phiy*pow(dphi[1][1],2) + 2*phiz*dphi[2][2]*dphi[2][1] + phiz*dphi[1][1]*(dphi[1][2] + dphi[2][1]) + 
        phiy*dphi[2][1]*(dphi[1][2] + dphi[2][1]) + pow(phiz,2)*ddphi[2][1][2] + pow(phiy,2)*ddphi[1][1][1] + 
        phiy*phiz*(ddphi[1][1][2] + ddphi[2][1][1]) + 2*phix*dphi[0][1]*dphi[0][0] + 
        dphi[0][1]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][1] + 
        phix*(dphi[1][1]*(dphi[0][1] + dphi[1][0]) + dphi[2][1]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][1][1] + ddphi[1][0][1]) + phiz*(ddphi[0][1][2] + ddphi[2][0][1]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (8*(dphi[1][2] - dphi[2][1])*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + 
        phiy*phiz*(dphi[1][2] + dphi[2][1]) + pow(phix,2)*dphi[0][0] + 
        phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      (2*phiy*pow(dphi[1][1],2) + 2*phiz*dphi[2][2]*dphi[2][1] + phiz*dphi[1][1]*(dphi[1][2] + dphi[2][1]) + 
        phiy*dphi[2][1]*(dphi[1][2] + dphi[2][1]) + pow(phiz,2)*ddphi[2][1][2] + pow(phiy,2)*ddphi[1][1][1] + 
        phiy*phiz*(ddphi[1][1][2] + ddphi[2][1][1]) + 2*phix*dphi[0][1]*dphi[0][0] + 
        dphi[0][1]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phix,2)*ddphi[0][0][1] + 
        phix*(dphi[1][1]*(dphi[0][1] + dphi[1][0]) + dphi[2][1]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][1][1] + ddphi[1][0][1]) + phiz*(ddphi[0][1][2] + ddphi[2][0][1]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) - 
   (8*phix*phiz*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      ((dphi[1][2] - dphi[2][1])*(ddphi[1][0][2] - ddphi[2][0][1]) + 
        (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][0][2] + ddphi[1][0][1] + ddphi[0][0][0]) + 
        (dphi[0][1] - dphi[1][0])*(ddphi[0][0][1] - ddphi[1][0][0]) + 
        (dphi[0][2] - dphi[2][0])*(ddphi[0][0][2] - ddphi[2][0][0])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (8*phiz*(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((dphi[1][2] - dphi[2][1])*(ddphi[1][0][2] - ddphi[2][0][1]) + 
        (dphi[2][2] + dphi[1][1] + dphi[0][0])*(ddphi[2][0][2] + ddphi[1][0][1] + ddphi[0][0][0]) + 
        (dphi[0][1] - dphi[1][0])*(ddphi[0][0][1] - ddphi[1][0][0]) + 
        (dphi[0][2] - dphi[2][0])*(ddphi[0][0][2] - ddphi[2][0][0])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (2*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      (ddphi[0][0][2] - ddphi[2][0][0]))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (4*phiz*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      (pow(dphi[0][0],2) + pow(dphi[1][0],2) + pow(dphi[2][0],2) + phix*ddphi[0][0][0] + phiy*ddphi[1][0][0] + 
        phiz*ddphi[2][0][0]))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (16*phix*phiz*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
           phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
           phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (pow(dphi[0][0],2) + pow(dphi[1][0],2) + pow(dphi[2][0],2) + phix*ddphi[0][0][0] + phiy*ddphi[1][0][0] + 
           phiz*ddphi[2][0][0])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (24*phiz*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
           phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
           phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (pow(dphi[0][0],2) + pow(dphi[1][0],2) + pow(dphi[2][0],2) + phix*ddphi[0][0][0] + phiy*ddphi[1][0][0] + 
           phiz*ddphi[2][0][0])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),4)) - 
   (8*(dphi[0][2] - dphi[2][0])*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + 
        phiy*phiz*(dphi[1][2] + dphi[2][1]) + pow(phix,2)*dphi[0][0] + 
        phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),4)*
      ((phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2])*
         (dphi[0][2]*dphi[0][0] + dphi[1][2]*dphi[1][0] + dphi[2][2]*dphi[2][0] + phix*ddphi[0][0][2] + 
           phiy*ddphi[1][0][2] + phiz*ddphi[2][0][2]) + (phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1])*
         (dphi[0][1]*dphi[0][0] + dphi[1][1]*dphi[1][0] + dphi[2][1]*dphi[2][0] + phix*ddphi[0][0][1] + 
           phiy*ddphi[1][0][1] + phiz*ddphi[2][0][1]) + (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
         (pow(dphi[0][0],2) + pow(dphi[1][0],2) + pow(dphi[2][0],2) + phix*ddphi[0][0][0] + phiy*ddphi[1][0][0] + 
           phiz*ddphi[2][0][0])))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) - 
   (12*phix*phiz*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),2)*
      (2*phix*pow(dphi[0][0],2) + 2*phiy*dphi[1][1]*dphi[1][0] + phiz*(dphi[1][2] + dphi[2][1])*dphi[1][0] + 
        2*phiz*dphi[2][2]*dphi[2][0] + phiy*(dphi[1][2] + dphi[2][1])*dphi[2][0] + 
        dphi[0][0]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phiz,2)*ddphi[2][0][2] + 
        pow(phiy,2)*ddphi[1][0][1] + phiy*phiz*(ddphi[1][0][2] + ddphi[2][0][1]) + pow(phix,2)*ddphi[0][0][0] + 
        phix*(dphi[1][0]*(dphi[0][1] + dphi[1][0]) + dphi[2][0]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][0][1] + ddphi[1][0][0]) + phiz*(ddphi[0][0][2] + ddphi[2][0][0]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2)) + 
   (16*phiz*(pow(dphi[1][2] - dphi[2][1],2) + pow(dphi[2][2] + dphi[1][1] + dphi[0][0],2) + 
        pow(dphi[0][1] - dphi[1][0],2) + pow(dphi[0][2] - dphi[2][0],2))*
      (phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0])*
      pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + phiy*phiz*(dphi[1][2] + dphi[2][1]) + 
        pow(phix,2)*dphi[0][0] + phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      (2*phix*pow(dphi[0][0],2) + 2*phiy*dphi[1][1]*dphi[1][0] + phiz*(dphi[1][2] + dphi[2][1])*dphi[1][0] + 
        2*phiz*dphi[2][2]*dphi[2][0] + phiy*(dphi[1][2] + dphi[2][1])*dphi[2][0] + 
        dphi[0][0]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phiz,2)*ddphi[2][0][2] + 
        pow(phiy,2)*ddphi[1][0][1] + phiy*phiz*(ddphi[1][0][2] + ddphi[2][0][1]) + pow(phix,2)*ddphi[0][0][0] + 
        phix*(dphi[1][0]*(dphi[0][1] + dphi[1][0]) + dphi[2][0]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][0][1] + ddphi[1][0][0]) + phiz*(ddphi[0][0][2] + ddphi[2][0][0]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),3)) + 
   (8*(dphi[0][2] - dphi[2][0])*pow(pow(phiz,2)*dphi[2][2] + pow(phiy,2)*dphi[1][1] + 
        phiy*phiz*(dphi[1][2] + dphi[2][1]) + pow(phix,2)*dphi[0][0] + 
        phix*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])),3)*
      (2*phix*pow(dphi[0][0],2) + 2*phiy*dphi[1][1]*dphi[1][0] + phiz*(dphi[1][2] + dphi[2][1])*dphi[1][0] + 
        2*phiz*dphi[2][2]*dphi[2][0] + phiy*(dphi[1][2] + dphi[2][1])*dphi[2][0] + 
        dphi[0][0]*(phiy*(dphi[0][1] + dphi[1][0]) + phiz*(dphi[0][2] + dphi[2][0])) + pow(phiz,2)*ddphi[2][0][2] + 
        pow(phiy,2)*ddphi[1][0][1] + phiy*phiz*(ddphi[1][0][2] + ddphi[2][0][1]) + pow(phix,2)*ddphi[0][0][0] + 
        phix*(dphi[1][0]*(dphi[0][1] + dphi[1][0]) + dphi[2][0]*(dphi[0][2] + dphi[2][0]) + 
           phiy*(ddphi[0][0][1] + ddphi[1][0][0]) + phiz*(ddphi[0][0][2] + ddphi[2][0][0]))))/
    (pow(pow(phix,2) + pow(phiy,2) + pow(phiz,2),2)*pow(pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
        pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
        pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2),2));

	return f4z;
}
