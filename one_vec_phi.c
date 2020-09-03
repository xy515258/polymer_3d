/*
 * Yang Xia, 3d version, parallel
 * This is now 2d due to the anisotropic interfacial energy
 */
#include <mpi.h>
#include <complex.h>
#include <fftw3-mpi.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include "aniso_interfacial.h"

/* Simulation size parameters.
 * Since N can be large, dynamically allocate arrays to prevent stack overflow.*/
#define NX 100
#define NY 100
#define NZ 100
#define Nt 200000

/* These parameters may not have been reduced to the most minimal set of parameters
 * possible to describe the system. */
#define DELTA_T 3.e-7
#define mobility1_max 1.e4          /* polymorph 1 */
#define DELTA_X 0.01
#define DELTA_Y DELTA_X             /* Assume for now DELTA_X = DELTA_Y = DELTA_Z. */
#define DELTA_Z DELTA_X
#define alpha_1 0.9

/* Interfacial */
#define epsc 2.e-5
#define eps_ani_2 0.
#define eps_ani_4 0.
#define eps_higher_order 0.e-9
#define eps_cs -2.5e-5
#define eps_cm -0.e-5
#define eps_m 4.e-5
#define eps_ms -3.e-4
#define gamma_cs -0.
#define gamma_cm -0.

/* Rotation angle anisotropy. */
#define deltabeta 2
#define nbeta 4

/* Anisotropy. */
#define deltaM 200.

/* u-related parameters */ 
#define D_U 20.
#define ku 0.25
#define k1_phi 0.01
#define k1_zeta 0.1

/* Miscellaneous parameters to tune. */
#define n 24
#define neps 2.
#define b1 -0.0025
#define b2 -0.01
#define gamma 10.
#define INIT_U -1. /* Actually I think this always has to be -1. */
#define INIT_THETA (0.*M_PI/180.)
#define INIT_PHI (90.*M_PI/180.) /* To recover 2d case, set this to pi/2. */
#define INIT_PHI_AMP 0.00
#define kBT (0.e-7)
#define cbar1 0.95
#define ctilde1 0.01
#define sub_int 20.
#define rm 20.
#define theta_0 M_PI/6 /* The angle of epitaxial orientation, range in -60 to 60 degree */
#define k2 0.2

/* Mobility in the z direction. */
#define m_z1 (-0.0 + 1.5*b1/(1.+b1)) /* positive values favor edge on (phi=90). */
#define n_phi 4. 

/* corr_length is relative to a number scaled from 0 to 1. */
#define corr_length_sq (1.e-4)

#define N_TO_SAVE 10 /* Number of snapshots during the process to save. */
/* This should divide evenly into Nt. */
#define when_to_save ((int) (((double) Nt)/((double) N_TO_SAVE)))

/* Define modulo to handle negative numbers properly. */
#define mod(a,b) ((((a)%(b))+(b))%(b))               //handle with periodic boundary
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )

double*** Make3DDoubleArray(int arraySizeX, int arraySizeY, int arraySizeZ);
void write_to_file(char *filename, double ***data, int local_NX, int local_0_start);       
void read_from_file(char *filename, double ***data, int local_NX, int local_0_start);   //Still 2D
double gaussrand(void);
void generate_noise_matrices(double ****eta, fftw_complex* sqrt_correlation_func_q, int local_NX); 
void Fill_boundary(double *** data, int local_NX);
void Comunicate_boundary(double *** data, int local_NX);
void initialize(double **** phi, double *** u, double *** zeta, double *** corr_func_data, int local_NX, int local_0_start);
void initialize_u(double * temp_u);
void initialize_zeta(double * temp_zeta);
int main(int argc, char **argv)
{
    int nproc, myrank;
    /* Initialize MPI */
    MPI_Init(&argc,&argv);
    fftw_mpi_init();
    /* Get the number of processes */
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    /* Get my process number (rank) */
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    
    ptrdiff_t local_NX;
    ptrdiff_t local_0_start;
    ptrdiff_t alloc_local = fftw_mpi_local_size_3d(NX, NY, NZ, MPI_COMM_WORLD, &local_NX, &local_0_start);
    
    double **** fphi = (double****) malloc(3*sizeof(double***));
    double **** phi = (double****) malloc(3*sizeof(double***));
    double **** eta = (double****) malloc(3*sizeof(double***));
    double **** Laplacian = (double****) malloc(3*sizeof(double***));
    for( int p=0; p<3; p++)
    {
        fphi[p] = Make3DDoubleArray(local_NX, NY, NZ);          //Change of phi in one step
        phi[p] = Make3DDoubleArray(local_NX, NY, NZ);
        eta[p] = Make3DDoubleArray(local_NX, NY, NZ);           //Fluctuation  
        Laplacian[p] = Make3DDoubleArray(local_NX, NY, NZ);
    }
    double*** u = Make3DDoubleArray(local_NX, NY, NZ);
    double*** fu = Make3DDoubleArray(local_NX, NY, NZ);              //Change of u in one step
    double*** corr_func_data = Make3DDoubleArray(local_NX, NY, NZ);
    double*** zeta = Make3DDoubleArray(local_NX, NY, NZ);            //Order parameter for substrate
    double*** dFdu = Make3DDoubleArray(local_NX, NY, NZ);            //delta F/ delta u
    double*** Laplacian_u = Make3DDoubleArray(local_NX, NY, NZ);     //nabla^2 u
    double*** Laplacian_dFdu = Make3DDoubleArray(local_NX, NY, NZ);  //nabla^2 dFdu
    double factors, morefactors;
    double absphi, absphisq1;
    double theta, phi_s;
    double Mphi;
    double m_value;
    double fluct_amp;
    double newphi1, newphi2, newphi3, dabsphi;
    double absphixy, absphixyinv;
    double dphix, dphiy, dphiz;
    double last_terms_factor, first_term, second_term_factor, second_term1, second_term2;
    double cosphi, cosphisq, sinphi, sinphisq;
    double b1factor, b2factor, hex_factor, hex_factor_p;
    double l1, l2, l3;
    double deltatheta;
    double absphipx, absphipy, absphimx, absphimy, absphipz, absphimz;
    double dphi_r_x, dphi_r_y, dphi_r_z, dphi_r_norm, v;
    double isotropic_mobility;
    double c1, DELTA;
    double phix, phiy, phiz;
    double dphi[3][3];
    double ddphi[3][3][3];
    double f2x, f4x, f2y, f4y, f2z, f4z;
    double dzeta_x, dzeta_y, dzeta_z, grad_zeta_sq;
    double grad_u_sq;
    double du[3];
    double ddu[3][3];
    double f_cs1, f_cs2, f_cm1, f_cm2, f_mc1, f_mc2, f_ms, f_m, f_u, f_ob1_u, f_ob1_phi, f_ob2;
    int snapshot_counter = 0; /* XXX Change this when reading from file! */
    
    fftw_complex* correlation_func_x = fftw_alloc_complex(alloc_local);      //Still 2D noise
    fftw_complex* correlation_func_q = fftw_alloc_complex(alloc_local);
    fftw_complex* sqrt_correlation_func_q = fftw_alloc_complex(alloc_local);
    fftw_plan plan_correlation_func;
    
    char filename[sizeof "snapshot1_000.txt"];
    
    srand(time(0)+myrank); /* Set random seed. */
    
    /*initialize the system*/
    initialize(phi, u, zeta, corr_func_data, local_NX, local_0_start);
    Comunicate_boundary(zeta, local_NX);
    Fill_boundary(zeta, local_NX);

    /* Uncomment below when instead of initializing from the start, 
     * reading initial state from a file. */
    /*read_from_file("snapshot1_010.txt", phi[0],local_NX, local_0_start);
    read_from_file("snapshot2_010.txt", phi[1],local_NX, local_0_start);
    read_from_file("snapshot3_010.txt", phi[2],local_NX, local_0_start);
    read_from_file("snapshotu_010.txt", u,local_NX, local_0_start);*/

    sprintf(filename, "snapshot1_%03d.txt", snapshot_counter);
    write_to_file(filename, phi[0], local_NX, local_0_start);
    sprintf(filename, "snapshot2_%03d.txt", snapshot_counter);
    write_to_file(filename, phi[1], local_NX, local_0_start);
    sprintf(filename, "snapshot3_%03d.txt", snapshot_counter);
    write_to_file(filename, phi[2], local_NX, local_0_start);
    sprintf(filename, "snapshotu_%03d.txt", snapshot_counter);
    write_to_file(filename, u, local_NX, local_0_start);

    write_to_file("Initial_zeta.txt", zeta, local_NX, local_0_start);

    /*----------------------------------------------------------------------*/
    /* Fourier transform the correlation function, in this case a Gaussian. */
    plan_correlation_func = fftw_mpi_plan_dft_3d(NX, NY, NZ, correlation_func_x, correlation_func_q, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE);
    
    for(int x=1; x<local_NX+1; x++)
    for(int y=1; y<NY+1; y++)
    for(int z=1; z<NZ+1; z++)
    {
        correlation_func_x[(x-1)*NY*NZ+(y-1)*NZ+(z-1)] = corr_func_data[x][y][z];
    }
    fftw_execute(plan_correlation_func);
    for(int i=0; i<local_NX*NY*NZ; i++)
    {
        sqrt_correlation_func_q[i] = csqrt(correlation_func_q[i]);
    }
    
    free(correlation_func_q);
    free(correlation_func_x);
    fftw_destroy_plan(plan_correlation_func);

    clock_t before = clock();
    /*----------------------------Simulation start----------------------------*/
    for (int i = 1; i < Nt; i++)
    {
        double phi_max=0;
        double dphi_max=0;
        double vol_frac=0;
        double u_max=0;
        double du_max=0;
        double vol_frac_u=0;
        double d[3][3][3];
        /* Communicate the information for neighbor CPUs */
        for(int p=0; p<3; p++)
            Comunicate_boundary(phi[p], local_NX);
        Comunicate_boundary(u, local_NX);
        for(int p=0; p<3; p++)
            Fill_boundary(phi[p], local_NX);
        Fill_boundary(u, local_NX);

        /* Calculate lap of phi first */
        for(int x=1; x<local_NX+1; x++)
        for(int y=1; y<NY+1; y++)
        for(int z=1; z<NZ+1; z++)
        {
            for(int ii=0;ii<3;ii++)
            for(int jj=0;jj<3;jj++)
            for(int kk=0;kk<3;kk++)
                d[ii][jj][kk]=1;
            /* Need to check every set of neighbors and if the difference
             * is more tha0n 90 degrees, flip it 180 degrees. */
            absphisq1 = phi[0][x][y][z]*phi[0][x][y][z] + phi[1][x][y][z]*phi[1][x][y][z] + phi[2][x][y][z]*phi[2][x][y][z];
            absphi = sqrt(absphisq1);

            if ( absphi >= DBL_EPSILON)
            {
                for(int ii=0;ii<3;ii++)
                for(int jj=0;jj<3;jj++)
                for(int kk=0;kk<3;kk++)
                {
                    deltatheta = (phi[0][x][y][z]*phi[0][x+ii-1][y+jj-1][z+kk-1] + phi[1][x][y][z]*phi[1][x+ii-1][y+jj-1][z+kk-1] + phi[2][x][y][z]*phi[2][x+ii-1][y+jj-1][z+kk-1]);
                    if (deltatheta < 0.)
                        d[ii][jj][kk] = -1.;
                    else d[ii][jj][kk] = 1.;
                }
            }
            for(int pp = 0; pp<3; pp++)
            {
                double sum_face=0;
                double sum_edge=0;
                double sum_corner=0;
                for(int ii=0;ii<2;ii++)
                    sum_face += d[2*ii][1][1]*phi[pp][x-1+2*ii][y][z]+d[1][2*ii][1]*phi[pp][x][y-1+2*ii][z]+d[1][1][2*ii]*phi[pp][x][y][z-1+2*ii];
                for(int ii=0;ii<2;ii++)
                for(int jj=0;jj<2;jj++)
                    sum_edge += d[2*ii][2*jj][1]*phi[pp][x-1+2*ii][y-1+2*jj][z]+d[2*ii][1][2*jj]*phi[pp][x-1+2*ii][y][z-1+2*jj]+d[1][2*ii][2*jj]*phi[pp][x][y-1+2*ii][z-1+2*jj];
                for(int ii=0;ii<2;ii++)
                for(int jj=0;jj<2;jj++)
                for(int kk=0;kk<2;kk++)
                    sum_corner += d[2*ii][2*jj][2*kk]*phi[pp][x-1+2*ii][y-1+2*jj][z-1+2*kk];
                
                Laplacian[pp][x][y][z] = 3/(13.*DELTA_X*DELTA_X)*(sum_face+1/2.*sum_edge+1/3.*sum_corner-44/3.*phi[pp][x][y][z]);
            }
        }

        /* Communicate the information for neighbor CPUs */
        for(int p=0; p<3; p++)
        {
            Comunicate_boundary(Laplacian[p], local_NX);
            Fill_boundary(Laplacian[p], local_NX);
        }
        
        /* Save a snapshot every once in a while. */
        if (i % when_to_save == 0)
        {
            snapshot_counter += 1;
            sprintf(filename, "snapshot1_%03d.txt", snapshot_counter);
            write_to_file(filename, phi[0], local_NX, local_0_start);
            sprintf(filename, "snapshot2_%03d.txt", snapshot_counter);
            write_to_file(filename, phi[1], local_NX, local_0_start);
            sprintf(filename, "snapshot3_%03d.txt", snapshot_counter);
            write_to_file(filename, phi[2], local_NX, local_0_start);
            sprintf(filename, "snapshotu_%03d.txt", snapshot_counter);
            write_to_file(filename, u, local_NX, local_0_start);
        }
        
        generate_noise_matrices(eta, sqrt_correlation_func_q, local_NX);

        //For conserved phi field
        double fphi_sum = 0;
        int V_0 = 0;
        
        /*-------------------------Calculation for every grid point-------------------------------*/
        for(int x=1; x<local_NX+1; x++)
        for(int y=1; y<NY+1; y++)
        for(int z=1; z<NZ+1; z++)
        {
            /*--------Calculate orientation angles-------*/
            absphisq1 = phi[0][x][y][z]*phi[0][x][y][z] + phi[1][x][y][z]*phi[1][x][y][z] + phi[2][x][y][z]*phi[2][x][y][z];
            absphi = sqrt(absphisq1);
            absphixy = sqrt(phi[0][x][y][z]*phi[0][x][y][z] + phi[1][x][y][z]*phi[1][x][y][z]);
            
            if (absphisq1 >= DBL_EPSILON)
            {
                phi_s = acos(phi[2][x][y][z]/absphi);
            }
            else
            {
                phi_s = 0.;
            }
            if (absphixy >= DBL_EPSILON)
            {
                theta = atan2(phi[1][x][y][z], phi[0][x][y][z]);
                absphixyinv = 1./absphixy;
            }
            else
            {
                theta = 0.;
                absphixyinv = 0.;
            }
            
            cosphi = cos(phi_s);
            sinphi = sin(phi_s);
            cosphisq = cosphi*cosphi;
            sinphisq = sinphi*sinphi;

            /*----------------Calculate bulk free energy terms-----------------*/
            hex_factor = 1*(exp(-(theta-theta_0)*(theta-theta_0)/0.005) + exp(-(theta-M_PI/3-theta_0)*(theta-M_PI/3-theta_0)/0.005) + exp(-(theta-2*M_PI/3-theta_0)*(theta-2*M_PI/3-theta_0)/0.005) + exp(-(theta-M_PI-theta_0)*(theta-M_PI-theta_0)/0.005) + exp(-(theta+M_PI/3-theta_0)*(theta+M_PI/3-theta_0)/0.005) + exp(-(theta+2*M_PI/3-theta_0)*(theta+2*M_PI/3-theta_0)/0.005) + exp(-(theta+M_PI-theta_0)*(theta+M_PI-theta_0)/0.005));
            hex_factor_p = -2/0.005*((theta-theta_0)*exp(-(theta-theta_0)*(theta-theta_0)/0.005)+(theta-M_PI/3-theta_0)*exp(-(theta-M_PI/3-theta_0)*(theta-M_PI/3-theta_0)/0.005)+(theta-2*M_PI/3-theta_0)*exp(-(theta-2*M_PI/3-theta_0)*(theta-2*M_PI/3-theta_0)/0.005)+(theta-M_PI-theta_0)*exp(-(theta-M_PI-theta_0)*(theta-M_PI-theta_0)/0.005)+(theta+M_PI/3-theta_0)*exp(-(theta+M_PI/3-theta_0)*(theta+M_PI/3-theta_0)/0.005)+(theta+2*M_PI/3-theta_0)*exp(-(theta+2*M_PI/3-theta_0)*(theta+2*M_PI/3-theta_0)/0.005) + (theta+M_PI-theta_0)*exp(-(theta+M_PI-theta_0)*(theta+M_PI-theta_0)/0.005));

            c1 = 0;//-sub_int*(pow(rm/(z+rm-1),12)-2*pow(rm/(z+rm-1),6));
            b1factor = (1. + b1*((1-1/sub_int*c1)*cos(n*theta)+c1*hex_factor)*sinphisq)*(1./(1. + b1));
            b2factor = (1. + b2*cos(n_phi*phi_s))*(1./(1. + b2));

            m_value = (alpha_1/M_PI)*atan(-gamma*(-0.5-0.5*u[x][y][z]));
            factors = -b1factor*b2factor*absphisq1 - \
            (m_value - 1.5 - m_z1*cosphisq)*absphi + (m_value-0.5-m_z1*cosphisq);
            morefactors = (1. - 2./3.*absphi)*m_z1*cosphi*sinphi;

            last_terms_factor = -0.25*absphisq1;
            first_term = -n_phi*b2/(1. + b2)*b1factor*sin(n_phi*phi_s);
            second_term_factor = b1/(1. + b1)*b2factor;
            second_term1 = 2.*((1-1/sub_int*c1)*cos(n*theta)+c1*hex_factor)*sinphi*cosphi;
            second_term2 = ((1-1/sub_int*c1)*n*sin(n*theta)-c1*hex_factor_p);

            dphix = phi[0][x][y][z]*phi[2][x][y][z]*absphixyinv;
            dphiy = phi[1][x][y][z]*phi[2][x][y][z]*absphixyinv;
            dphiz = -absphixy;

            /* Calculate obstacle energy */
            f_ob1_u = (1-u[x][y][z]*u[x][y][z])*(k1_zeta*pow(zeta[x][y][z],2)+k1_phi*absphisq1);
            f_ob1_phi = 2*k1_phi*(u[x][y][z]-pow(u[x][y][z],3)/3.+2./3.);
            if (absphisq1 >= 0.02)
                f_ob2 = k2*pow(zeta[x][y][z],2)*(1./absphi - absphi) + 4./cosh((absphisq1-1.1)/0.05)/cosh((absphisq1-1.1)/0.05);
            else
                f_ob2 = 4./cosh((absphisq1-1.1)/0.05)/cosh((absphisq1-1.1)/0.05);

            /*------------------Compute derivatives with finite difference---------------*/
            absphipx = sqrt(phi[0][x+1][y][z]*phi[0][x+1][y][z] + phi[1][x+1][y][z]*phi[1][x+1][y][z] + phi[2][x+1][y][z]*phi[2][x+1][y][z]);
            absphimx = sqrt(phi[0][x-1][y][z]*phi[0][x-1][y][z] + phi[1][x-1][y][z]*phi[1][x-1][y][z] + phi[2][x-1][y][z]*phi[2][x-1][y][z]);
            dphi_r_x = (absphipx - absphimx)*(1./(2.*DELTA_X));
            absphipy = sqrt(phi[0][x][y+1][z]*phi[0][x][y+1][z] + phi[1][x][y+1][z]*phi[1][x][y+1][z] + phi[2][x][y+1][z]*phi[2][x][y+1][z]);
            absphimy = sqrt(phi[0][x][y-1][z]*phi[0][x][y-1][z] + phi[1][x][y-1][z]*phi[1][x][y-1][z] + phi[2][x][y-1][z]*phi[2][x][y-1][z]);
            dphi_r_y = (absphipy - absphimy)*(1./(2.*DELTA_Y));
            absphipz = sqrt(phi[0][x][y][z+1]*phi[0][x][y][z+1] + phi[1][x][y][z+1]*phi[1][x][y][z+1] + phi[2][x][y][z+1]*phi[2][x][y][z+1]);
            absphimz = sqrt(phi[0][x][y][z-1]*phi[0][x][y][z-1] + phi[1][x][y][z-1]*phi[1][x][y][z-1] + phi[2][x][y][z-1]*phi[2][x][y][z-1]);
            dphi_r_z = (absphipz - absphimz)*(1./(2.*DELTA_Z));
            dphi_r_norm = sqrt(dphi_r_x*dphi_r_x+dphi_r_y*dphi_r_y+dphi_r_z*dphi_r_z);

            dzeta_x = (zeta[x+1][y][z] - zeta[x-1][y][z])*(1./(2.*DELTA_X));
            dzeta_y = (zeta[x][y+1][z] - zeta[x][y-1][z])*(1./(2.*DELTA_X));
            dzeta_z = (zeta[x][y][z+1] - zeta[x][y][z-1])*(1./(2.*DELTA_X));
            grad_zeta_sq = dzeta_x*dzeta_x + dzeta_y*dzeta_y + dzeta_z*dzeta_z;

            du[0] = (u[x+1][y][z] - u[x-1][y][z])*(1./(2.*DELTA_X));
            du[1] = (u[x][y+1][z] - u[x][y-1][z])*(1./(2.*DELTA_X));
            du[2] = (u[x][y][z+1] - u[x][y][z-1])*(1./(2.*DELTA_X)); 
            ddu[0][0] = (u[x+1][y][z] - 2*u[x][y][z] + u[x-1][y][z])/DELTA_X/DELTA_X;
            ddu[1][1] = (u[x][y+1][z] - 2*u[x][y][z] + u[x][y-1][z])/DELTA_X/DELTA_X;
            ddu[2][2] = (u[x][y][z+1] - 2*u[x][y][z] + u[x][y][z-1])/DELTA_X/DELTA_X;
            ddu[0][1] = (u[x+1][y+1][z] - u[x+1][y-1][z] - u[x-1][y+1][z] + u[x-1][y-1][z])/4/DELTA_X/DELTA_X;
            ddu[0][2] = (u[x+1][y][z+1] - u[x+1][y][z-1] - u[x-1][y][z+1] + u[x-1][y][z-1])/4/DELTA_X/DELTA_X;
            ddu[1][2] = (u[x][y+1][z+1] - u[x][y-1][z+1] - u[x][y+1][z-1] + u[x][y-1][z-1])/4/DELTA_X/DELTA_X;
            ddu[1][0] = ddu[0][1];
            ddu[2][0] = ddu[0][2];
            ddu[2][1] = ddu[1][2];
            grad_u_sq = du[0]*du[0] + du[1]*du[1] + du[2]*du[2];

            for(int ii=0; ii<3; ii++)
            {
                dphi[ii][0] = (d[2][1][1]*phi[ii][x+1][y][z] - d[0][1][1]*phi[ii][x-1][y][z])/2/DELTA_X;
                dphi[ii][1] = (d[1][2][1]*phi[ii][x][y+1][z] - d[1][0][1]*phi[ii][x][y-1][z])/2/DELTA_X;
                dphi[ii][2] = (d[1][1][2]*phi[ii][x][y][z+1] - d[1][1][0]*phi[ii][x][y][z-1])/2/DELTA_X;
                ddphi[ii][0][0] = (d[2][1][1]*phi[ii][x+1][y][z] - 2*phi[ii][x][y][z] + d[0][1][1]*phi[ii][x-1][y][z])/DELTA_X/DELTA_X;
                ddphi[ii][1][1] = (d[1][2][1]*phi[ii][x][y+1][z] - 2*phi[ii][x][y][z] + d[1][0][1]*phi[ii][x][y-1][z])/DELTA_X/DELTA_X;
                ddphi[ii][2][2] = (d[1][1][2]*phi[ii][x][y][z+1] - 2*phi[ii][x][y][z] + d[1][1][0]*phi[ii][x][y][z-1])/DELTA_X/DELTA_X;
                ddphi[ii][0][1] = (d[2][2][1]*phi[ii][x+1][y+1][z] - d[2][0][1]*phi[ii][x+1][y-1][z] - d[0][2][1]*phi[ii][x-1][y+1][z] + d[0][0][1]*phi[ii][x-1][y-1][z])/4/DELTA_X/DELTA_X;
                ddphi[ii][0][2] = (d[2][1][2]*phi[ii][x+1][y][z+1] - d[2][1][0]*phi[ii][x+1][y][z-1] - d[0][1][2]*phi[ii][x-1][y][z+1] + d[0][1][0]*phi[ii][x-1][y][z-1])/4/DELTA_X/DELTA_X;
                ddphi[ii][1][2] = (d[1][2][2]*phi[ii][x][y+1][z+1] - d[1][0][2]*phi[ii][x][y-1][z+1] - d[1][2][0]*phi[ii][x][y+1][z-1] + d[1][0][0]*phi[ii][x][y-1][z-1])/4/DELTA_X/DELTA_X;
                ddphi[ii][1][0] = ddphi[ii][0][1];
                ddphi[ii][2][0] = ddphi[ii][0][2];
                ddphi[ii][2][1] = ddphi[ii][1][2];
            }

            if (dphi_r_x != 0)
                v = atan2(dphi_r_y, dphi_r_x);
            else
                v = 0;
            
            double d[3][3][3];
            for(int ii=0;ii<3;ii++)
            for(int jj=0;jj<3;jj++)
            for(int kk=0;kk<3;kk++)
                d[ii][jj][kk]=1;
            /* Need to check every set of neighbors and if the difference
             * is more than 90 degrees, flip it 180 degrees. */
            if (absphi >= DBL_EPSILON)
            {
                for(int ii=0;ii<3;ii++)
                for(int jj=0;jj<3;jj++)
                for(int kk=0;kk<3;kk++)
                {
                    deltatheta = (phi[0][x][y][z]*phi[0][x+ii-1][y+jj-1][z+kk-1] + phi[1][x][y][z]*phi[1][x+ii-1][y+jj-1][z+kk-1] + phi[2][x][y][z]*phi[2][x+ii-1][y+jj-1][z+kk-1]);
                    if (deltatheta < 0.)
                        d[ii][jj][kk] = -1.;
                    else d[ii][jj][kk] = 1.;
                }
            }
            
            /* Calculate Laplacian using 27 points stencil*/
            double sum_face_u=0;
            double sum_edge_u=0;
            double sum_corner_u=0;
            for(int ii=0;ii<2;ii++)
                sum_face_u += u[x-1+2*ii][y][z]+u[x][y-1+2*ii][z]+u[x][y][z-1+2*ii];
            for(int ii=0;ii<2;ii++)
            for(int jj=0;jj<2;jj++)
                sum_edge_u += u[x-1+2*ii][y-1+2*jj][z]+u[x-1+2*ii][y][z-1+2*jj]+u[x][y-1+2*ii][z-1+2*jj];
            for(int ii=0;ii<2;ii++)
            for(int jj=0;jj<2;jj++)
            for(int kk=0;kk<2;kk++)
                sum_corner_u += u[x-1+2*ii][y-1+2*jj][z-1+2*kk];
            Laplacian_u[x][y][z] = 3/(13.*DELTA_X*DELTA_X)*(sum_face_u+1/2.*sum_edge_u+1/3.*sum_corner_u-44/3.*u[x][y][z]);

            /* Calculate the lap of lap of phi */
            double LapLaplacian[3];
            for(int pp = 0; pp<3; pp++)
            {
                double sum_face=0;
                double sum_edge=0;
                double sum_corner=0;
                for(int ii=0;ii<2;ii++)
                    sum_face += d[2*ii][1][1]*Laplacian[pp][x-1+2*ii][y][z]+d[1][2*ii][1]*Laplacian[pp][x][y-1+2*ii][z]+d[1][1][2*ii]*Laplacian[pp][x][y][z-1+2*ii];
                for(int ii=0;ii<2;ii++)
                for(int jj=0;jj<2;jj++)
                    sum_edge += d[2*ii][2*jj][1]*Laplacian[pp][x-1+2*ii][y-1+2*jj][z]+d[2*ii][1][2*jj]*Laplacian[pp][x-1+2*ii][y][z-1+2*jj]+d[1][2*ii][2*jj]*Laplacian[pp][x][y-1+2*ii][z-1+2*jj];
                for(int ii=0;ii<2;ii++)
                for(int jj=0;jj<2;jj++)
                for(int kk=0;kk<2;kk++)
                    sum_corner += d[2*ii][2*jj][2*kk]*Laplacian[pp][x-1+2*ii][y-1+2*jj][z-1+2*kk];
                
                LapLaplacian[pp] = 3/(13.*DELTA_X*DELTA_X)*(sum_face+1/2.*sum_edge+1/3.*sum_corner-44/3.*Laplacian[pp][x][y][z]);
            }
            
            l1 =2*(epsc)*Laplacian[0][x][y][z];
            l2 =2*(epsc)*Laplacian[1][x][y][z];
            l3 =2*(epsc)*Laplacian[2][x][y][z];

            /*-----------Calculate all the variations of the anisotropic interfacial energy.-------------*/
            phix = phi[0][x][y][z];
            phiy = phi[1][x][y][z];
            phiz = phi[2][x][y][z];
            if(phix*phix+phiy*phiy+phiz*phiz <= 0.00001 || (pow(phix*dphi[0][2] + phiy*dphi[1][2] + phiz*dphi[2][2],2) + 
       pow(phix*dphi[0][1] + phiy*dphi[1][1] + phiz*dphi[2][1],2) + 
       pow(phix*dphi[0][0] + phiy*dphi[1][0] + phiz*dphi[2][0],2)) <= 0.00001)
            {
                f2x = 0;
                f4x = 0;
                f2y = 0;
                f4y = 0;
                f2z = 0;
                f4z = 0;
            }
            else
            {
                f2x = calc_f2x(phix, phiy, phiz, dphi, ddphi);
                f2y = calc_f2y(phix, phiy, phiz, dphi, ddphi);
                f2z = calc_f2z(phix, phiy, phiz, dphi, ddphi);
                f4x = calc_f4x(phix, phiy, phiz, dphi, ddphi);
                f4y = calc_f4y(phix, phiy, phiz, dphi, ddphi);
                f4z = calc_f4z(phix, phiy, phiz, dphi, ddphi);
            }
            /* Set threshold to avoid numerical divergence */
	        double thre = 1.e4;
            f2x = min(f2x,thre);
            f2x = max(f2x,-thre);
            f2y = min(f2y,thre);
            f2y = max(f2y,-thre);
            f2z = min(f2z,thre);
            f2z = max(f2z,-thre);
            f4x = min(f4x,thre);
            f4x = max(f4x,-thre);
            f4y = min(f4y,thre);
            f4y = max(f4y,-thre);
            f4z = min(f4z,thre);
            f4z = max(f4z,-thre);

            /*---------------Calculate the first variatio of the anisotropic interfacial energy of phi-------------*/
            double f_anix, f_aniy, f_aniz;
            f_anix = (epsc)*(eps_ani_2*f2x-eps_ani_4*f4x);
            f_aniy = (epsc)*(eps_ani_2*f2y-eps_ani_4*f4y);
            f_aniz = (epsc)*(eps_ani_2*f2z-eps_ani_4*f4z);

            /*--------------Calculate the first variation of the interfacial energy of phi|sub and phi|u w.r.t. phi----------*/
            f_cs1 = 2*eps_cs*grad_zeta_sq;
            f_cs2 = 2*eps_cs*gamma_cs*(dzeta_x*phix+dzeta_y*phiy+dzeta_z*phiz);
            f_cm1 = 2*eps_cm*grad_u_sq;
            f_cm2 = 2*eps_cm*gamma_cm*(du[0]*phix+du[1]*phiy+du[2]*phiz); 

            /*-------------Calculate the first variation of the bulk energy of u------------*/
            f_u = -4*ku*u[x][y][z]*(1-u[x][y][z]*u[x][y][z]);
            /*-------------Calculate the first variation of the interfacial energy of u------------*/
            f_m = -2*eps_m*Laplacian_u[x][y][z];
            /*-------------Calculate the first variation of the interfacial energy of u|phi and u|sub w.r.t. u----------*/
            f_mc1 = -2*eps_cm*absphisq1*Laplacian_u[x][y][z];
            f_mc2 = -2*eps_cm*gamma_cm*(pow(phiz,2)*ddu[2][2] + pow(phiy,2)*ddu[1][1] + phiz*(2*dphi[2][2]*du[2] + dphi[1][2]*du[1] + 2*phiy*ddu[1][2] + du[2]*(dphi[1][1] + dphi[0][0]) + dphi[0][2]*du[0] + 2*phix*ddu[0][2]) + phiy*(du[2]*dphi[2][1] + dphi[2][2]*du[1] + 2*dphi[1][1]*du[1] + du[1]*dphi[0][0] + dphi[0][1]*du[0] + 2*phix*ddu[0][1]) + phix*(du[1]*dphi[1][0] + du[2]*dphi[2][0] + dphi[2][2]*du[0] + dphi[1][1]*du[0] + 2*dphi[0][0]*du[0] + phix*ddu[0][0]));
            f_ms = 2*eps_ms*grad_zeta_sq;
            /*-------------Calculate the first variation of the total free energy w.r.t. u--------------*/
            dFdu[x][y][z] = f_u + f_ob1_u + f_m + f_mc1 + f_mc2 + f_ms;

            /*-----------------Calculate the anisotropic mobility for phi---------------*/
            isotropic_mobility = (DELTA_T*mobility1_max)*0.5*(1. + tanh(-(absphi-cbar1)/ctilde1));
            if (absphi >= 0.001 && dphi_r_norm >= 0.001)
            {
                Mphi = isotropic_mobility;//*(1+deltaM*\
            (fabs(1-pow((dphi_r_x*phi[0][x][y][z]+dphi_r_y*phi[1][x][y][z]+dphi_r_z*phi[2][x][y][z])/(dphi_r_norm*absphi),2))))/deltaM*(sinphisq + (1.+cos(nbeta*v))*cosphisq/2.);
            }else
            {
                Mphi = isotropic_mobility;
            }
            
            /* Fluctuations */
            if (absphi >= DBL_EPSILON && dphi_r_norm >=DBL_EPSILON)
                fluct_amp = absphi*(1.-absphi)*sqrt(32.*kBT*isotropic_mobility)*(pow(1-pow((dphi_r_x*phi[0][x][y][z]+dphi_r_y*phi[1][x][y][z]+dphi_r_z*phi[2][x][y][z])/(dphi_r_norm*absphi),2),2)*sinphisq+cosphisq);
            else
                fluct_amp = 0;
            
            /*------------------Calculate the change of phi----------------------*/
            fphi[0][x][y][z] = Mphi*(-f_anix - eps_higher_order*2*LapLaplacian[0] - (f_ob1_phi+f_ob2+f_cs1+f_cm1)*phi[0][x][y][z] - f_cs2*dzeta_x - f_cm2*du[0] + l1 + factors*phi[0][x][y][z] + morefactors*dphix + last_terms_factor*(first_term*dphix + second_term_factor*(second_term1*dphix + second_term2*phi[1][x][y][z]))) + eta[0][x][y][z]*fluct_amp;
            fphi[1][x][y][z] = Mphi*(-f_aniy - eps_higher_order*2*LapLaplacian[1] - (f_ob1_phi+f_ob2+f_cs1+f_cm1)*phi[1][x][y][z] - f_cs2*dzeta_y - f_cm2*du[1] + l2 + factors*phi[1][x][y][z] + morefactors*dphiy + last_terms_factor*(first_term*dphiy + second_term_factor*(second_term1*dphiy - second_term2*phi[0][x][y][z]))) + eta[1][x][y][z]*fluct_amp;
            fphi[2][x][y][z] = Mphi*(-f_aniz - eps_higher_order*2*LapLaplacian[2] - (f_ob1_phi+f_ob2+f_cs1+f_cm1)*phi[2][x][y][z] - f_cs2*dzeta_z - f_cm2*du[2] + l3 + factors*phi[2][x][y][z] + morefactors*dphiz + last_terms_factor*(first_term*dphiz + second_term_factor*second_term1*dphiz)) + eta[2][x][y][z]*fluct_amp;

            // For conserved phi field
            /*if(absphi >= 0.02)
            {
                fphi_sum += (fphi[0][x][y][z]*phi[0][x][y][z]+fphi[1][x][y][z]*phi[1][x][y][z]+fphi[2][x][y][z]*phi[2][x][y][z])/absphi;
                V_0 += 1;
            }*/

        }

        Comunicate_boundary(dFdu, local_NX);
        Fill_boundary(dFdu, local_NX);

        /*-------------Calculate the lap of dFdu--------------*/
        for(int x=1; x<local_NX+1; x++)
        for(int y=1; y<NY+1; y++)
        for(int z=1; z<NZ+1; z++)
        {
            double sum_face_dFdu=0;
            double sum_edge_dFdu=0;
            double sum_corner_dFdu=0;
            for(int ii=0;ii<2;ii++)
                sum_face_dFdu += dFdu[x-1+2*ii][y][z]+dFdu[x][y-1+2*ii][z]+dFdu[x][y][z-1+2*ii];
            for(int ii=0;ii<2;ii++)
            for(int jj=0;jj<2;jj++)
                sum_edge_dFdu += dFdu[x-1+2*ii][y-1+2*jj][z]+dFdu[x-1+2*ii][y][z-1+2*jj]+dFdu[x][y-1+2*ii][z-1+2*jj];
            for(int ii=0;ii<2;ii++)
            for(int jj=0;jj<2;jj++)
            for(int kk=0;kk<2;kk++)
                sum_corner_dFdu += dFdu[x-1+2*ii][y-1+2*jj][z-1+2*kk];
            Laplacian_dFdu[x][y][z] = 3/(13.*DELTA_X*DELTA_X)*(sum_face_dFdu+1/2.*sum_edge_dFdu+1/3.*sum_corner_dFdu-44/3.*dFdu[x][y][z]);
        }

        /*-----------Update all the variables-------------*/
        // For conserved phi field
        /*MPI_Allreduce(MPI_IN_PLACE, &fphi_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &V_0, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);*/

        for(int x=1; x<local_NX+1; x++)
        for(int y=1; y<NY+1; y++)
        for(int z=1; z<NZ+1; z++)
        {
            absphisq1 = phi[0][x][y][z]*phi[0][x][y][z] + phi[1][x][y][z]*phi[1][x][y][z] + phi[2][x][y][z]*phi[2][x][y][z];
            absphi = sqrt(absphisq1);
            // For conserved phi field
            /*for(int p=0; p<3; p++)
            {
                if(absphi >= DBL_EPSILON)
                    fphi[p][x][y][z] += -fphi_sum*phi[p][x][y][z]/absphi/V_0;
            }*/

            newphi1 = phi[0][x][y][z] + fphi[0][x][y][z];
            newphi2 = phi[1][x][y][z] + fphi[1][x][y][z];
            newphi3 = phi[2][x][y][z] + fphi[2][x][y][z];
            dabsphi = sqrt(newphi1*newphi1 + newphi2*newphi2 + newphi3*newphi3) - absphi;
            
            DELTA = 0.8;
            /*------------------Calculate the change of u------------------------*/
            fu[x][y][z] = DELTA_T*D_U*Laplacian_dFdu[x][y][z] - 2*dabsphi*(1./DELTA);

            /* Update phi and u for the next time step */
            phi[0][x][y][z] += fphi[0][x][y][z];
            phi[1][x][y][z] += fphi[1][x][y][z];
            phi[2][x][y][z] += fphi[2][x][y][z];
            u[x][y][z] += fu[x][y][z];

            if(pow(phi[0][x][y][z],2)+pow(phi[1][x][y][z],2)+pow(phi[2][x][y][z],2)>0.5)
                vol_frac += 1.;
            for(int ii=0; ii<3; ii++)
            {
                phi_max = max(phi_max, fabs(phi[ii][x][y][z]));
                dphi_max = max(dphi_max, fabs(fphi[ii][x][y][z]));
            }
            if(u[x][y][z]>0.5)
                vol_frac_u += 1.;
            u_max = max(u_max, fabs(u[x][y][z]));
            du_max = max(du_max, fabs(fu[x][y][z]));
        }
        MPI_Allreduce(MPI_IN_PLACE, &phi_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &dphi_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &vol_frac, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &u_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &du_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &vol_frac_u, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        /* Print information at every timestep */
        if(myrank ==0 /*&& i%10==0*/)
        {
            printf("step= %d, phi_max= %f, dphi_max= %f, vol_frac= %f\n", i, phi_max, dphi_max, vol_frac/NX/NY/NZ);
            printf("step= %d, u_max= %f, du_max= %f, vol_frac_u= %f\n", i, u_max, du_max, vol_frac_u/NX/NY/NZ);
        }
    }
    clock_t difference = clock()-before;
    int time = difference/CLOCKS_PER_SEC;
    if(myrank ==0)
        printf("Time elapse = %d\n",time);
    
    /* Final snapshot. */
    snapshot_counter += 1;
    sprintf(filename, "snapshot1_%03d.txt", snapshot_counter);
    write_to_file(filename, phi[0], local_NX, local_0_start);
    sprintf(filename, "snapshot2_%03d.txt", snapshot_counter);
    write_to_file(filename, phi[1], local_NX, local_0_start);
    sprintf(filename, "snapshot3_%03d.txt", snapshot_counter);
    write_to_file(filename, phi[2], local_NX, local_0_start);
    sprintf(filename, "snapshotu_%03d.txt", snapshot_counter);
    write_to_file(filename, u, local_NX, local_0_start);
    
    /* Free allocated memory */
    for(int x=0; x<local_NX+2; x++)
    {
        for(int y=0; y<NY+2; y++)
        {
            free(phi[0][x][y]);
            free(phi[1][x][y]);
            free(phi[2][x][y]);
            free(fphi[0][x][y]);
            free(fphi[1][x][y]);
            free(fphi[2][x][y]);
            free(eta[0][x][y]);
            free(eta[1][x][y]);
            free(eta[2][x][y]);
            free(u[x][y]);
            free(fu[x][y]);
            free(zeta[x][y]);
            free(Laplacian[0][x][y]);
            free(Laplacian[1][x][y]);
            free(Laplacian[2][x][y]);
            free(Laplacian_u[x][y]);
            free(Laplacian_dFdu[x][y]);
            free(dFdu[x][y]);
        }
        free(phi[0][x]);
        free(phi[1][x]);
        free(phi[2][x]);
        free(fphi[0][x]);
        free(fphi[1][x]);
        free(fphi[2][x]);
        free(eta[0][x]);
        free(eta[1][x]);
        free(eta[2][x]);
        free(u[x]);
        free(fu[x]);
        free(zeta[x]);
        free(Laplacian[0][x]);
        free(Laplacian[1][x]);
        free(Laplacian[2][x]);
        free(Laplacian_u[x]);
        free(Laplacian_dFdu[x]);
        free(dFdu[x]);
    }
    for(int x=0; x<local_NX; x++)
        free(corr_func_data[x]);
    free(phi[0]);
    free(phi[1]);
    free(phi[2]);
    free(fphi[0]);
    free(fphi[1]);
    free(fphi[2]);
    free(eta[0]);
    free(eta[1]);
    free(eta[2]);
    free(u);
    free(fu);
    free(corr_func_data);
    free(sqrt_correlation_func_q);
    free(phi);
    free(fphi);
    free(eta);
    free(zeta);
    free(Laplacian[0]);
    free(Laplacian[1]);
    free(Laplacian[2]);
    free(Laplacian_u);
    free(Laplacian_dFdu);
    free(dFdu);
    fftw_cleanup();
    
    MPI_Finalize();
    return 0;
}

void generate_noise_matrices(double **** eta, fftw_complex* sqrt_correlation_func_q, int local_NX)
{
    fftw_complex* in_x = (fftw_complex*) malloc(local_NX*NY*NZ*sizeof(fftw_complex));
    fftw_complex* in_q = (fftw_complex*) malloc(local_NX*NY*NZ*sizeof(fftw_complex));
    fftw_complex* correlated_q = (fftw_complex*) malloc(local_NX*NY*NZ*sizeof(fftw_complex));
    fftw_complex* correlated_x = (fftw_complex*) malloc(local_NX*NY*NZ*sizeof(fftw_complex));
    fftw_plan plan_f;
    fftw_plan plan_b;
    int j, k, l;
    
    /* Some of the below code may be particularly inelegant and involve copy-pasting
     * instead of better coding practies. */
    
    plan_f = fftw_mpi_plan_dft_3d(NX, NY, NZ, in_x, in_q, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE);
    
    /* Generate uncorrelated noise. */
    for (j = 0; j < local_NX; j++)
    for (k = 0; k < NY; k++)
    for (l = 0; l < NZ; l++)
        in_x[j*NY*NZ + k*NZ + l] = gaussrand();
    
    /*-------------------------------------------*/
    /* Fourier transform the uncorrelated noise. */
    fftw_execute(plan_f);
    
    /* Now multiply by sqrt(correlation_func_q). */
    for (j = 0; j < local_NX; j++)
    for (k = 0; k < NY; k++)
    for (l = 0; l < NZ; l++)
        correlated_q[j*NY*NZ + k*NZ + l] = sqrt_correlation_func_q[j*NY*NZ + k*NZ + l]*in_q[j*NY*NZ + k*NZ + l];
    
    /* Now transform it back. */
    plan_b = fftw_mpi_plan_dft_3d(NX, NY, NZ, correlated_q, correlated_x, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan_b);
    
    /* Rescale by factor of NX*NY, convert to a 2d matrix,
     and scale appropriately. */
    for (j = 0; j < local_NX; j++)
    for (k = 0; k < NY; k++)
    for (l = 0; l < NZ; l++)
    {
        eta[0][j+1][k+1][l+1] = creal(correlated_x[j*NY*NZ + k*NZ + l])*(1./((double)NX*NY*NZ));
        /* Generate uncorrelated noise. */
        in_x[j*NY*NZ + k*NZ + l] = gaussrand();
    }
    
    /*-------------------------------------------*/
    /* Fourier transform the uncorrelated noise. */
    fftw_execute(plan_f);
    
    /* Now multiply by sqrt(correlation_func_q). */
    for (j = 0; j < local_NX; j++)
    for (k = 0; k < NY; k++)
    for (l = 0; l < NZ; l++)
        correlated_q[j*NY*NZ + k*NZ + l] = sqrt_correlation_func_q[j*NY*NZ + k*NZ + l]*in_q[j*NY*NZ + k*NZ + l];
    
    /* Now transform it back. */
    fftw_execute(plan_b);
    
    /* Rescale by factor of NX*NY, convert to a 2d matrix,
     and scale appropriately. */
    for (j = 0; j < local_NX; j++)
    for (k = 0; k < NY; k++)
    for (l = 0; l < NZ; l++)
    {
        eta[1][j+1][k+1][l+1] = creal(correlated_x[j*NY*NZ + k*NZ + l])*(1./((double)NX*NY*NZ));
        /* Generate uncorrelated noise. */
        in_x[j*NY*NZ + k*NZ + l] = gaussrand();
    }
    
    /*-------------------------------------------*/
    /* Fourier transform the uncorrelated noise. */
    fftw_execute(plan_f);
    
    /* Now multiply by sqrt(correlation_func_q). */
    for (j = 0; j < local_NX; j++)
    for (k = 0; k < NY; k++)
    for (l = 0; l < NZ; l++)
        correlated_q[j*NY*NZ + k*NZ + l] = sqrt_correlation_func_q[j*NY*NZ + k*NZ + l]*in_q[j*NY*NZ + k*NZ + l];
    
    /* Now transform it back. */
    fftw_execute(plan_b);
    
    /* Rescale by factor of NX*NY, convert to a 2d matrix,
     and scale appropriately. */
    for (j = 0; j < local_NX; j++)
    for (k = 0; k < NY; k++)
    for (l = 0; l < NZ; l++)
    {
        eta[1][j+1][k+1][l+1] = creal(correlated_x[j*NY*NZ + k*NZ + l])*(1./((double)NX*NY*NZ));
    }
    
    free(in_x);
    free(in_q);
    free(correlated_q);
    free(correlated_x);
    
    fftw_destroy_plan(plan_f);
    fftw_destroy_plan(plan_b);
    
    return;
}

/* Write results to text file output. */
void write_to_file(char *filename, double ***data, int local_NX, int local_0_start)    
{
    int tag = 0;
    int np, rank;
    MPI_Status status;
    double * resize = (double*) malloc((local_NX*NY*NZ)*sizeof(double));
    
    for(int z=1; z<NZ+1; z++)
    for(int x=1; x<local_NX+1; x++)
    for(int y=1; y<NY+1; y++)
        resize[(x-1)*NY*NZ+(y-1)*NZ+(z-1)] = data[x][y][z];
    
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if ( rank == 0 )
    {
        double * buffer = (double*) malloc((NX*NY*NZ)*sizeof(double));
        memcpy(buffer, resize, local_NX*NY*NZ*sizeof(double));
        
        for (int i=1; i<np; i++)
        {
            int count, offset;
            MPI_Recv(&offset, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
            MPI_Recv(&count, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
            MPI_Recv(buffer + offset, count, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);
        }
        FILE *fp;
        fp = fopen(filename, "w");
        for(int z=0; z<NZ; z++)
        {
            for(int y=0; y<NY; y++)
            {
                for(int x=0; x<NX; x++)
                    fprintf(fp, "%0.3f ", buffer[x*NY*NZ+y*NZ+z]);
                fprintf(fp, "\n");
            }
            fprintf(fp, "\n");
        }
        free(buffer);
	fclose(fp);
    }
    else
    {
        int offset = local_0_start*NY*NZ;
        int count = local_NX*NY*NZ;
        MPI_Send(&offset, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
        MPI_Send(&count, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
        MPI_Send(resize, count, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
    }
    free(resize);
}

/* Read from a text file a matrix. */
void read_from_file(char *filename, double ***data, int local_NX, int local_0_start)
{
    int np, rank;
    MPI_Status status;
    int tag = 0;
    
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double * buffer_data = (double*) malloc((local_NX*NY*NZ)*sizeof(double));

    if(rank==0)
    {
        double * temp_data = (double*) malloc((NX*NY*NZ)*sizeof(double));
        FILE *fp;
        int ndx = 0;
        float number = 0.;

        // Read from file
        fp = fopen(filename, "r");
        if (fp != NULL)
        {
            for(int z=0; z<NZ; z++)
            for(int y=0; y<NY; y++)
            for(int x=0; x<NX; x++)
            {
                ndx = x*NY*NZ+y*NZ+z;
                fscanf(fp, "%f ", &number);
                temp_data[ndx] = number;
            }
        }
        else
        {
            printf("error! file does not exist. \n");
        }
        fclose(fp);

        // Transfer to self
        memcpy(buffer_data, temp_data, local_NX*NY*NZ*sizeof(double));

        // Transfer to slaves
        for (int i=1; i<np; i++)
        {
            int offset, count;
            // master asks for local information from slave
            MPI_Recv(&offset, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
            MPI_Recv(&count, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);

            // master sends data
            MPI_Send(temp_data+offset, count, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
        }

        free(temp_data);
    }
    else
    {
        // slaves send local information to master
        int offset = local_0_start*NY*NZ;
        int count = local_NX*NY*NZ;
        MPI_Send(&offset, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
        MPI_Send(&count, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);

        // slave recv data
        MPI_Recv(buffer_data, count, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
    }
    
    for(int x=1; x<local_NX+1; x++)
    for(int y=1; y<NY+1; y++)
    for(int z=1; z<NZ+1; z++)
    {
        int ndx = (x-1)*NY*NZ+(y-1)*NZ+(z-1);
        data[x][y][z] = buffer_data[ndx];
    }

    free(buffer_data);
}

/* Modified version of: http://c-faq.com/lib/gaussian.html */
double gaussrand(void)
{
	static double V1, V2, S;
	static int phase = 0;
	double X;
    
	if(phase == 0) {
		do {
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;
            
			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);
        
		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);
    
	phase = 1 - phase;
    /*X = X*sd + mean;*/ /* Assume mean=0, sd=1. */
    
	return X;
}

/* Create a 3D array. */
double*** Make3DDoubleArray(int arraySizeX, int arraySizeY, int arraySizeZ)
{
    double*** theArray;
    int i;
    
    theArray = (double***) malloc((arraySizeX+2)*sizeof(double**));
    
    for (i = 0; i < arraySizeX+2; i++)
    {
        theArray[i] = (double**) malloc((arraySizeY+2)*sizeof(double*));
        for (int j=0; j<arraySizeY+2; j++)
            theArray[i][j] = (double*) malloc((arraySizeZ+2)*sizeof(double));
    }
    return theArray;
}

void Comunicate_boundary(double *** data, int local_NX)
{
    int np, rank;
    MPI_Status status;
    
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    double * datam_s = (double*) malloc((NY*NZ)*sizeof(double));
    double * datap_s = (double*) malloc((NY*NZ)*sizeof(double));
    double * datam_r = (double*) malloc((NY*NZ)*sizeof(double));
    double * datap_r = (double*) malloc((NY*NZ)*sizeof(double));
    
    // Prepare boundary data
    for(int y=1; y<NY+1; y++)
    for(int z=1; z<NZ+1; z++)
    {
        datam_s[(y-1)*NZ+(z-1)] = data[1][y][z];
        datap_s[(y-1)*NZ+(z-1)] = data[local_NX][y][z];
    }
    
    // Prepare dest proc
    int rankm = rank - 1;
    int rankp = rank + 1;
    if(rankm<0) rankm = np-1;
    if(rankp>=np) rankp = 0;
   
    // Send and receive boundary data
    MPI_Sendrecv(datam_s, NY*NZ, MPI_DOUBLE, rankm, 1, datap_r, NY*NZ, MPI_DOUBLE, rankp, 1, MPI_COMM_WORLD, &status);
    MPI_Sendrecv(datap_s, NY*NZ, MPI_DOUBLE, rankp, 1, datam_r, NY*NZ, MPI_DOUBLE, rankm, 1, MPI_COMM_WORLD, &status); 
    
    for(int y=1; y<NY+1; y++)
    for(int z=1; z<NZ+1; z++)
    {
        data[0][y][z] = datam_r[(y-1)*NZ+(z-1)];
        data[local_NX+1][y][z] = datap_r[(y-1)*NZ+(z-1)];
    }

    //Fill the boundary in x direction as no-flow boundary
    /*if(rank==0)
    {
        for(int y=1; y<NY+1; y++)
        for(int z=1; z<NZ+1; z++)
            data[0][y][z] = data[1][y][z];
    }
    if(rank==np-1)
    {
	for(int y=1; y<NY+1; y++)
        for(int z=1; z<NZ+1; z++)
            data[local_NX+1][y][z] = data[local_NX][y][z];
    }*/
    free(datam_s);
    free(datap_s);
    free(datam_r);
    free(datap_r);
}

void Fill_boundary(double *** data, int local_NX)
{

    for(int x=0; x<local_NX+2; x++)
    for(int y=0; y<NY+2; y++)
    for(int z=0; z<NZ+2; z++)
    {
        int x1=x;
        int y1=y;
        int z1=z;
        if(y==0)              // Fill the boundary in y direction as periodic boundary
            y1=NY;
        else if(y==NY+1)
            y1=1;
        if(z==0)              // Fill the boundary in z direction as no-flow boundary
            z1=1;
        else if(z==NZ+1)
            z1=NZ;
        data[x][y][z] = data[x1][y1][z1];
    }
}

void initialize(double **** phi, double *** u, double *** zeta, double *** corr_func_data, int local_NX, int local_0_start)
{
    int tag = 0;
    double * buffer_phi0 = (double*) malloc((local_NX*NY*NZ)*sizeof(double));
    double * buffer_phi1 = (double*) malloc((local_NX*NY*NZ)*sizeof(double));
    double * buffer_phi2 = (double*) malloc((local_NX*NY*NZ)*sizeof(double));
    double * buffer_u = (double*) malloc((local_NX*NY*NZ)*sizeof(double));
    double * buffer_zeta = (double*) malloc((local_NX*NY*NZ)*sizeof(double));
    double * buffer_corr_func_data = (double*) malloc((local_NX*NY*NZ)*sizeof(double));
    int np, rank;
    MPI_Status status;
    
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if ( rank == 0 )
    {
        // master prepares and sends data
        double * temp_phi0 = (double*) malloc((NX*NY*NZ)*sizeof(double));
        double * temp_phi1 = (double*) malloc((NX*NY*NZ)*sizeof(double));
        double * temp_phi2 = (double*) malloc((NX*NY*NZ)*sizeof(double));
        double * temp_u = (double*) malloc((NX*NY*NZ)*sizeof(double));
        double * temp_zeta = (double*) malloc((NX*NY*NZ)*sizeof(double));
        double * temp_corr_func_data = (double*) malloc((NX*NY*NZ)*sizeof(double));
        
        /* Generate initial state. */
        /* First fill the arrays with zeros. */
        for(int x=0; x<NX; x++)
            for(int y=0; y<NY; y++)
                for(int z=0; z<NZ; z++)
                {
                    int ndx = x*NY*NZ+y*NZ+z;
                    temp_phi0[ndx] = INIT_PHI_AMP;
                    temp_phi1[ndx] = INIT_PHI_AMP;
                    temp_phi2[ndx] = INIT_PHI_AMP;
                    temp_corr_func_data[ndx] = (1./(2.*M_PI*corr_length_sq))*exp(-0.5/corr_length_sq*((double)((x-NX/2)*(x-NX/2) + (y-NY/2)*(y-NY/2) + (z)*(z)))*(1./((double)NX*NY)));
                }

        initialize_u(temp_u);
        initialize_zeta(temp_zeta);
        
        /* Add the initial nucleus in the center. */
        /* the angle phi goes from 0 to pi, as measured from the z axis. */
        for(int x=0; x<NX; x++)
        for(int y=0; y<NY; y++)
        for(int z=0; z<NZ; z++)
        {
            int ndx = x*NY*NZ+y*NZ+z;
            //double Radius = max(fabs(x-NX/2),fabs(y-NY/2));
            //Radius = max(Radius,fabs(z-40));
            double Radius = sqrt((x-NX/2)*(x-NX/2)+(y-NY/2)*(y-NY/2)+(z-NZ/2)*(z-NZ/2));//sqrt((x-NX/2)*(x-NX/2)+(y-NY/2)*(y-NY/2)+(z-NY/2)*(z-NY/2));
            double diffuse = -0.5*(tanh((Radius-5)))+0.5;
            temp_phi0[ndx] = diffuse*cos(INIT_THETA)*sin(INIT_PHI);
            temp_phi1[ndx] = diffuse*sin(INIT_THETA)*sin(INIT_PHI);
            temp_phi2[ndx] = diffuse*cos(INIT_PHI);
        }

        // transfer to self
        memcpy(buffer_phi0, temp_phi0, local_NX*NY*NZ*sizeof(double));
        memcpy(buffer_phi1, temp_phi1, local_NX*NY*NZ*sizeof(double));
        memcpy(buffer_phi2, temp_phi2, local_NX*NY*NZ*sizeof(double));
        memcpy(buffer_u, temp_u, local_NX*NY*NZ*sizeof(double));
        memcpy(buffer_zeta, temp_zeta, local_NX*NY*NZ*sizeof(double));
        memcpy(buffer_corr_func_data, temp_corr_func_data, local_NX*NY*NZ*sizeof(double));
        // transfer to slaves
        for (int i=1; i<np; i++)
        {
            int offset, count;
            // master asks for local information from slave
            MPI_Recv(&offset, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
            MPI_Recv(&count, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
            printf("%d %d\n", i, offset);
            
            // master sends data
            MPI_Send(temp_phi0+offset, count, MPI_DOUBLE, i, 3, MPI_COMM_WORLD);
            MPI_Send(temp_phi1+offset, count, MPI_DOUBLE, i, 4, MPI_COMM_WORLD);
            MPI_Send(temp_phi2+offset, count, MPI_DOUBLE, i, 5, MPI_COMM_WORLD);
            MPI_Send(temp_u+offset, count, MPI_DOUBLE, i, 6, MPI_COMM_WORLD);
            MPI_Send(temp_zeta+offset, count, MPI_DOUBLE, i, 8, MPI_COMM_WORLD);
            MPI_Send(temp_corr_func_data+offset, count, MPI_DOUBLE, i, 9, MPI_COMM_WORLD);
        }
        
        free(temp_phi0);
        free(temp_phi1);
        free(temp_phi2);
        free(temp_u);
        free(temp_zeta);
        free(temp_corr_func_data);
    }
    else
    {
        // slaves send local information to master
        int offset = local_0_start*NY*NZ;
        int count = local_NX*NY*NZ;
        MPI_Send(&offset, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
        MPI_Send(&count, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
        
        // slave recv data
        MPI_Recv(buffer_phi0, count, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &status);
        MPI_Recv(buffer_phi1, count, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, &status);
        MPI_Recv(buffer_phi2, count, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD, &status);
        MPI_Recv(buffer_u, count, MPI_DOUBLE, 0, 6, MPI_COMM_WORLD, &status);
        MPI_Recv(buffer_zeta, count, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD, &status);
        MPI_Recv(buffer_corr_func_data, count, MPI_DOUBLE, 0, 9, MPI_COMM_WORLD, &status);
    }
    
    for(int x=1; x<local_NX+1; x++)
    for(int y=1; y<NY+1; y++)
    for(int z=1; z<NZ+1; z++)
    {
        int ndx = (x-1)*NY*NZ+(y-1)*NZ+(z-1);
        phi[0][x][y][z] = buffer_phi0[ndx];
        phi[1][x][y][z] = buffer_phi1[ndx];
        phi[2][x][y][z] = buffer_phi2[ndx];
        u[x][y][z] = buffer_u[ndx];
        zeta[x][y][z] = buffer_zeta[ndx];
        corr_func_data[x][y][z] = buffer_corr_func_data[ndx];
    }
    
    free(buffer_phi0);
    free(buffer_phi1);
    free(buffer_phi2);
    free(buffer_u);
    free(buffer_zeta);
    free(buffer_corr_func_data);
        
}

void initialize_u(double * temp_u)
{
    double Radius;
    for(int x=0; x<NX; x++)
    for(int y=0; y<NY; y++)
    for(int z=0; z<NZ; z++)
    {
        int ndx = x*NY*NZ+y*NZ+z;
        temp_u[ndx] = 1.;
        /*Radius = max(fabs(x-NX/2),fabs(y-NY/2));
        Radius = max(Radius,fabs(z-30));*/
        //Radius = sqrt((x-NX/2)*(x-NX/2)+(y-NY/2)*(y-NY/2)+(z-70)*(z-70));
        //temp_u[ndx] += -(tanh((-z+10.)/2.));
        //temp_u[ndx] += -(tanh((Radius-20)/2.));
    }
}

void initialize_zeta(double * temp_zeta)
{
    double Radius;
    for(int x=0; x<NX; x++)
    for(int y=0; y<NY; y++)
    for(int z=0; z<NZ; z++)
    {
        int ndx = x*NY*NZ+y*NZ+z;
        temp_zeta[ndx] = 0;
        //temp_zeta[ndx] += -0.5*(tanh((z-10.)))+0.5;
        //temp_zeta[ndx] = max(temp_zeta[ndx],0);
        
    }
}

