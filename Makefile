fftw = /Users/yang/Downloads/fftw-3.3.6-pl2/

default:
	mpicc -Wall -O2 -c aniso_interfacial.c
	mpicc -Wall -O2 -c one_vec_phi.c -I$(fftw)/include 
	mpicc -Wall -o one_vec_phi one_vec_phi.o aniso_interfacial.o -L/usr/local/lib -lfftw3 -lfftw3_mpi

