all: transporter

build:
	make C=ifort

debug:
	make C=ifortDEBUG

zeus:
	make C=ifortZEUS

ifeq ($(C),ifort)
FC=ifort
BASEDIR=/home/mkk/libs
FCFLAGS= -c -O3  -132  -I$(BASEDIR)/SuperLU_4.3/SRC
FCCFLAGS= -c -O3 -I$(BASEDIR)/SuperLU_4.3/SRC
FBFLAGS=  -O3  -132
FLFLAGS= -L/opt/intel/mkl/lib/intel64/ -L/opt/intel/composer_xe_2013_sp1.1.106/mkl/lib/intel64/
#LIBS= $(BASEDIR)/lapack-3.4.2/liblapack.a  $(BASEDIR)/SuperLU_4.3/lib/libsuperlu_4.3.a  $(BASEDIR)/BLAS/blas_IFORT.a
LIBS=  $(BASEDIR)/libsuperlu_4.3.a -mkl -static-intel

else ifeq ($(C),ifortDEBUG)
FC=ifort
BASEDIR=/home/mkk/libs
FCFLAGS= -c -132 -traceback -O0 -check all -fpe0 -warn -traceback -debug extended -I$(BASEDIR)/SuperLU_4.3/SRC
FCCFLAGS= -c -O0 -Wall -g -I$(BASEDIR)/SuperLU_4.3/SRC
FBFLAGS=  -O0 -132
#LIBS=  $(BASEDIR)/lapack-3.4.2/liblapack.a  $(BASEDIR)/SuperLU_4.3/lib/libsuperlu_4.3.a  $(BASEDIR)/BLAS/blas_IFORT.a  -lm
LIBS=  $(BASEDIR)/libsuperlu_4.3.a -mkl -static-intel
else ifeq ($(C),ifortZEUS)
BASEDIR=/people/gjkolasi
FC=ifort
FCFLAGS= -c -O3  -132  -I$(BASEDIR)/SuperLU_4.3/SRC
FCCFLAGS= -c -O3  -I$(BASEDIR)/SuperLU_4.3/SRC
FBFLAGS=  -O3  -132
LIBS= -mkl $(BASEDIR)/libsuperlu_4.3.a
FLFLAGS=
endif



transporter: main.f90 modutils.o modinip.o modjed.o modpop.o spinmodpop.o modspinzrodlo.o modzrodlo.o zgssv.o modsystem.o spinmodsystem.o
	$(FC) $(FBFLAGS)  main.f90 *.o $(LIBS)   -o $@
modinip.o: modinip.F90
	$(FC) $(FCFLAGS) modinip.F90 -o $@
modutils.o: modutils.F90
	$(FC) $(FCFLAGS) modutils.F90 -o $@
modjed.o: modjed.F90
	$(FC) $(FCFLAGS) modjed.F90 -o $@
modpop.o: modpop.F90
	$(FC) $(FCFLAGS) modpop.F90 -o $@
spinmodpop.o: spinmodpop.f90
	$(FC) $(FCFLAGS) spinmodpop.f90 -o $@

modzrodlo.o: modzrodlo.f90
	$(FC) $(FCFLAGS) modzrodlo.f90 -o $@
modspinzrodlo.o: modspinzrodlo.f90
	$(FC) $(FCFLAGS) modspinzrodlo.f90 -o $@

modsystem.o: modsystem.f90
	$(FC) $(FCFLAGS) modsystem.f90 -o $@
spinmodsystem.o: spinmodsystem.f90
	$(FC) $(FCFLAGS) spinmodsystem.f90 -o $@
zgssv.o:c_fortran_zgssv.c
	gcc   $(FCCFLAGS) c_fortran_zgssv.c -o $@



clean:
	rm *.o *.mod
