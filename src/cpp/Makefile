#---------------------------#
#DO_MPI= no
DO_MPI= yes

With_Elemental=yes

ifeq ($(DO_MPI), yes)
  CC= mpicc
  CXX= mpic++
else
  CC= gcc
  CXX= g++
endif
#---------------------------#

ifeq ($(With_Elemental),yes)
  include ElVars_GPC
endif



# C main() call fortran subroutine, for ifort.
COMFFLAGS = -nofor_main    
FLIB = -L/usr/local/lib  -lm   




#-----------------------------------------------------------------#
#------------------------- My NOTEBOOKS---------------------------#
#-----------------------------------------------------------------#

ifeq ($(MACHINE),mypro13)
    MYLIBPATH= /Users/wangxin/workspace/lib/mylib
    INIPATH = /Users/wangxin/workspace/lib/iniparser
    MPIPATH = /Users/wangxin/software/work/mpich2/mpich2-1.2.1p1

    MPIINCLUDE= $(MPIPATH)/include

    GENLIB   = lmylib
    FISHLIB  = lfish

    GSLIB = -L/Users/wangxin/workspace/lib/gsl/lib -lgsl -lgslcblas -lm
    MYLIB = -L$(MYLIBPATH)  -$(GENLIB) -$(FISHLIB) 

    #ifeq ($(DO_MPI), yes)
    #    CFLAGS = -g -arch x86_64 -D_MPI_ -I /Users/wangxin/workspace/lib/gsl/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
    #else
    #    CFLAGS = -g -arch i386 -arch x86_64 -I /Users/wangxin/workspace/lib/gsl/include \
    #    	 -I $(MYLIBPATH) -I $(INIPATH)  

    #endif

    ifeq ($(DO_MPI), yes)
        CFLAGS = -g  -D_MPI_ -I /Users/wangxin/workspace/lib/gsl/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
    else
        CFLAGS = -g -I /Users/wangxin/workspace/lib/gsl/include \
        	 -I $(MYLIBPATH) -I $(INIPATH)  

    endif


endif





#---- CITA VOLE ----#
ifeq ($(MACHINE),cita)
    MYLIBPATH= /home/xwang/workspace/lib/mylib
    INIPATH = /home/xwang/workspace/lib/iniparser
    #MPIPATH = 
    #FFTWPATH = 
    CUBAPATH=/home/xwang/workspace/software/cuba/current

    MPIINCLUDE= $(MPIPATH)/include
    #FFTWINCLUDE = $(FFTWPATH)/

    CUBAINCLUDE=$(CUBAPATH)/include

    GENLIB  = lmylib
    FISHLIB = lfish
    CUBALIB = lcuba

    #FFTWLIB = -lfftw3 -lm
    #FFTWLIB = -lfftw3f -lfftw3 -lm

    GSLIB = -L/home/xwang/workspace/lib/gsl/lib -lgsl -lgslcblas -lm
    MYLIB = -L$(MYLIBPATH)  -$(GENLIB) -$(FISHLIB)  $(FFTWLIB) -L$(CUBAPATH)/lib -$(CUBALIB)


    #ifeq ($(DO_OMP),yes)
    #    CFLAGS = -fopenmp -g  -D _OMP_ -I /home/xwang/workspace/lib/gsl/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
    #    FFTWLIB =  -lfftw3f_threads -lfftw3_threads -lfftw3f -lfftw3 -lm

    #else
    #    CFLAGS = -g -I /home/xwang/workspace/lib/gsl/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
    #    FFTWLIB = -lfftw3f -lfftw3 -lm
    #endif


    ifeq ($(DO_MPI), yes)
	CFLAGS = -g -D _MPI_ -I /home/xwang/workspace/lib/gsl/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)   -I $(CUBAINCLUDE)
    else
	CFLAGS = -g -I /home/xwang/workspace/lib/gsl/include -I $(MYLIBPATH) -I $(INIPATH)   -I $(CUBAINCLUDE)
    endif


endif


#---
GSLPATH=/home/p/pen/wangxin/workspace/lib
CFLAGS= -g -D _MPI_ -I $(GSLPATH)/include
GSLIB= -L$(GSLPATH)/lib -lgsl -lgslcblas -lm



#################################################################
#LIB  = $(GSLIB) $(MYLIB)




.SUFFIXES : .o .c .cpp
.c.o :
	$(CC) $(CFLAGS) $(EL_COMPILE_FLAGS) -c $<

.f.o :
	$(FC) $(FFLAGS) -c $<

.cpp.o :
	$(CXX) $(CFLAGS)  $(EL_COMPILE_FLAGS) -c $<




# workspace path
WORKSPACE = ../../workspace
AR=ar
ARFLAGS=-r
RANLIB=ranlib


default: quadestxx


#--------------------------------------------------#
#quadestxx: io.o mpinit.o quadest.o driver.o
quadestxx: io.o mpinit.o cov.o quadest_init.o quadest.o driver.o
	$(CXX) $(EL_LINK_FLAGS)  $^ -o $(WORKSPACE)/$@  $(EL_LIBS) $(GSLIB)



clean:
	rm -f *.o *.c~




