#---------------------------#
DO_MPI= no
#DO_MPI= yes


ifeq ($(DO_MPI), yes)
  CC= mpicc
  #CC= mpic++
else
  CC= gcc
  #CC= g++
endif
MPICC= mpicc
FC= ifort
#---------------------------#



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

    INILIB = -L$(INIPATH) -liniparser
    MPIINCLUDE= $(MPIPATH)/include

    GENLIB   = lmylib
    FISHLIB  = lfish

    GSLIB = -L/Users/wangxin/workspace/lib/gsl/lib -lgsl -lgslcblas -lm
    MYLIB = -L$(MYLIBPATH)  -$(GENLIB) -$(FISHLIB) 

    ifeq ($(DO_MPI), yes)
	CFLAGS = -g -arch x86_64 -D_MPI_ -I /Users/wangxin/workspace/lib/gsl/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
    else
	CFLAGS = -g -arch i386 -arch x86_64 -I /Users/wangxin/workspace/lib/gsl/include \
		 -I $(MYLIBPATH) -I $(INIPATH)  

	#CFLAGS = -g -I /Users/wangxin/workspace/lib/gsl/include \
	#	 -I $(MYLIBPATH) -I $(INIPATH)  -I $(MPIINCLUDE)
    endif

endif


ifeq ($(MACHINE),myair)
    MYLIBPATH= /Users/wangxin/workspace/lib/mylib
    INIPATH = /Users/wangxin/workspace/lib/iniparser
    MPIPATH = /Users/wangxin/software/work/mpich2/mpich2-1.2.1p1

    INILIB = -L$(INIPATH) -liniparser
    MPIINCLUDE= $(MPIPATH)/include

    GENLIB   = lmylib
    FISHLIB  = lfish

    GSLIB = -L/Users/wangxin/workspace/lib/gsl/lib -lgsl -lgslcblas -lm
    MYLIB = -L$(MYLIBPATH)  -$(GENLIB) -$(FISHLIB) 

    ifeq ($(DO_MPI), yes)
	CFLAGS = -g -arch x86_64 -D_MPI_ -I /Users/wangxin/workspace/lib/gsl/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
    else
	CFLAGS = -g -arch i386 -arch x86_64 -I /Users/wangxin/workspace/lib/gsl/include \
		 -I $(MYLIBPATH) -I $(INIPATH)  -I $(MPIINCLUDE)

	#CFLAGS = -g -I /Users/wangxin/workspace/lib/gsl/include \
	#	 -I $(MYLIBPATH) -I $(INIPATH)  -I $(MPIINCLUDE)
    endif

endif


ifeq ($(MACHINE),mylaptop)
    MYLIBPATH= /Users/harrygun/workspace/lib/mylib
    INIPATH = /Users/harrygun/workspace/lib/iniparser
    MPIPATH = /Users/harrygun/software/work/mpich2/mpich2-1.2.1p1

    INILIB = -L$(INIPATH) -liniparser
    MPIINCLUDE= $(MPIPATH)/include

    GENLIB   = lmylib
    FISHLIB  = lfish

    GSLIB = -L/Users/harrygun/workspace/lib/gsl/lib -lgsl -lgslcblas -lm
    MYLIB = -L$(MYLIBPATH)  -$(GENLIB) -$(FISHLIB) 

    ifeq ($(DO_MPI), yes)
	CFLAGS = -g -arch x86_64 -D_MPI_ -I /Users/harrygun/workspace/lib/gsl/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
    else
	CFLAGS = -g -arch i386 -arch x86_64 -I /Users/harrygun/workspace/lib/gsl/include \
		 -I $(MYLIBPATH) -I $(INIPATH)  -I $(MPIINCLUDE)

	#CFLAGS = -g -I /Users/harrygun/workspace/lib/gsl/include \
	#	 -I $(MYLIBPATH) -I $(INIPATH)  -I $(MPIINCLUDE)
    endif

endif

#---- gwln ----
ifeq ($(MACHINE),gwln)
    MYLIBPATH= /home/wangxin/workspace/lib/mylib
    INIPATH = /home/wangxin/workspace/lib/iniparser/$(HOSTNAME)
    MPIPATH = /home/wangxin/software/mpich2/gwln1

    INILIB = -L$(INIPATH) -liniparser
    MPIINCLUDE= $(MPIPATH)/include

    GENLIB   = lmylib
    FISHLIB  = lfish

    GSLIB = -L /home/wangxin/workspace/lib/gsl/$(HOSTNAME)/lib -lgsl -lgslcblas -lm
    MYLIB = -L$(MYLIBPATH)/$(HOSTNAME)  -$(GENLIB) -$(FISHLIB) 

    ifeq ($(DO_MPI), yes)
	CFLAGS = -g -D_MPI_ -I /home/wangxin/workspace/lib/gsl/$(HOSTNAME)/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
    else
	CFLAGS = -g -fPIC -I /home/wangxin/workspace/lib/gsl/$(HOSTNAME)/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
    endif

endif

#---- HHPC ----
ifeq ($(MACHINE),HHPC)
    MYLIBPATH= /home/wangxin/lib/mylib
    INIPATH = /home/wangxin/lib/iniparser
    MPIPATH = /usr/mpi/intel/openmpi-1.4.3

    INILIB = -L$(INIPATH) -liniparser
    MPIINCLUDE= $(MPIPATH)/include

    GENLIB   = lmylib
    FISHLIB  = lfish

    GSLIB = -L/home/wangxin/lib/gsl/lib -lgsl -lgslcblas -lm
    MYLIB = -L$(MYLIBPATH)  -$(GENLIB) -$(FISHLIB) 


    ifeq ($(DO_MPI), yes)
	CFLAGS = -g -D_MPI_ -I /home/wangxin/lib/gsl/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
    else
	CFLAGS = -g -I /home/wangxin/lib/gsl/include -I $(MYLIBPATH) -I $(INIPATH) -I $(MPIINCLUDE)  
    endif

endif


#---- CITA VOLE ----#
ifeq ($(MACHINE),cita)
    MYLIBPATH= /home/xwang/workspace/lib/mylib
    INIPATH = /home/xwang/workspace/lib/iniparser
    #MPIPATH = 
    #FFTWPATH = 
    CUBAPATH=/home/xwang/workspace/software/cuba/current

    INILIB = -L$(INIPATH) -liniparser
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







#################################################################
LIB  = $(INILIB) $(GSLIB) $(MYLIB)




.SUFFIXES : .o .c
.c.o :
	$(CC) $(CFLAGS) -c $<

.f.o :
	$(FC) $(FFLAGS) -c $<





# workspace path
WORKSPACE = ../../workspace
AR=ar
ARFLAGS=-r
RANLIB=ranlib


default: quadest


#--------------------------------------------------#
quadest: quadest.o io.o mpinit.o misc.o driver.o
	$(CC)  $^ -o $(WORKSPACE)/$@  $(LIB)



clean:
	rm -f *.o *.c~




