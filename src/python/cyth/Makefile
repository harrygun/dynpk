ifeq ($(MACHINE),mylaptop)
    cflag= "-arch i386 -arch x86_64 -I /Users/harrygun/workspace/software/works/cuba/include -I/Users/harrygun/workspace/lib/gsl/include -I../../clib/src/c -I /Library/Frameworks/Python.framework/Versions/Current/lib/python2.7/site-packages/mpi4py/include"
    libflag= "-L /Users/harrygun/workspace/software/works/cuba/lib -L /Users/harrygun/workspace/lib/mylib/script/clib -L ./  -L/Users/harrygun/workspace/lib/gsl/lib -lgsl -lgslcblas -lm -L../../clib/lib -lmysf"

    python='python'
endif


ifeq ($(MACHINE),myair)
#    cflag= "-arch i386 -arch x86_64 -I /Users/wangxin/workspace/software/works/cuba/include -I /Users/wangxin/workspace/lib/gsl/include -I../../clib/src/c -I /Library/Frameworks/Python.framework/Versions/Current/lib/python2.7/site-packages/mpi4py/include"
    cflag= "-I /Users/wangxin/workspace/software/works/cuba/include -I /Users/wangxin/workspace/lib/gsl/include -I../../clib/src/c -I /Library/Frameworks/Python.framework/Versions/Current/lib/python2.7/site-packages/mpi4py/include"
    libflag= "-L /Users/wangxin/workspace/software/works/cuba/lib -L /Users/wangxin/workspace/lib/mylib/python/clib -L ./  -L/Users/wangxin/workspace/lib/gsl/lib -lgsl -lgslcblas -lm -L ../../clib/lib -lmysf"

    python='python'
endif


ifeq ($(MACHINE),mypro13)
    cflag= "-I /Users/wangxin/workspace/software/works/cuba/include -I /Users/wangxin/workspace/lib/gsl/include -I../../clib/src/c -I /Library/Frameworks/Python.framework/Versions/Current/lib/python2.7/site-packages/mpi4py/include"
    libflag= "-L /Users/wangxin/workspace/software/cuba/lib -L /Users/wangxin/workspace/lib/mylib/python/clib -L ./  -L/Users/wangxin/workspace/lib/gsl/lib -lgsl -lgslcblas -lm "

    python=python2.7
endif


ifeq ($(MACHINE),gwln)
    cflag="-I /home/wangxin/workspace/software/cuba/include -I /home/wangxin/software/mpich2/gwln1/include -I/home/wangxin/workspace/lib/include  -I ../../clib/src/c"
    libflag="-L /home/wangxin/workspace/software/cuba/lib -L ./ -L/home/wangxin/workspace/lib/lib -lgsl -lgslcblas -lm -L ../../clib/lib -lmysf"

    python='python'
endif


ifeq ($(MACHINE),HHPC)
    cflag="-I/home/wangxin/software/cuba/include -I/usr/mpi/intel/openmpi-1.4.3/include -I/home/wangxin/workspace/lib/include -I/home/wangxin/lib/gsl/include -I ../../clib/src/c"
    libflag="-L/home/wangxin/software/cuba/lib -L ./ -L/home/wangxin/lib/gsl/lib -lgsl -lgslcblas -lm -L ../../clib/lib -lmysf"
    python='python'
endif


ifeq ($(MACHINE),France_bingo)
    cflag="-I /home/xwang/software/cuba/include -I /usr/local/include  -I ../../clib/src/c"
    libflag="-L/home/xwang/software/cuba/lib  -L./ -L/usr/local/lib -lgsl -lgslcblas -lm -L ../../clib/lib -lmysf"

    python='python2.7'
endif

ifeq ($(MACHINE),cita)
    cflag="-fPIC -I /home/xwang/workspace/software/cuba/current/include -I /opt/openmpi/1.6.1-intel-12.1.5/include -I /home/xwang/workspace/lib/gsl/include  -I ../../clib/src/c"
    libflag="-L /home/xwang/workspace/software/cuba/current/lib -L ./ -L /home/xwang/workspace/lib/gsl/lib -lgsl -lgslcblas -lm -L ../../clib/lib -lmysf"

    python='python'
endif


all:
	$(python) -V
	CFLAGS=$(cflag) LDFLAGS=$(libflag) \
	$(python) setup.py build_ext --inplace


clean:
	@echo Cleaning
	@rm -f *.c *.o *.so *~
	@rm -rf build
