from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext, extension
import numpy as np
import os
import cython_gsl


machine=os.environ['MACHINE']
print 'compiling on', machine


if machine=='mypro13':
    mpi4py_include='/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/mpi4py/include' 
    #LIB_PATH=os.environ['LD_LIBRARY_PATH'].split(':')

elif machine=='mylaptop':
    mpi4py_include='/Library/Frameworks/Python.framework/Versions/Current/lib/python2.7/site-packages/mpi4py/include' 
    LIB_PATH=os.environ['LD_LIBRARY_PATH'].split(':')

elif machine=='myair':
    mpi4py_include='/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/mpi4py' 
    LIB_PATH=[' '] #os.environ['LD_LIBRARY_PATH'].split(':')

elif machine=='gwln':
    mpi4py_include='/home/wangxin/software/python/lib/python2.7/site-packages/mpi4py/include'
    LIB_PATH=os.environ['LD_LIBRARY_PATH'].split(':')

elif machine=='HHPC':
    mpi4py_include='/home/wangxin/software/python/python/2.7.3/lib/python2.7/site-packages/mpi4py/include'
    LIB_PATH=os.environ['LD_LIBRARY_PATH'].split(':')


src_dir= './'

#src_dir= '../c/'
#src_c = ['misc.c', 'kernel.c', 'abkernel.c']
#c_lib = #'libwigner.a'

kernel_src=['covm.pyx'] 
print kernel_src

module = [ Extension("covm", sources=[kernel_src[0]], 
           include_dirs= [np.get_include(), mpi4py_include, src_dir, cython_gsl.get_cython_include_dir()], 
           libraries=cython_gsl.get_libraries()+['m'],
           library_dirs=[cython_gsl.get_library_dir(), './', src_dir]),
          ] 


setup( name = 'Quadratic_Estimator_Cython', cmdclass = {'build_ext': build_ext},
       ext_modules = module)



