if [ -f /etc/bashrc ]; then
       . /etc/bashrc
fi

export PATH=$PATH:/scinet/bgq/bin/:/scinet/gpc/tools/elemental/0.86-rc1-1056-g4a97eb8/lib:/scinet/gpc/intel/compilers_and_libraries_2016.3.210/linux/mkl/lib/intel64/:/scinet/gpc/tools/elemental/0.86-rc1-1056-g4a97eb8/lib/:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/scinet/gpc/toolbin:/scinet/gpc/toolbin/vnc:/scinet/gpc/toolbin/x11/bin:/usr/lpp/mmfs/bin:/opt/torque/bin:/opt/torque/sbin

export MKLROOT=/scinet/gpc/intel/compilers_and_libraries_2016.3.210/linux/mkl 
export SCINET_ELEMENTAL_LIB=/scinet/gpc/tools/elemental/0.86-rc1-1056-g4a97eb8/lib 
export SCINET_ELEMENTAL_INC=/scinet/gpc/tools/elemental/0.86-rc1-1056-g4a97eb8/include 
export SCINET_ELEMENTAL_BASE=/scinet/gpc/tools/elemental/0.86-rc1-1056-g4a97eb8 

#module purge 
#module load gcc/5.2.0
#module load openmpi/gcc/1.8.3
#module load use.experimental elemental/0.86-rc1-1056-g4a97eb8

# commands which work for both GPC and TCS can go here

HOST=$(uname)

if [ "${HOST}" == "AIX" ]; then
	# do things for the TCS machine
	# user environment for all shells goes here
	# replace colon with your own commands
	:
else
	# do things for the GPC machine

	# user environment for all shells goes here
	# replace colon with your own commands
	:
fi

