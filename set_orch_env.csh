#
# Copy path to external libs here (or comment and set environment variables otherwise)
#
# module load ... # load some modules here if needed on your system
#
# "lib" or "include" end dummy path if required at end of path
setenv MPI_LIB /path/to/mpi/lib
setenv BOOSTPATH /path/to/boost/include
setenv HDF5_INC /path/to/hdf5/include
setenv HDF5_LIB /path/to/hdf5/lib
setenv SUNDIALS_INC /path/to/sundials/include
setenv SUNDIALS_LIB /path/to/sundials/lib
setenv EIGEN_INC /path/to/eigen
setenv TENSORFLOW_CAPI_PATH /path/to/tensorflow
setenv OPENCV_PATH /path/to/opencv
#
# Do not touch
#
set script_path = `ls -l /proc/$$/fd | sed -e 's/^[^/]*//' | grep "/set_orch_env\.csh" | sed 's/set_orch_env\.csh//g'` 
echo $script_path
setenv ORCH_BASE $script_path/
setenv ORCH_MASTER ${ORCH_BASE}
setenv ORCH ${ORCH_MASTER}ORCh
setenv ORCH_TEST ${ORCH_MASTER}Tests
setenv GTCOMB_CT_HOME ${ORCH_BASE}Cantera
setenv GTCOMB_CT_HOSTTYPE ${GTCOMB_CT_HOME}/user/linux
setenv GTCOMB_CT_DATA ${GTCOMB_CT_HOME}/data
#
sed -i "/^prefix =/c\prefix = '"${ORCH_BASE}"Cantera_lib'" ${GTCOMB_CT_HOME}/cantera.conf
sed -i "/^python_prefix =/c\python_prefix = '"${ORCH_BASE}"dummy_python'" ${GTCOMB_CT_HOME}/cantera.conf
sed -i "/^python3_prefix =/c\python3_prefix = '"${ORCH_BASE}"dummy_python'" ${GTCOMB_CT_HOME}/cantera.conf
sed -i "/^sundials_include =/c\sundials_include = '"${SUNDIALS_INC}"'" ${GTCOMB_CT_HOME}/cantera.conf
sed -i "/^sundials_libdir =/c\sundials_libdir = '"${SUNDIALS_LIB}"'" ${GTCOMB_CT_HOME}/cantera.conf
sed -i "/^boost_inc_dir =/c\boost_inc_dir = '"${BOOSTPATH}"'" ${GTCOMB_CT_HOME}/cantera.conf
sed -i "/^extra_inc_dirs =/c\extra_inc_dirs = '"${EIGEN_INC}"'" ${GTCOMB_CT_HOME}/cantera.conf
#
setenv CT_INC ${ORCH_BASE}/Cantera_lib/include
setenv CT_LIB ${ORCH_BASE}/Cantera_lib/lib
#
setenv LD_LIBRARY_PATH ${ORCH_BASE}:${LD_LIBRARY_PATH}
