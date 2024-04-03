CMAKE_MINIMUM_REQUIRED(VERSION 3.10 FATAL_ERROR)
PROJECT(mizuroute_openwq LANGUAGES CXX Fortran)

if(NOT CMAKE_BUILD_TYPE)
    set(COMPILE_MODE fast)
else()
    SET(COMPILE_MODE ${CMAKE_BUILD_TYPE})
endif()

# Declare use of multiple Languages
ENABLE_LANGUAGE(C)
INCLUDE(FortranCInterface)
FortranCInterface_VERIFY(CXX)
SET(CMAKE_CXX_STANDARD 17)

find_package(HDF5 REQUIRED)
find_package(Boost REQUIRED)
find_package(LAPACK REQUIRED)

# Compile openWQ
SET(exec_name "mizuroute_openwq_${COMPILE_MODE}")
# Default executable directory (Change with -DEXEC_DIR=...)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../../bin)
MESSAGE(STATUS "Executable directory: ${EXEC_DIR}")

SET(OPENWQ_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../openwq/openwq)
SET(OPENWQ_INTERFACE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../openwq)

# File definitions for subcomponents - For OpenWQ
FILE(GLOB OPENWQ_SOURCES 
    ${OPENWQ_DIR}/src/*.cpp 
    ${OPENWQ_DIR}/src/global/*.cpp
    ${OPENWQ_DIR}/src/extwatflux_ss/*.cpp
    ${OPENWQ_DIR}/src/initiate/*.cpp
    ${OPENWQ_DIR}/src/chem/*.cpp
    ${OPENWQ_DIR}/src/watertransp/*.cpp
    ${OPENWQ_DIR}/src/readjson/*.cpp
    ${OPENWQ_DIR}/src/couplercalls/*.cpp
    ${OPENWQ_DIR}/src/output/*.cpp
    ${OPENWQ_DIR}/src/solver/*.cpp
    ${OPENWQ_DIR}/src/units/*.cpp
    ${OPENWQ_DIR}/src/utils/*.cpp
    ${OPENWQ_INTERFACE_DIR}/OpenWQ_hydrolink.cpp
    ${OPENWQ_INTERFACE_DIR}/OpenWQ_interface.cpp)
SET(OPENWQ_INCLUDES ${OPENWQ_DIR}/src/)

IF("${COMPILE_MODE}" STREQUAL "debug")
    SET(CMAKE_Fortran_FLAGS "-g3 -Wall -fmax-errors=0 -fbacktrace -fcheck=all -ffpe-trap=zero -ffree-line-length-none -lstdc++")
    SET(CMAKE_CXX_FLAGS "-g3 -O0 -Wall -pedantic")
ELSEIF("${COMPILE_MODE}" STREQUAL "fast")
    SET(CMAKE_Fortran_FLAGS "-O3 -fmax-errors=0 -ffree-line-length-none -lstdc++ -fopenmp")
    SET(CMAKE_CXX_FLAGS "-O3 ${OpenMP_CXX_FLAGS}")
    FIND_PACKAGE(OpenMP)
ENDIF()


SET(MIZUROUTE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../src)
# Setting mizuroute
SET(MIZUROUTE_UTILS
    ${MIZUROUTE_DIR}/nrtype.f90
    ${MIZUROUTE_DIR}/nr_utility.f90
    ${MIZUROUTE_DIR}/ascii_util.f90
    ${MIZUROUTE_DIR}/ncio_utils.f90
    ${MIZUROUTE_DIR}/gamma_func.f90
)
SET(MIZUROUTE_DATATYPES
    ${OPENWQ_INTERFACE_DIR}/openWQ.f90 # needs to be compiled before globalData
    ${MIZUROUTE_DIR}/public_var.f90
    ${MIZUROUTE_DIR}/dataTypes.f90
    ${MIZUROUTE_DIR}/var_lookup.f90
    ${MIZUROUTE_DIR}/time_utils.f90
    ${MIZUROUTE_DIR}/datetime_data.f90
    ${MIZUROUTE_DIR}/globalData.f90
    ${MIZUROUTE_DIR}/popMetadat.f90
    ${MIZUROUTE_DIR}/allocation.f90
)
SET(MIZUROUTE_INIT
    ${MIZUROUTE_DIR}/network_topo.f90
    ${MIZUROUTE_DIR}/process_param.f90
    ${MIZUROUTE_DIR}/process_ntopo.f90
    ${MIZUROUTE_DIR}/pfafstetter.f90
    ${MIZUROUTE_DIR}/domain_decomposition.f90
)
SET(MIZUROUTE_IO
    ${MIZUROUTE_DIR}/remap.f90
    ${MIZUROUTE_DIR}/read_control.f90
    ${MIZUROUTE_DIR}/read_param.f90
    ${MIZUROUTE_DIR}/read_streamSeg.f90
    ${MIZUROUTE_DIR}/write_streamSeg.f90
    ${MIZUROUTE_DIR}/read_runoff.f90
    ${MIZUROUTE_DIR}/get_basin_runoff.f90
    ${MIZUROUTE_DIR}/read_remap.f90
    ${MIZUROUTE_DIR}/read_restart.f90
    ${MIZUROUTE_DIR}/write_restart.f90
    ${MIZUROUTE_DIR}/write_simoutput.f90
)
SET(MIZUROUTE_CORE
    ${MIZUROUTE_DIR}/model_finalize.f90
    ${MIZUROUTE_DIR}/accum_runoff.f90
    ${MIZUROUTE_DIR}/basinUH.f90
    ${MIZUROUTE_DIR}/irf_route.f90
    ${MIZUROUTE_DIR}/kwt_route.f90
    ${MIZUROUTE_DIR}/dfw_route.f90
    ${MIZUROUTE_DIR}/kw_route.f90
    ${MIZUROUTE_DIR}/mc_route.f90
    ${MIZUROUTE_DIR}/main_route.f90
    ${MIZUROUTE_DIR}/model_setup.f90
    ${OPENWQ_INTERFACE_DIR}/mizuroute_openwq.f90
)
SET(MIZUROUTE_DRIVER 
    ${MIZUROUTE_DIR}/route_runoff.f90
)
# Stitch together MIZUROUTE programs
SET(MIZUROUTE_MOD_ALL
    ${MIZUROUTE_UTILS}
    ${MIZUROUTE_DATATYPES}
    ${MIZUROUTE_INIT}
    ${MIZUROUTE_IO}
    ${MIZUROUTE_CORE})

# compile openWQ into an object file
ADD_LIBRARY(openWQ OBJECT ${OPENWQ_SOURCES})
TARGET_INCLUDE_DIRECTORIES(openWQ PRIVATE 
    $ENV{ARMA_INCLUDES} ${HDF5_INCLUDE_DIRS} ${OPENWQ_INCLUDES})
TARGET_LINK_LIBRARIES(openWQ PUBLIC $ENV{ARMA_LIB} ${HDF5_C_LIBRARY_hdf5})

IF(OpenMP_CXX_FOUND)
    MESSAGE("** OpenWQ **: OpenMP enabled!")
    TARGET_LINK_LIBRARIES(openWQ PUBLIC OpenMP::OpenMP_CXX)
ENDIF()


ADD_EXECUTABLE(${exec_name} ${MIZUROUTE_MOD_ALL} ${MIZUROUTE_DRIVER})
SET_PROPERTY(TARGET ${exec_name} PROPERTY LINKER_LANGUAGE Fortran)
TARGET_INCLUDE_DIRECTORIES(${exec_name} PUBLIC
    ${LAPACK_INCLUDES} $ENV{EBROOTNETCDFMINFORTRAN}/include)
TARGET_LINK_LIBRARIES(${exec_name} PUBLIC -lnetcdff ${LAPACK_LIBRARIES} 
    openWQ)
