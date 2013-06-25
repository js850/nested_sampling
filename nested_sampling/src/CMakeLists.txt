# So far build is with cmake.
# run run_cython.sh to generate all interface files
cmake_minimum_required(VERSION 2.8)
find_package(PythonLibs)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})
#GSL
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GSL_LIBRARIES})
message("sources ${sources}")

include_directories(${PYTHON_INCLUDE_PATH})

enable_language(C)
enable_language(Fortran)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")

add_library(runmc SHARED runmc.c mc.c lj.c ${PYTHON_LIBS})
target_link_libraries(runmc -L/usr/lib -lgsl -lgslcblas -lm)
set_target_properties(runmc PROPERTIES PREFIX "")

add_library(run_ising_mc SHARED run_ising_mc.c mcising.c ${PYTHON_LIBS})
target_link_libraries(run_ising_mc -L/usr/lib -lgsl -lgslcblas -lm)
set_target_properties(run_ising_mc PROPERTIES PREFIX "")

add_library(cv_trapezoidal SHARED cv_trapezoidal.c cv.c ${PYTHON_LIBS})
target_link_libraries(cv_trapezoidal -L/usr/lib -lgsl -lgslcblas -lm)
set_target_properties(cv_trapezoidal PROPERTIES PREFIX "")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
#add_executable(LJ_ns_plot LJ_ns_plot.cpp)
#install (TARGETS LJ_ns_plot DESTINATION bin)
