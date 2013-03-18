# So far build is with cmake.
# run run_cython.sh to generate all interface files
cmake_minimum_required(VERSION 2.8)
find_package(PythonLibs)
message("sources ${sources}")

include_directories(${PYTHON_INCLUDE_PATH})

enable_language(CXX)
enable_language(Fortran)

add_library(runmc SHARED runmc.c mc.c lj.c ${PYTHON_LIBS})
target_link_libraries(runmc  -lm -lgsl -lgslcblas)
set_target_properties(runmc PROPERTIES PREFIX "")
