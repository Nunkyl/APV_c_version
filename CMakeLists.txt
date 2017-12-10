cmake_minimum_required(VERSION 3.7)
project(APV_c_version)

set(CMAKE_C_STANDARD 11)

# added gsl
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lgsl -lgslcblas")
link_directories(/usr/local/Cellar/gsl) # headers
link_libraries(/usr/local/Cellar/gsl/2.4/lib/libgsl.a)


set(SOURCE_FILES main.c PValue.c PValue.h)
add_executable(APV_c_version ${SOURCE_FILES} PValue.c PValue.h)