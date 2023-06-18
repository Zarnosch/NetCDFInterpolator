# At LEAST 3.2 but newer is better
cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(GridInterpolation VERSION 0.1 LANGUAGES CXX)

# find_package(Matplot++ REQUIRED)

# Set the path where cmake modules are stored
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)
message("Module path: " ${CMAKE_MODULE_PATH})

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
    message("CCache found")
endif()


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)


find_package(NetCDF REQUIRED)

find_package(OpenMP)

find_package (Eigen3 3.3 REQUIRED)

find_package(Boost REQUIRED COMPONENTS 
		iostreams python38 filesystem
)

message(STATUS ${NETCDF_INCLUDES})
message(STATUS ${NETCDF_LIBRARIES})


# Define an executable
add_executable(grid_interpol
               src/main.cc
               src/util.cc
               )
               
add_library(grid_interpol_lib SHARED
       	src/main.cc
		src/util.cc
       	)

target_include_directories(
    grid_interpol
    PRIVATE ./inc
    PRIVATE ./
    PRIVATE ${NETCDF_INCLUDES}
    /usr/local/inc
    /usr/local/include
    /usr/local/include/eigen3
    )
    
target_include_directories(
    grid_interpol_lib
    PRIVATE ./inc
    PRIVATE ./
    PRIVATE ${NETCDF_INCLUDES}
    /usr/local/inc
    /usr/local/include
    /usr/local/include/eigen3
    )


# dependencies for libraries
target_link_libraries(
    grid_interpol
    PUBLIC ${Eigen_INCLUDE_DIRS}
    /usr/local/lib/libnetcdf-cxx4.so
    )
    
# dependencies for libraries
target_link_libraries(
    grid_interpol_lib
    PUBLIC ${Eigen_INCLUDE_DIRS}
    /usr/local/lib/libnetcdf-cxx4.so
    )
    
# additional compiler flags for more warnings
list(APPEND CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic -std=c++2a ${OpenMP_CXX_FLAGS}") # 
set(CMAKE_CXX_FLAGS_DEBUG "-g3 -llz4")
set(CMAKE_CXX_FLAGS_RELEASE "-o3 -llz4") #-o3
set(CMAKE_EXE_LINKER_FLAGS "-L/usr/local -lz -lm")



