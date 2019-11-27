# Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.1)

# Project's name
project(xg)
# We build using c++14
set(CMAKE_CXX_STANDARD 14)

# Use all standard-compliant optimizations
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g")

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

  # assumes clang build
  # we can't reliably detect when we're using clang, so for the time being we assume
  # TODO: can't we though?
  
  # adapted from https://stackoverflow.com/questions/46414660/macos-cmake-and-openmp
  # find_package(OpenMP) does not work reliably on macOS, so we do its work ourselves
  set (OpenMP_C "${CMAKE_C_COMPILER}")
  set (OpenMP_C_FLAGS " -Xpreprocessor -fopenmp -I/opt/local/include/libomp -I/usr/local/include -L/opt/local/lib/libomp -L/usr/local/lib")
  set (OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
  set (OpenMP_CXX "${CMAKE_CXX_COMPILER}")
  set (OpenMP_CXX_FLAGS " -Xpreprocessor -fopenmp -I/opt/local/include/libomp -I/usr/local/include -L/opt/local/lib/libomp -L/usr/local/lib")
  set (OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
  set (OpenMP_libomp_LIBRARY "omp")
  set (OpenMP_libgomp_LIBRARY "gomp")
  set (OpenMP_libiomp5_LIBRARY "iomp5")
  
  # and now add the OpenMP parameters to the compile flags
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} -lomp")
  
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")

  find_package(OpenMP REQUIRED)
  
  # add the flags it detects to the compile flags
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fopenmp")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  
endif()

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

# The following folder will be included
include_directories("${PROJECT_SOURCE_DIR}")

# Add external projects
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

# TODO: We're using INSTALL_DIR very wrong. We *should* be actually installing
# the external projects into their prefixes and working with the installed
# files. Instead we're building but not installing them and trying to work with
# the non-installed build trees.
# 
# Hence the blanked out INSTALL_COMMANDs to suppress the install step.
#
# By blanking out UPDATE_COMMAND, we have to blow away (rm -rf) the build directories
# for projects when we update their git tag.
# The alternative is that we have to reconfigure/rebuild certain repos every single time.

# libhandlegraph (full build using its cmake config)
ExternalProject_Add(handlegraph
  GIT_REPOSITORY "https://github.com/vgteam/libhandlegraph.git"
  GIT_TAG "7b788bfaea21d976e372b130a2353deff99757f1"
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR}
  #BUILD_COMMAND ""
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(handlegraph INSTALL_DIR)
set(handlegraph_INCLUDE "${INSTALL_DIR}/src/handlegraph/src/include")
set(handlegraph_LIB "${INSTALL_DIR}/src/handlegraph-build")

# sdsl-lite (full build using its cmake config)
ExternalProject_Add(sdsl-lite
  GIT_REPOSITORY "https://github.com/simongog/sdsl-lite.git"
  GIT_TAG "ddb0fbbc33bb183baa616f17eb48e261ac2a3672"
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR} # TODO ADD static build flag
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(sdsl-lite INSTALL_DIR)
set(sdsl-lite_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/include")
set(sdsl-lite-divsufsort_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/include")
set(sdsl-lite_LIB "${INSTALL_DIR}/src/sdsl-lite-build/lib")
set(sdsl-lite-divsufsort_LIB "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/lib")

# mmmultimap (memory mapped multimap)
ExternalProject_Add(mmmultimap
  GIT_REPOSITORY "https://github.com/ekg/mmmultimap.git"
  GIT_TAG "38d5d70f4f9e10aae4548b022c2a7d55024c5e99"
  BUILD_COMMAND ""
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(mmmultimap SOURCE_DIR)
set(mmmultimap_INCLUDE "${SOURCE_DIR}/src")

# taywee's C++ args library, header only
ExternalProject_Add(tayweeargs
  GIT_REPOSITORY "https://github.com/Taywee/args.git"
  GIT_TAG "3de44ec671db452cc0c4ef86399b108939768abb"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(tayweeargs SOURCE_DIR)
set(tayweeargs_INCLUDE "${SOURCE_DIR}")

# gfakluge (header only)
ExternalProject_Add(gfakluge
  GIT_REPOSITORY "https://github.com/vgteam/gfakluge.git"
  GIT_TAG "3b9faa69fc27f34b9f65d78dd66bc4c97f11014b"
  BUILD_IN_SOURCE TRUE
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR} # TODO ADD static build flag
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(gfakluge INSTALL_DIR)
set(gfakluge_INCLUDE "${INSTALL_DIR}/src/gfakluge/src")
set(gfakluge_tinyFA_INCLUDE "${INSTALL_DIR}/src/gfakluge/src/tinyFA")

# In-place Parallel Super Scalar Samplesort (IPS⁴o), header only
ExternalProject_Add(ips4o
  GIT_REPOSITORY "https://github.com/vgteam/ips4o.git"
  GIT_TAG "22069381cc1bf2df07ee1ff47f6b6073fcfb4508"
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(ips4o SOURCE_DIR)
set(ips4o_INCLUDE "${SOURCE_DIR}")

#set(CMAKE_BUILD_TYPE Release)

add_library(xg_objs OBJECT src/xg.cpp)

add_dependencies(xg_objs handlegraph)
add_dependencies(xg_objs sdsl-lite)
add_dependencies(xg_objs mmmultimap)
add_dependencies(xg_objs tayweeargs)
add_dependencies(xg_objs gfakluge)
add_dependencies(xg_objs ips4o)

target_include_directories(xg_objs PUBLIC
  "${CMAKE_SOURCE_DIR}/src"
  "${handlegraph_INCLUDE}"
  "${sdsl-lite_INCLUDE}"
  "${sdsl-lite-divsufsort_INCLUDE}"
  "${tayweeargs_INCLUDE}"
  "${ips4o_INCLUDE}"
  "${mmmultimap_INCLUDE}"
  "${gfakluge_INCLUDE}"
  "${gfakluge_tinyFA_INCLUDE}")

set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

# The following folder will be included
include_directories("${PROJECT_SOURCE_DIR}")

# build a static library
#add_library(libxg STATIC ${CMAKE_SOURCE_DIR}/src/xg.cpp)
add_library(libxg STATIC $<TARGET_OBJECTS:xg_objs>)
set_target_properties(libxg PROPERTIES OUTPUT_NAME "xg")

# set up our target executable and specify its dependencies and includes
add_executable(xg
  $<TARGET_OBJECTS:xg_objs>
  ${CMAKE_SOURCE_DIR}/src/main.cpp
  )
target_include_directories(xg PUBLIC
  "${CMAKE_SOURCE_DIR}/src"
  "${handlegraph_INCLUDE}"
  "${sdsl-lite_INCLUDE}"
  "${sdsl-lite-divsufsort_INCLUDE}"
  "${tayweeargs_INCLUDE}"
  "${ips4o_INCLUDE}"
  "${mmmultimap_INCLUDE}"
  "${gfakluge_INCLUDE}"
  "${gfakluge_tinyFA_INCLUDE}")


# macOS doesn't want you to link in libatomic this way
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  target_link_libraries(xg
    "${sdsl-lite_LIB}/libsdsl.a"
    "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
    "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a"
    "${handlegraph_LIB}/libhandlegraph.a"
    )
elseif (TRUE)
  target_link_libraries(xg
    "${sdsl-lite_LIB}/libsdsl.a"
    "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
    "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a"
    "${handlegraph_LIB}/libhandlegraph.a"
    "-latomic"
    )
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

