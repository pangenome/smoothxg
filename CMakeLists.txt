# Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.1)

# Project's name
project(smoothxg)
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

# libhandlegraph (full build using its cmake config)
ExternalProject_Add(handlegraph
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/libhandlegraph"
  CMAKE_ARGS "${CMAKE_ARGS};-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>")
ExternalProject_Get_property(handlegraph INSTALL_DIR)
set(handlegraph_INCLUDE "${INSTALL_DIR}/include")
set(handlegraph_LIB "${INSTALL_DIR}/lib")

# sdsl-lite (full build using its cmake config)
ExternalProject_Add(sdsl-lite
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/sdsl-lite"
  CMAKE_ARGS "${CMAKE_ARGS};-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(sdsl-lite INSTALL_DIR)
set(sdsl-lite_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/include")
set(sdsl-lite-divsufsort_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/include")
set(sdsl-lite_LIB "${INSTALL_DIR}/src/sdsl-lite-build/lib")
set(sdsl-lite-divsufsort_LIB "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/lib")

# odgi
ExternalProject_Add(odgi
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/odgi"
  UPDATE_COMMAND ""
  CMAKE_ARGS "${CMAKE_ARGS};-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>")
#  INSTALL_COMMAND ""
#  CONFIGURE_COMMAND "")
ExternalProject_Get_property(odgi INSTALL_DIR)
set(odgi_INCLUDE "${INSTALL_DIR}/include")
set(odgi_LIB "${INSTALL_DIR}/lib")
set(dynamic_INCLUDE ${CMAKE_SOURCE_DIR}/deps/odgi/deps/DYNAMIC/include)
set(hopscotch_map_INCLUDE ${CMAKE_SOURCE_DIR}/deps/odgi/deps/hopscotch-map/include)
set(sparsepp_INCLUDE ${CMAKE_SOURCE_DIR}/deps/odgi/deps/sparsepp/sparsepp)
set(flathashmap_INCLUDE ${CMAKE_SOURCE_DIR}/deps/odgi/deps/flat_hash_map)
set(random_dist_INCLUDE ${CMAKE_SOURCE_DIR}/deps/odgi/deps/cpp_random_distributions)

file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/deps/odgi/include)
execute_process(COMMAND bash scripts/generate_git_version.sh include WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/deps/odgi)

# mmmulti (memory mapped multimap)
ExternalProject_Add(mmmulti
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/mmmulti"
  BUILD_COMMAND ""
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(mmmulti SOURCE_DIR)
set(mmmulti_INCLUDE "${SOURCE_DIR}/src")

# In-place Parallel Super Scalar Samplesort (IPS⁴o), header only
ExternalProject_Add(ips4o
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/ips4o"
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(ips4o SOURCE_DIR)
set(ips4o_INCLUDE "${SOURCE_DIR}")

# atomic queue library
ExternalProject_Add(atomicqueue
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/atomic_queue"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(atomicqueue SOURCE_DIR)
set(atomicqueue_INCLUDE "${SOURCE_DIR}")

# atomic bitvector class
ExternalProject_Add(atomicbitvector
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/atomicbitvector"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(atomicbitvector SOURCE_DIR)
set(atomicbitvector_INCLUDE "${SOURCE_DIR}/include")

# taywee's C++ args library, header only
ExternalProject_Add(tayweeargs
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/args"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(tayweeargs SOURCE_DIR)
set(tayweeargs_INCLUDE "${SOURCE_DIR}")

# gfakluge (now header only)
ExternalProject_Add(gfakluge
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/gfakluge"
  CMAKE_ARGS "${CMAKE_ARGS};-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(gfakluge SOURCE_DIR)
set(gfakluge_INCLUDE "${SOURCE_DIR}/src")
set(gfakluge_tinyFA_INCLUDE "${SOURCE_DIR}/src/tinyFA")
ExternalProject_Get_property(gfakluge INSTALL_DIR)
set(gfakluge_LIB "${INSTALL_DIR}/src/gfakluge")

# cgranges, header only
ExternalProject_Add(cgranges
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/cgranges"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(cgranges SOURCE_DIR)
set(cgranges_INCLUDE "${SOURCE_DIR}/cpp")

# paryfor parallel_for
ExternalProject_Add(paryfor
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/paryfor"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(paryfor SOURCE_DIR)
set(paryfor_INCLUDE "${SOURCE_DIR}")

ExternalProject_Add(libbf
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/libbf"
  CMAKE_ARGS "${CMAKE_ARGS};-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>")
ExternalProject_Get_property(libbf INSTALL_DIR)
#message(STATUS "libbf target " ${INSTALL_DIR})
set(libbf_INCLUDE "${INSTALL_DIR}/include")
set(libbf_LIB "${INSTALL_DIR}/lib")

add_subdirectory(deps/spoa EXCLUDE_FROM_ALL)
set(spoa_INCLUDE "${CMAKE_SOURCE_DIR}/deps/spoa/include")

set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE Release)

add_library(smoothxg_objs OBJECT
  src/xg.cpp
  src/chain.cpp
  src/prep.cpp
  src/blocks.cpp
  src/smooth.cpp)

add_dependencies(smoothxg_objs handlegraph)
add_dependencies(smoothxg_objs sdsl-lite)
add_dependencies(smoothxg_objs odgi)
add_dependencies(smoothxg_objs mmmulti)
add_dependencies(smoothxg_objs atomicqueue)
add_dependencies(smoothxg_objs atomicbitvector)
add_dependencies(smoothxg_objs tayweeargs)
add_dependencies(smoothxg_objs gfakluge)
add_dependencies(smoothxg_objs ips4o)
add_dependencies(smoothxg_objs cgranges)
add_dependencies(smoothxg_objs paryfor)
add_dependencies(smoothxg_objs libbf)

set_target_properties(smoothxg_objs PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

#set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

# The following folder will be included
include_directories("${PROJECT_SOURCE_DIR}")

# build a static library
#add_library(libsmoothxg STATIC ${CMAKE_SOURCE_DIR}/src/smoothxg.cpp)
#add_library(libsmoothxg STATIC $<TARGET_OBJECTS:smoothxg_objs>)
#set_target_properties(libsmoothxg PROPERTIES OUTPUT_NAME "smoothxg")

set(smoothxg_HEADERS
  ${CMAKE_SOURCE_DIR}/src/*.hpp)

# set up our target executable and specify its dependencies and includes
#add_executable(smoothxg
#  $<TARGET_OBJECTS:smoothxg_objs>
#  ${CMAKE_SOURCE_DIR}/src/main.cpp
#)

#target_include_directories(smoothxg PUBLIC
set(smoothxg_INCLUDES
  "${CMAKE_SOURCE_DIR}/src"
  "${odgi_INCLUDE}"
  "${handlegraph_INCLUDE}"
  "${sdsl-lite_INCLUDE}"
  "${sdsl-lite-divsufsort_INCLUDE}"
  "${tayweeargs_INCLUDE}"
  "${ips4o_INCLUDE}"
  "${mmmulti_INCLUDE}"
  "${atomicqueue_INCLUDE}"
  "${atomicbitvector_INCLUDE}"
  "${gfakluge_INCLUDE}"
  "${gfakluge_tinyFA_INCLUDE}"
  "${cgranges_INCLUDE}"
  "${spoa_INCLUDE}"
  "${dynamic_INCLUDE}"
  "${hopscotch_map_INCLUDE}"
  "${sparsepp_INCLUDE}"
  "${flathashmap_INCLUDE}"
  "${paryfor_INCLUDE}"
  "${libbf_INCLUDE}"
  "${random_dist_INCLUDE}")

set(smoothxg_LIBS
  "${sdsl-lite_LIB}/libsdsl.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a"
  "${handlegraph_LIB}/libhandlegraph.a"
  "${odgi_LIB}/libodgi.a"
  "-latomic")

target_include_directories(smoothxg_objs PUBLIC ${smoothxg_INCLUDES})

add_executable(smoothxg
  $<TARGET_OBJECTS:smoothxg_objs>
  ${CMAKE_SOURCE_DIR}/src/main.cpp)
target_link_libraries(smoothxg spoa)
target_link_libraries(smoothxg ${smoothxg_LIBS})
set_target_properties(smoothxg PROPERTIES OUTPUT_NAME "smoothxg")
target_include_directories(smoothxg PUBLIC ${smoothxg_INCLUDES})

#add_library(libsmoothxg_static STATIC $<TARGET_OBJECTS:smoothxg_objs>)
#set_target_properties(libsmoothxg_static PROPERTIES OUTPUT_NAME "smoothxg")
#set_target_properties(libsmoothxg_static PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
#add_library(libsmoothxg_shared SHARED $<TARGET_OBJECTS:smoothxg_objs>)
#set_target_properties(libsmoothxg_shared PROPERTIES OUTPUT_NAME "smoothxg")
#set_target_properties(libsmoothxg_shared PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

#set_target_properties(libsmoothxg_static PROPERTIES PUBLIC_HEADER "${smoothxg_HEADERS}")
#set_target_properties(libsmoothxg_shared PROPERTIES PUBLIC_HEADER "${smoothxg_HEADERS}")

install(TARGETS smoothxg DESTINATION bin)
#install(TARGETS libsmoothxg_static ARCHIVE DESTINATION lib LIBRARY DESTINATION lib PUBLIC_HEADER DESTINATION include/smoothxg)
#install(TARGETS libsmoothxg_shared ARCHIVE DESTINATION lib LIBRARY DESTINATION lib PUBLIC_HEADER DESTINATION include/smoothxg)

if (APPLE)
elseif (TRUE)
  if (BUILD_STATIC)
    set(CMAKE_EXE_LINKER_FLAGS "-static")
  endif()
endif()
