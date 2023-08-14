# Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

# Project's name
project(smoothxg)

# We build using c++17
set(CMAKE_CXX_STANDARD 17)

enable_testing()

include(CheckIPOSupported) # adds lto
check_ipo_supported(RESULT ipo_supported OUTPUT output)
SET(CMAKE_RANLIB "gcc-ranlib") # too try lto with older runtimes

# This builds a static version of ./bin/smoothxg
option(BUILD_STATIC "Build static binary" OFF)
if (BUILD_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(BUILD_SHARED_LIBS OFF)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

find_package(OpenMP)
# if(OpenMP_CXX_FOUND)
#     target_link_libraries(MyTarget PUBLIC OpenMP::OpenMP_CXX)
# endif()

find_package(PkgConfig REQUIRED)
find_package(ZLIB REQUIRED)
# find_package(ZSTD REQUIRED) - needs recent versions of cmake

# Preload the following libraries before running tests
set(PRELOAD "libasan.so:libjemalloc.so.2")

if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Choose the type of build, options are: Release Debug Generic!")
endif()

# set(CMAKE_BUILD_TYPE Debug) -- don't uncomment this, instead run
#   cmake -DCMAKE_BUILD_TYPE=Debug ..

if (${CMAKE_BUILD_TYPE} MATCHES Release)
  if(NOT EXTRA_FLAGS)
    set(EXTRA_FLAGS "-Ofast -march=native -flto -fno-fat-lto-objects")
  endif()
  set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG") # reset CXX_FLAGS to replace -O3 with -Ofast

  # Increase SPOA's performance
  set(spoa_optimize_for_portability OFF)
  set(spoa_optimize_for_native ON)
else()
  # Increase SPOA's portability (for Docker/Singularity)
  set(spoa_optimize_for_portability ON)
  set(spoa_optimize_for_native OFF)
endif ()

if (${CMAKE_BUILD_TYPE} MATCHES Debug)
  # Debug use the defaults
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O -g -fsanitize=address")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O -g -fsanitize=address")
else()
  # Use all standard-compliant optimizations - always add these:
  set (CMAKE_C_FLAGS "${OpenMP_C_FLAGS} ${PIC_FLAG} ${EXTRA_FLAGS}")
  set (CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} ${PIC_FLAG} ${EXTRA_FLAGS}")
endif ()

if(POA_DEBUG)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPOA_DEBUG")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPOA_DEBUG")
endif ()

#set(THREADS_PREFER_PTHREAD_FLAG ON)
#find_package(Threads REQUIRED)

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

message(STATUS "SMOOTHXG CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
message(STATUS "SMOOTHXG CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
message(STATUS "SMOOTHXG CMAKE_C_FLAGS: ${CMAKE_C_FLAGS}")

# Set the output folder
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
  CMAKE_ARGS "${CMAKE_ARGS};-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE};-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>")
#  INSTALL_COMMAND ""
#  CONFIGURE_COMMAND "")
ExternalProject_Get_property(odgi INSTALL_DIR)
set(odgi_INCLUDE "${INSTALL_DIR}/include")
set(odgi_LIB "${INSTALL_DIR}/lib")
set(dynamic_INCLUDE ${CMAKE_SOURCE_DIR}/deps/odgi/deps/DYNAMIC/include)
set(bbhash_INCLUDE ${CMAKE_SOURCE_DIR}/deps/odgi/deps/BBHash)
set(hopscotch_map_INCLUDE ${CMAKE_SOURCE_DIR}/deps/odgi/deps/hopscotch-map/include)
set(sparsepp_INCLUDE ${CMAKE_SOURCE_DIR}/deps/odgi/deps/sparsepp/sparsepp)
set(flathashmap_INCLUDE ${CMAKE_SOURCE_DIR}/deps/odgi/deps/flat_hash_map)
set(random_dist_INCLUDE ${CMAKE_SOURCE_DIR}/deps/odgi/deps/cpp_random_distributions)

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
set(atomicqueue_INCLUDE "${SOURCE_DIR}/include/atomic_queue")

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

ExternalProject_Add(libbf
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/libbf"
  CMAKE_ARGS "${CMAKE_ARGS};-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>")
ExternalProject_Get_property(libbf INSTALL_DIR)
#message(STATUS "libbf target " ${INSTALL_DIR})
set(libbf_INCLUDE "${INSTALL_DIR}/include")
set(libbf_LIB "${INSTALL_DIR}/lib")

# dirtyzipf pow-approximate Zipf distribution
ExternalProject_Add(dirtyzipf
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/odgi/deps/dirtyzipf"
        UPDATE_COMMAND ""
        INSTALL_COMMAND ""
        BUILD_COMMAND ""
        CONFIGURE_COMMAND "")
ExternalProject_Get_property(dirtyzipf SOURCE_DIR)
set(dirtyzipf_INCLUDE "${SOURCE_DIR}")

# header-only pseudorandom number generator library
ExternalProject_Add(xoshiro
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/odgi/deps/Xoshiro-cpp"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(xoshiro SOURCE_DIR)
set(xoshiro_INCLUDE "${SOURCE_DIR}")

ExternalProject_Add(sautocorr
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/sautocorr"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(sautocorr SOURCE_DIR)
set(sautocorr_INCLUDE "${SOURCE_DIR}")

add_subdirectory(${CMAKE_SOURCE_DIR}/deps/spoa EXCLUDE_FROM_ALL)
set(spoa_INCLUDE "${CMAKE_SOURCE_DIR}/deps/spoa/include")

# deps/abPOA sources are built inline because the upstream CMakeLists.txt file
# is hard coding build parameters and fails on finding zlib.h, for example.
#
# We may realign with upstream and convince them to improve their build system

set(abpoa_SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/abPOA")
add_library(abpoa STATIC
    ${abpoa_SOURCE_DIR}/src/abpoa_align.c
    ${abpoa_SOURCE_DIR}/src/abpoa_graph.c
    ${abpoa_SOURCE_DIR}/src/abpoa_plot.c
    ${abpoa_SOURCE_DIR}/src/abpoa_seed.c
    ${abpoa_SOURCE_DIR}/src/abpoa_seq.c
    ${abpoa_SOURCE_DIR}/src/kalloc.c
    ${abpoa_SOURCE_DIR}/src/kstring.c
    ${abpoa_SOURCE_DIR}/src/simd_abpoa_align.c
    ${abpoa_SOURCE_DIR}/src/simd_check.c
    ${abpoa_SOURCE_DIR}/src/utils.c)

#add_subdirectory(deps/mmmulti/deps/mio)
ExternalProject_Add(mio
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/mmmulti/deps/mio"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(mio SOURCE_DIR)
set(mio_INCLUDE "${SOURCE_DIR}/include")

# mkmh
ExternalProject_Add(mkmh
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/mkmh"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(mkmh SOURCE_DIR)
set(mkmh_INCLUDE "${SOURCE_DIR}")

ExternalProject_Add(edlib
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/edlib"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  TEST_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(edlib SOURCE_DIR)
set(edlib_INCLUDE "${SOURCE_DIR}/edlib/include")

ExternalProject_Add(xxHash
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/xxHash"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  TEST_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(xxHash SOURCE_DIR)
set(xxHash_INCLUDE "${SOURCE_DIR}")

add_subdirectory(deps/WFA EXCLUDE_FROM_ALL)

add_library(smoothxg_objs OBJECT
  src/tempfile.cpp
  deps/xxHash/xxhash.c
  src/xg.cpp
  src/chain.cpp
  src/prep.cpp
  src/cleanup.cpp
  src/blocks.cpp
  ${mkmh_INCLUDE}/murmur3/murmur3.cpp
  src/breaks.cpp
  src/smooth.cpp
  src/utils.cpp
  src/zstdutil.cpp
  ${sautocorr_INCLUDE}/sautocorr.cpp
  src/consensus_graph.cpp)

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
add_dependencies(smoothxg_objs libbf)
add_dependencies(smoothxg_objs dirtyzipf)
add_dependencies(smoothxg_objs xoshiro)
add_dependencies(smoothxg_objs sautocorr)
add_dependencies(smoothxg_objs abpoa)
add_dependencies(smoothxg_objs mio)
add_dependencies(smoothxg_objs mkmh)

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
  "${mio_INCLUDE}"
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
  "${libbf_INCLUDE}"
  "${random_dist_INCLUDE}"
  "${dirtyzipf_INCLUDE}"
  "${xoshiro_INCLUDE}"
  "${sautocorr_INCLUDE}"
  "${abPOA_INCLUDE}"
  "${abPOA_src_INCLUDE}"
  "${bbhash_INCLUDE}"
  "${mkmh_INCLUDE}"
  "${edlib_INCLUDE}"
  "deps/WFA"
  "deps/patchmap")

set(smoothxg_LIBS
  "${sdsl-lite_LIB}/libsdsl.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a"
  "${odgi_LIB}/libodgi.a"
  "${handlegraph_LIB}/libhandlegraph.a"
  "-latomic"
  wfa
  )

target_include_directories(smoothxg_objs PUBLIC ${smoothxg_INCLUDES})

add_executable(smoothxg
  $<TARGET_OBJECTS:smoothxg_objs>
  ${CMAKE_SOURCE_DIR}/src/main.cpp)
target_link_libraries(smoothxg spoa abpoa)
target_link_libraries(smoothxg ${smoothxg_LIBS})

target_link_libraries(smoothxg z zstd jemalloc)

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

# version stuff
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/include)
execute_process(COMMAND bash ${CMAKE_SOURCE_DIR}/scripts/generate_git_version.sh ${CMAKE_SOURCE_DIR}/include)

add_test(
  NAME smoothxg-test
  # 2 cores because of CI limitations
  COMMAND bin/smoothxg -t 2 -g test/data/DRB1-3123.fa.gz.pggb-s3000-p70-n10-a70-K16-k8-w10000-j5000-e5000.seqwish.gfa -j 5k -e 5k -l 700,900,1100 -m test/data/DRB1-3123.fa.gz.pggb-s3000-p70-n10-a70-K16-k8-w10000-j5000-e5000.smooth.maf -C "consensus,10,100:test/data/gi_568815592_32578768-32589835.txt:y,1000:test/data/gi_568815592_32578768-32589835.txt:n,10000" -o test/data/DRB1-3123.fa.gz.pggb-s3000-p70-n10-a70-K16-k8-w10000-j5000-e5000.smooth.gfa -r 12
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set_tests_properties(smoothxg-test PROPERTIES ENVIRONMENT "ASAN_OPTIONS=detect_leaks=1:symbolize=1;LSAN_OPTIONS=verbosity=0:log_threads=1")

if (APPLE)
elseif (TRUE)
  if (BUILD_STATIC)
    set(CMAKE_EXE_LINKER_FLAGS "-static")
  endif()
endif()
