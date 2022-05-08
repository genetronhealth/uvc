cmake_minimum_required(VERSION 3.13)
set(PROJECT_NAME uvc)
project(${PROJECT_NAME} C CXX)

set(CMAKE_CXX_STANDARD 14)
add_definitions(-DCOMMIT_VERSION=":latest" -DCOMMIT_DIFF_SH="Schaudge King" -DCOMMIT_DIFF_FULL="")
#add_definitions(-DUSE_STDLIB_THREAD)
add_compile_options(-std=c++14 -static-libstdc++ -Wall)

### add all required header files
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/CLI11-1.7.1)

### discover samtools and htslib
find_library(SAMTOOLS_BAM_LIB NAMES bam HINTS ${SAMTOOLS_DIR} STATIC IMPORTED)
find_library(HTS_LIB NAMES hts HINTS ${HTSLIB_DIR} STATIC IMPORTED)

### add all source file to the SRC_LIST *variable*
# omit

add_executable(${PROJECT_NAME} main.cpp CmdLineArgs.cpp common.cpp grouping.cpp logging.cpp version.cpp)
target_link_libraries(${PROJECT_NAME} -fopenmp hts pthread m z bz2 lzma)

add_executable(debarcode debarcode_main.c)
target_link_libraries(debarcode hts pthread m z bz2 lzma)
