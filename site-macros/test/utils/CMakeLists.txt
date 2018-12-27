# ==============================================================================
#   PROJECT TEST utils.cmake 
# ==============================================================================
project( test_utils C CXX Fortran )
cmake_minimum_required(VERSION 2.6)

# ==============================================================================
#   LOAD
# ==============================================================================
# TODO: currently manual test. convert to automatic somehow
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../Modules)
include(utils)

list(APPEND x  "/program files/am/cik")
list(APPEND x  "/program files/am/micro soft/cik")

list_stringize(y ${x} "a b")
file(WRITE tester.txt " 
 set(x ${y})
 ")
