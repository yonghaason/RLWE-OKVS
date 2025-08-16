cmake_minimum_required (VERSION 3.15)

project(RLWEOKVSPSU)

set(CMAKE_CXX_STANDARD 20)

find_package(SEAL 4.1 REQUIRED HINTS "${CMAKE_SOURCE_DIR}/thirdparty/install/SEAL")
find_package(libOTe REQUIRED HINTS "${CMAKE_SOURCE_DIR}/thirdparty/install/libOTe")

add_subdirectory(rlwe-okvs)
add_subdirectory(tests)
add_subdirectory(exec)