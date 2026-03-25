cmake_minimum_required (VERSION 3.15)

project(RLWEOKVS)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-O3 -march=native -Wall -Wfatal-errors -flto -mavx512f -mavx512dq -mavx512vl")


find_package(SEAL 4.1 REQUIRED HINTS "${CMAKE_SOURCE_DIR}/thirdparty/install/SEAL")
find_package(libOTe REQUIRED HINTS "${CMAKE_SOURCE_DIR}/thirdparty/install/libOTe")
set(XXHASH_ROOT "${CMAKE_SOURCE_DIR}/thirdparty/xxHash")
add_library(xxhash STATIC IMPORTED GLOBAL)
set_target_properties(xxhash PROPERTIES
    IMPORTED_LOCATION "${XXHASH_ROOT}/libxxhash.a"
    INTERFACE_INCLUDE_DIRECTORIES "${XXHASH_ROOT}"
)

add_subdirectory(rlwe-okvs)
add_subdirectory(tests)

add_executable(run exec/main.cpp)

target_link_libraries(run RLWEOKVS RLWEOKVS_Tests)
target_include_directories(run PRIVATE
    ${CMAKE_SOURCE_DIR}
)
