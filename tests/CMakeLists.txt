set(SRCS
	OKVS_tests.cpp
    OPRF_Tests.cpp
    SEAL_tests.cpp
    RPMT_tests.cpp
    PSU_tests.cpp
    GMW_tests.cpp
    UnitTests.cpp
)

add_library(RLWEOKVS_Tests ${SRCS})
target_link_libraries(RLWEOKVS_Tests RLWEOKVS)
target_include_directories(RLWEOKVS_Tests PRIVATE
    ${CMAKE_SOURCE_DIR}/rlwe-okvs
)