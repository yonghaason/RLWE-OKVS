set(SRCS
	OKVS_tests.cpp
    UnitTests.cpp
)

add_library(RLWEOKVS_Tests ${SRCS})
target_link_libraries(RLWEOKVS_Tests RLWEOKVS)
target_include_directories(RLWEOKVS_Tests PRIVATE
    ${CMAKE_SOURCE_DIR}/rlwe-okvs
)