set(SRCS
    "rlwe-okvs/okvs.cpp"
    "rlwe-okvs/rpmt.cpp"
)

add_library(RLWEOKVS STATIC ${SRCS})
target_link_libraries(RLWEOKVS SEAL::seal oc::libOTe)