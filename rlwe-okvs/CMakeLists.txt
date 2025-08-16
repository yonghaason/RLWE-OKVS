set(SRCS
    "rlwe-okvs/okvs.cpp"
)

add_library(RLWEOKVS STATIC ${SRCS})
target_link_libraries(RLWEOKVS SEAL::seal oc::libOTe)