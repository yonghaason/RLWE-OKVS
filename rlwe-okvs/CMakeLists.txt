set(SRCS
    "rlwe-okvs/okvs.cpp"
    "rlwe-okvs/rpmt.cpp"
    "band_okvs/band_okvs.cc"
    "band_okvs/band.cc"
    "psu.cpp"
)

add_library(RLWEOKVS STATIC ${SRCS})
target_link_libraries(RLWEOKVS SEAL::seal oc::libOTe)