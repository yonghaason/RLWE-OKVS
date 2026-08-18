set(SRCS
    "rlwe-okvs/okvs.cpp"
    "rlwe-okvs/sspmt.cpp"
    "rlwe-okvs/layer_budget.cpp"
    "band_okvs/band_okvs.cc"
    "band_okvs/band.cc"
    "band_okvs/oprf.cpp"
    "GMW/Circuit.cpp"
    "GMW/Gmw.cpp"
    "GMW/SilentTripleGen.cpp"
    "PS/benes.cpp"
    "PS/PermuteShare.cpp"
    "pso.cpp"
)

add_library(RLWEOKVS STATIC ${SRCS})
target_link_libraries(RLWEOKVS SEAL::seal oc::libOTe xxhash)
