set(SRCS
    "rlwe-okvs/okvs.cpp"
    "rlwe-okvs/sspmt.cpp"
    "band_okvs/band_okvs.cc"
    "band_okvs/band.cc"
    "GMW/Circuit.cpp"
    "GMW/Gmw.cpp"
    "GMW/SilentTripleGen.cpp"
    "GMW/TripleGen.cpp"
    "psu.cpp"
)

add_library(RLWEOKVS STATIC ${SRCS})
target_link_libraries(RLWEOKVS SEAL::seal oc::libOTe)