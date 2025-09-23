#pragma once
#include "cryptoTools/Common/Defines.h"
#include "cryptoTools/Common/block.h"
#include "cryptoTools/Common/BitVector.h"
#include "cryptoTools/Crypto/PRNG.h"
#include "cryptoTools/Common/Timer.h"
#include "cryptoTools/Network/Channel.h"
#include "coproto/coproto.h"

#include "seal/seal.h"

namespace rlweOkvs
{
    using Proto = coproto::task<>;
    using Socket = coproto::Socket;

    class OprfSender: public oc::TimerAdapter
    {
        uint32_t mN;
        uint32_t mNreceiver;
        oc::PRNG mPrng;

    public:
        void init(uint32_t n, uint32_t nReceiver, oc::block seed) {
            mN = n;
            mNreceiver = nReceiver;
            mPrng.SetSeed(seed);
        };

        Proto run(
            const std::vector<oc::block> &Y,
            std::vector<oc::block> &FY,
            Socket &chl);
    };

    class OprfReceiver: public oc::TimerAdapter
    {
        uint32_t mN;
        uint32_t mNsender;
        oc::PRNG mPrng;

    public:
        void init(uint32_t n, uint32_t nSender, oc::block seed) {
            mN = n;
            mNsender = nSender;
            mPrng.SetSeed(seed);
        };

        Proto run(
            const std::vector<oc::block> &X,
            std::vector<oc::block> &FX,
            Socket &chl);
    };
};
