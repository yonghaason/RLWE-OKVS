#pragma once
#include "cryptoTools/Common/Defines.h"
#include "cryptoTools/Common/block.h"
#include "cryptoTools/Common/BitVector.h"
#include "cryptoTools/Crypto/PRNG.h"
#include "cryptoTools/Common/Timer.h"
#include "cryptoTools/Network/Channel.h"
#include "coproto/coproto.h"
#include "rlwe-okvs/sspmt.h"

#include "seal/seal.h"

int psu(int n, int t_s, int m,  vector<int> vector);


namespace rlweOkvs
{
    using Proto = coproto::task<>;
    using Socket = coproto::Socket;

    class PsuSender: public oc::TimerAdapter
    {
        uint32_t mN;
        uint32_t mNreceiver;
        oc::PRNG mPrng;        
        sspmtParams mSsParams;

    public:
        void init(uint32_t n, uint32_t nReceiver, 
            sspmtParams ssParams, oc::block seed) {
            mN = n;
            mNreceiver = nReceiver;
            mSsParams = ssParams;
            mPrng.SetSeed(seed);
        };

        Proto run(
            const std::vector<oc::block> &Y, 
            Socket &chl);
    };

    class PsuReceiver: public oc::TimerAdapter
    {
        uint32_t mN;
        uint32_t mNsender;
        oc::PRNG mPrng;
        sspmtParams mSsParams;

    public:
        void init(uint32_t n, uint32_t nSender, 
            sspmtParams ssParams, oc::block seed) {
            mN = n;
            mNsender = nSender;
            mSsParams = ssParams;
            mPrng.SetSeed(seed);
        };

        Proto run(
            const std::vector<oc::block> &X, 
            std::vector<oc::block>& D, 
            Socket &chl);
    };
};
