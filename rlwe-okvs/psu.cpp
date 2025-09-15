#include "psu.h"
#include "rlwe-okvs/rpmt.h"
#include "libOTe/TwoChooseOne/Silent/SilentOtExtReceiver.h"
#include "libOTe/TwoChooseOne/Silent/SilentOtExtSender.h"
#include "libOTe/Vole/Silent/SilentVoleReceiver.h"
#include "libOTe/Vole/Silent/SilentVoleSender.h"
#include "band_okvs/band_okvs.h"
#include "band_okvs/band.h"
#include <memory>

#include <set>

using namespace std;
using namespace seal;
using namespace oc;

namespace rlweOkvs 
{
    Proto PsuSender::run(
        const std::vector<oc::block> &Y,
        Socket &chl)
    {
        vector<block> FY;
        co_await oprf(Y, FY, chl);
        
        RpmtSender rpmtSender;
        uint32_t logp = 60;
        uint32_t numSlots = 1 << 13;
        rpmtSender.init(mN, mNreceiver, logp, numSlots);
        rpmtSender.setTimer(getTimer());

        vector<uint32_t> ot_idx;
        co_await rpmtSender.run(FY, ot_idx, chl);

        auto comm = chl.bytesSent() + chl.bytesReceived();
        
        SilentOtExtSender otSender;
        vector<array<block, 2>> sendMsgs(mN);
        for (size_t i = 0; i < mN; i++) {
            sendMsgs[i][1] = oc::AllOneBlock;
            sendMsgs[i][0] = Y[ot_idx[i]];
        }
        co_await otSender.sendChosen(sendMsgs, mPrng, chl);

        comm = chl.bytesSent() + chl.bytesReceived() - comm;
        cout << "FinalOT takes " << comm << " bytes" << endl;

        setTimePoint("Sender::Final OT");
    };

    Proto PsuSender::oprf(
        const std::vector<oc::block> &Y,
        std::vector<oc::block> &FY,
        Socket &chl)
    {
        auto m = static_cast<size_t>((1.1 * mNreceiver));
        block Delta;

        auto comm = chl.bytesSent() + chl.bytesReceived();

        SilentVoleSender<block, block, oc::CoeffCtxGF128> voleSender;
        co_await voleSender.silentSendInplace(Delta, m, mPrng, chl);

        comm = chl.bytesSent() + chl.bytesReceived() - comm;
        cout << "VOLE takes " << comm << " bytes" << endl;

        vector<block> pp;
        co_await chl.recvResize(pp);

        cout << "Sender receives OPRF ("
             << pp.size() * sizeof(block) << " Bytes)" << endl;


        band_okvs::BandOkvs okvs;        
        auto band_length = 196;
        okvs.Init(mNreceiver, m, band_length, oc::ZeroBlock);
        
        vector<block> kk(m);
        for (size_t i = 0; i < m; i++) {
            kk[i] = voleSender.mB[i] ^ Delta.gf128Mul(pp[i]);
        }

        FY.resize(mN);
        okvs.Decode(Y.data(), kk.data(), FY.data());

        vector<block> h(mN);
        oc::mAesFixedKey.hashBlocks(Y, h);
        
        for (size_t i = 0; i < mN; i++) {
            FY[i] = FY[i] ^ Delta.gf128Mul(h[i]);
        }

        // Random Oracle
        oc::mAesFixedKey.hashBlocks(FY, FY);

        setTimePoint("Sender::OPRF");
    }

    Proto PsuReceiver::run(
        const std::vector<oc::block> &X, 
        std::vector<oc::block>& D, 
        Socket &chl)
    {
        vector<block> FX;
        co_await oprf(X, FX, chl);        
                
        RpmtReceiver rpmtReceiver;
        uint32_t logp = 60;
        uint32_t numSlots = 1 << 13;
        rpmtReceiver.init(mN, mNsender, logp, numSlots);
        rpmtReceiver.setTimer(getTimer());

        BitVector rpmt;
        co_await rpmtReceiver.run(FX, rpmt, chl);

        SilentOtExtReceiver otReceiver;
        vector<block> recvMsgs(mNsender);
        co_await otReceiver.receiveChosen(rpmt, recvMsgs, mPrng, chl);

        for (size_t i = 0; i < recvMsgs.size(); i++) {
            if (recvMsgs[i] != oc::AllOneBlock) {
                D.push_back(recvMsgs[i]);
            }
        }

        setTimePoint("Receiver::Final OT");
    };

    Proto PsuReceiver::oprf(
        const std::vector<oc::block> &X,
        std::vector<oc::block> &FX,
        Socket &chl)
    {
        band_okvs::BandOkvs okvs;        
        auto m = static_cast<size_t>((1.1 * mN));
        auto band_length = 196;
        okvs.Init(mN, m, band_length, oc::ZeroBlock);

        vector<block> pp(m);
        vector<block> val(mN);
        oc::mAesFixedKey.hashBlocks(X, val);

        okvs.Encode(X.data(), val.data(), pp.data());        

        SilentVoleReceiver<block, block, oc::CoeffCtxGF128> voleReceiver;
        co_await voleReceiver.silentReceiveInplace(m, mPrng, chl);

        for (size_t i = 0; i < m; i++) {
            pp[i] = pp[i] ^ voleReceiver.mC[i];
        }

        co_await chl.send(move(pp));

        FX.resize(mN);
        okvs.Decode(X.data(), voleReceiver.mA.data(), FX.data());   

        // Random Oracle
        oc::mAesFixedKey.hashBlocks(FX, FX);

        setTimePoint("Receiver::OPRF");
    }
}