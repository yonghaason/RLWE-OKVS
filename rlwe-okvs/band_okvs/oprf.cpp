#include "oprf.h"
#include "libOTe/TwoChooseOne/Silent/SilentOtExtReceiver.h"
#include "libOTe/TwoChooseOne/Silent/SilentOtExtSender.h"
#include "libOTe/Vole/Silent/SilentVoleReceiver.h"
#include "libOTe/Vole/Silent/SilentVoleSender.h"
#include "band_okvs.h"
#include "band.h"
#include <memory>

#include <set>

using namespace std;
using namespace seal;
using namespace oc;

namespace rlweOkvs 
{
    Proto OprfSender::run(
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

    Proto OprfReceiver::run(
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

    }
}