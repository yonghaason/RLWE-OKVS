#include "psu.h"

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
        rpmtSender.init(mN, mNreceiver, mPmtParams, mPrng.get());
        rpmtSender.setTimer(getTimer());

        co_await rpmtSender.run(Y, chl);

        auto ot_idx = rpmtSender.get_ot_idx();
        
        auto comm = chl.bytesSent() + chl.bytesReceived();
                
        vector<array<block, 2>> rotMsgs(mN);
        co_await otSender.send(rotMsgs, mPrng, chl);

        vector<block> otp(mN);
        
        for (size_t i = 0; i < mN; i++) {
            otp[i] = rotMsgs[i][0] ^ Y[ot_idx[i]];
        }

        co_await chl.send(std::move(otp));

        comm = chl.bytesSent() + chl.bytesReceived() - comm;
        cout << "FinalOT takes " << comm << " bytes" << endl;

        setTimePoint("Sender::Final OT");
    };

    Proto PsuReceiver::run(
        const std::vector<oc::block> &X, 
        std::vector<oc::block>& D, 
        Socket &chl)
    {
        rpmtReceiver.init(mN, mNsender, mPmtParams, mPrng.get());
        rpmtReceiver.setTimer(getTimer());

        BitVector rpmt;
        co_await rpmtReceiver.run(X, rpmt, chl);
        
        vector<block> rotMsgs(mNsender);
        co_await otReceiver.receive(rpmt, rotMsgs, mPrng, chl);

        vector<block> otp;
        co_await chl.recvResize(otp);

        for (size_t i = 0; i < otp.size(); i++) {
            auto tmp = otp[i] ^ rotMsgs[i];
            if (tmp.mData[1] == 0) {
                D.push_back(tmp);
            }
        }

        setTimePoint("Receiver::Final OT");
    };
}
