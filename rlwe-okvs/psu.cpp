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
        vector<block> FY;
        
        oprfSender.init(mN, mNreceiver, mPrng.get());
        co_await oprfSender.run(Y, FY, chl);                    
        
        sspmtSender.init(mN, mNreceiver, mSsParams, mPrng.get());
        sspmtSender.setTimer(getTimer());

        BitVector sspmt;
        co_await sspmtSender.run(FY, sspmt, chl);
        auto ot_idx = sspmtSender.get_ot_idx();
        
        auto comm = chl.bytesSent() + chl.bytesReceived();
                
        vector<array<block, 2>> rotMsgs(mN);
        co_await otSender.send(rotMsgs, mPrng, chl);

        vector<block> otp(mN);
        
        for (size_t i = 0; i < mN; i++) {
            otp[i] = rotMsgs[i][sspmt[i]] ^ Y[ot_idx[i]];
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
        vector<block> FX;
        
        oprfreceiver.init(mN, mNsender, mPrng.get());
        co_await oprfreceiver.run(X, FX, chl);        
           
        sspmtReceiver.init(mN, mNsender, mSsParams, mPrng.get());
        sspmtReceiver.setTimer(getTimer());

        BitVector sspmt;
        co_await sspmtReceiver.run(FX, sspmt, chl);

        
        vector<block> rotMsgs(mNsender);
        co_await otReceiver.receive(sspmt, rotMsgs, mPrng, chl);

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