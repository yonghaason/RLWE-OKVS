#include "pso.h"

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

    Proto PsiCardSumSender::run(
        const std::vector<oc::block> &Y,
        const std::vector<oc::u32> &payloads,
        Socket &chl)
    {
        if (payloads.size() != Y.size()) {
            throw RTE_LOC;
        }

        rpmtSender.init(mN, mNreceiver, mPmtParams, mPrng.get());
        rpmtSender.setTimer(getTimer());
        co_await rpmtSender.run(Y, chl);

        auto ot_idx = rpmtSender.get_ot_idx();
        auto comm = chl.bytesSent() + chl.bytesReceived();

        vector<array<block, 2>> rotMsgs(mN);
        co_await otSender.send(rotMsgs, mPrng, chl);

        vector<u32> cotsend(mN);
        oc::u32 maskSum = 0;
        for (size_t i = 0; i < mN; i++) {
            maskSum += rotMsgs[i][0].mData[0];
            cotsend[i] = rotMsgs[i][0].mData[0] + rotMsgs[i][1].mData[0] + payloads[ot_idx[i]];
        }

        co_await chl.send(std::move(cotsend));
        co_await chl.send(maskSum);

        comm = chl.bytesSent() + chl.bytesReceived() - comm;
        cout << "FinalOT(card-sum) takes " << comm << " bytes" << endl;
        setTimePoint("Sender::Final OT Card Sum");
    };

    Proto PsiCardSumReceiver::run(
        const std::vector<oc::block> &X,
        oc::u64& psiCardSum,
        Socket &chl)
    {
        rpmtReceiver.init(mN, mNsender, mPmtParams, mPrng.get());
        rpmtReceiver.setTimer(getTimer());

        BitVector rpmt;
        co_await rpmtReceiver.run(X, rpmt, chl);

        vector<block> rotMsgs(mNsender);
        co_await otReceiver.receive(rpmt, rotMsgs, mPrng, chl);

        vector<u32> cotrecv(mN);
        co_await chl.recvResize(cotrecv);

        oc::u32 maskSum;
        co_await chl.recv(maskSum);

        oc::u32 maskedSum = 0;
        for (size_t i = 0; i < rotMsgs.size(); i++) {
            u32 recovered = rpmt[i] ? (cotrecv[i] - rotMsgs[i].mData[0]) : rotMsgs[i].mData[0];
            maskedSum += recovered;
        }

        psiCardSum = maskedSum - maskSum;
        setTimePoint("Receiver::Final OT Card Sum");
    };
}
