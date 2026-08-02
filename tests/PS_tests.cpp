#include "PS_tests.h"
#include "PS/PermuteShare.h"

#include "cryptoTools/Common/Defines.h"
#include "cryptoTools/Common/CLP.h"
#include "cryptoTools/Common/Timer.h"
#include "cryptoTools/Crypto/PRNG.h"
#include "coproto/Socket/LocalAsyncSock.h"

#include "macoro/sync_wait.h"
#include "macoro/when_all.h"
#include "macoro/thread_pool.h"

#include <iostream>
#include <numeric>
#include <vector>

using namespace std;
using namespace oc;
using namespace rlweOkvs;

// F_PS: the sender holds pi, the receiver holds x; they end with a, b such
// that a ^ b = pi(x). Checked against the permutation applied in the clear.
// The switch OTs are generated in a separate phase first, which is where they
// belong -- they are random and depend only on n.
void permute_share_test(const oc::CLP& cmd)
{
    u64 n = cmd.getOr("n", 1ull << cmd.getOr("nn", 12));
    u64 nt = cmd.getOr("nt", 1);

    PRNG prng(block(3128, 8712));

    // A random permutation for the sender, and a random input for the receiver.
    std::vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    for (u64 i = n - 1; i > 0; --i) std::swap(perm[i], perm[prng.get<u64>() % (i + 1)]);

    BitVector input(n);
    input.randomize(prng);

    macoro::thread_pool pool0, pool1;
    auto e0 = pool0.make_work();
    pool0.create_threads(nt);
    auto e1 = pool1.make_work();
    pool1.create_threads(nt);
    auto socket = coproto::LocalAsyncSocket::makePair();
    socket[0].setExecutor(pool0);
    socket[1].setExecutor(pool1);

    Timer timer_s, timer_r;
    PermuteShareSender sender;
    PermuteShareReceiver receiver;
    sender.setTimer(timer_s);
    receiver.setTimer(timer_r);
    sender.init(n);
    receiver.init(n);

    SilentOtExtReceiver otRecv;
    SilentOtExtSender otSend;
    PRNG prngS(block(1, 2)), prngR(block(3, 4));

    timer_s.setTimePoint("start");
    timer_r.setTimePoint("start");

    auto runBoth = [&](auto&& p0, auto&& p1) {
        auto r = macoro::sync_wait(macoro::when_all_ready(
            std::move(p0) | macoro::start_on(pool0),
            std::move(p1) | macoro::start_on(pool1)));
        std::get<0>(r).result();
        std::get<1>(r).result();
    };

    runBoth(sender.setup(otRecv, prngS, socket[0]),
            receiver.setup(otSend, prngR, socket[1]));
    const u64 setupBytes = socket[0].bytesSent() + socket[1].bytesSent();

    sender.setPermutation(perm);

    BitVector senderShare, receiverShare;
    runBoth(sender.run(senderShare, socket[0]),
            receiver.run(input, receiverShare, prngR, socket[1]));

    if (senderShare.size() != n || receiverShare.size() != n) throw RTE_LOC;

    for (u64 i = 0; i < n; ++i) {
        const u8 got = senderShare[i] ^ receiverShare[i];
        // benes routes input position perm[i] to output i
        const u8 want = input[perm[i]];
        if (got != want) {
            cout << "mismatch at " << i << ": got " << (int)got << ", want "
                 << (int)want << endl;
            throw RTE_LOC;
        }
    }

    if (cmd.isSet("v")) {
        const double MB = 1024 * 1024;
        const u64 total = socket[0].bytesSent() + socket[1].bytesSent();
        cout << "n = " << n << ", switches = " << permuteShareSwitchCount(n)
             << endl;
        cout << timer_s << endl << timer_r << endl;
        cout << "comm: offline " << setupBytes / MB << " MB + online "
             << (total - setupBytes) / MB << " MB" << endl;
    }
}
