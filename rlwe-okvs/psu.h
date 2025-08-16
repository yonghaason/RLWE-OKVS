#include "cryptoTools/Common/Timer.h"
#include "cryptoTools/Common/BitVector.h"
#include "band_okvs.h"
#include "uint.h"

namespace secJoin
{
    using Proto = coproto::task<>;
    using Socket = coproto::Socket;

    class AltModPsuSender : public oc::TimerAdapter
    {        
        oc::Timer timer;

    public:

        Proto run(span<block> Y, PRNG &prng, Socket &chl);
    };

    class AltModPsuReceiver :public oc::TimerAdapter
    {
        oc::Timer timer;
        
    public:

        Proto run(span<block> X, std::vector<block>& D, PRNG &prng, Socket &chl);

    };

}
