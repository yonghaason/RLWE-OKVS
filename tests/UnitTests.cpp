#include "UnitTests.h"
#include "OKVS_Tests.h"
#include "SEAL_tests.h"
#include "RPMT_tests.h"

#include <functional>

namespace rlweOkvsTests {
    oc::TestCollection Tests([](oc::TestCollection& t) {
    t.add("OKVS_encode_test                 ", encode_test),
    t.add("OKVS_decode_test                 ", decode_test),
    t.add("context_consistency_test         ", context_consistency_test),
    t.add("RPMT_receiver_test                 ", receiver_test),
    t.add("RPMT_sender_test                 ", sender_test);
    });
}