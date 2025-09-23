#include "UnitTests.h"
#include "OKVS_Tests.h"
#include "OPRF_Tests.h"
#include "SEAL_tests.h"
#include "RPMT_tests.h"
#include "PSU_tests.h"
#include "GMW_tests.h"

#include <functional>

namespace rlweOkvsTests {
    oc::TestCollection Tests([](oc::TestCollection& t) {
    t.add("OKVS_encode_test                 ", encode_test),
    t.add("OKVS_decode_test                 ", decode_test),
    t.add("context_consistency_test         ", context_consistency_test),
    t.add("RPMT_protocol_test               ", rpmt_protocol_test);
    t.add("OPRF_protocol_test               ", oprf_protocol_test);
    t.add("PSU_protocol_test                ", psu_protocol_test);
    t.add("GMW_iszero_test                  ", Gmw_iszero_test);
    });
}