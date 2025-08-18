#include "UnitTests.h"
#include "OKVS_Tests.h"

#include <functional>

namespace rlweOkvsTests {
    oc::TestCollection Tests([](oc::TestCollection& t) {
    t.add("OKVS_encode_test                 ", encode_test),
    t.add("OKVS_decode_test                 ", decode_test);
    });
}