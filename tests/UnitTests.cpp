#include "UnitTests.h"
#include "OKVS_tests.h"
#include "OPRF_tests.h"
#include "SEAL_tests.h"
#include "PSO_tests.h"
#include "GMW_tests.h"
#include "SEQ_tests.h"
#include "SSPMT_tests.h"
#include "PS_tests.h"

#include <functional>

namespace rlweOkvsTests {
    oc::TestCollection Tests([](oc::TestCollection& t) {
    t.add("OKVS_width_test                  ", width_test),
    // t.add("OKVS_encode_test                 ", encode_test),
    // t.add("OKVS_decode_test                 ", decode_test),
    // t.add("OPRF_protocol_test               ", oprf_protocol_test);
    t.add("Sequencing_test                  ", opti_sequencing_test);
    t.add("SSPMT_test                       ", sspmt_test);
    t.add("PSI_card_test                    ", psi_card_test);
    t.add("PSI_sum_test                     ", psi_sum_test);
    t.add("PSI_threshold_test               ", psi_threshold_test);
    t.add("PSU_test                         ", psu_test);
    t.add("SSPMT_net_test                   ", sspmt_net_test);
    // t.add("PermuteShare_test                ", permute_share_test);
    // t.add("GMW_iszero_test                  ", Gmw_iszero_test);
    // t.add("GMW_threshold_test               ", Gmw_threshold_test);
    });
}
