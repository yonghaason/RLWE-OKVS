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
    // Benchmark entry points (see benchmark.sh): indices 1-5, one party per
    // process over a real TCP connection.
    t.add("OKVS_width_test                  ", width_test);              // 0
    t.add("SSPMT_net_test                   ", sspmt_net_test);          // 1
    t.add("PSU_net_test                     ", psu_net_test);            // 2
    t.add("PSI_card_net_test                ", psi_card_net_test);       // 3
    t.add("PSI_sum_net_test                 ", psi_sum_net_test);        // 4
    t.add("PSI_threshold_net_test           ", psi_threshold_net_test);  // 5
    // In-process variants (both parties in one process; kept for debugging).
    t.add("Sequencing_test                  ", opti_sequencing_test);
    t.add("SSPMT_test                       ", sspmt_test);
    t.add("PSI_card_test                    ", psi_card_test);
    t.add("PSI_sum_test                     ", psi_sum_test);
    t.add("PSI_threshold_test               ", psi_threshold_test);
    t.add("PSU_test                         ", psu_test);
    // t.add("OKVS_encode_test                 ", encode_test);
    // t.add("OKVS_decode_test                 ", decode_test);
    // t.add("OPRF_protocol_test               ", oprf_protocol_test);
    // t.add("PermuteShare_test                ", permute_share_test);
    // t.add("GMW_iszero_test                  ", Gmw_iszero_test);
    // t.add("GMW_threshold_test               ", Gmw_threshold_test);
    });
}
