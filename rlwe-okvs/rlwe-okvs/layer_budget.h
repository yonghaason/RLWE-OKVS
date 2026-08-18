#pragma once

#include <cstdint>
#include <vector>

namespace rlweOkvs {

std::vector<uint64_t> hallCaps(uint64_t n, uint32_t numSlots,
                               uint64_t positionRange, uint32_t b,
                               uint32_t lambda);

uint64_t budgetForSpan(const std::vector<uint64_t> &caps, uint32_t b,
                       uint32_t W);

}
