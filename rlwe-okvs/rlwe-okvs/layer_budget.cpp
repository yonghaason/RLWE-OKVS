#include "layer_budget.h"
#include "sspmt.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace rlweOkvs {

namespace {

double logBinomUpperTail(uint64_t n, double logP, double logQ, double lgn1,
                         uint64_t k, double mean) {
  if (k == 0) return 0.0;
  if (k > n) return -std::numeric_limits<double>::infinity();

  double maxLog = -std::numeric_limits<double>::infinity();
  double sum = 0.0;
  for (uint64_t j = k; j <= n; ++j) {
    const double l = lgn1 - std::lgamma((double)j + 1.0) -
                     std::lgamma((double)(n - j) + 1.0) + (double)j * logP +
                     (double)(n - j) * logQ;
    if (l > maxLog) {
      sum = (maxLog == -std::numeric_limits<double>::infinity())
                ? 1.0
                : sum * std::exp(maxLog - l) + 1.0;
      maxLog = l;
    } else {
      sum += std::exp(l - maxLog);
    }
    if ((double)j > mean && l < maxLog - 80.0) break;
  }
  return maxLog + std::log(sum);
}

}

std::vector<uint64_t> hallCaps(uint64_t n, uint32_t numSlots,
                               uint64_t positionRange, uint32_t b,
                               uint32_t lambda) {
  std::vector<uint64_t> caps;
  if (n == 0 || b == 0) {
    return caps;
  }

  const double logEps = -(double)lambda * std::log(2.0) -
                        std::log((double)numSlots) -
                        2.0 * std::log((double)b);
  const double lgn1 = std::lgamma((double)n + 1.0);

  caps.resize(b);
  uint64_t k = 0;
  for (uint32_t g = 1; g <= b; ++g) {
    const double p = (double)g / (double)positionRange;
    const double logP = std::log((double)g) - std::log((double)positionRange);
    const double logQ = std::log1p(-p);
    const double mean = (double)n * p;

    while (logBinomUpperTail(n, logP, logQ, lgn1, k, mean) > logEps) {
      ++k;
    }
    caps[g - 1] = k + 1;
  }
  return caps;
}

uint64_t budgetForSpan(const std::vector<uint64_t> &caps, uint32_t b,
                       uint32_t W) {
  uint64_t budget = 0;
  for (uint32_t g = 1; g <= b; ++g) {
    budget = std::max(budget, oc::divCeil(caps[g - 1] * (uint64_t)(b + W - 1),
                                          (uint64_t)(g + W - 1)));
  }
  return budget;
}

u32 sspmtParams::resolveSpanBlocks(u64 n) const {
  if (span_blocks) {
    return span_blocks;
  }
  const uint32_t numSlots = heNumSlots;
  const uint64_t m = oc::roundUpTo(bandExpansion * n, numSlots);
  const uint32_t b = (uint32_t)(m / numSlots);
  const auto caps = hallCaps(n, numSlots, m - bandWidth + 1, b,
                             layerBudgetLambda);
  if (caps.empty()) {
    return 1;
  }

  const double floor = (double)caps.back();

  const uint32_t Wlo = std::min<uint32_t>(layerBudgetLambda, b);
  uint32_t span = Wlo;
  double best = std::numeric_limits<double>::max();
  for (uint32_t W = Wlo; W <= b; ++W) {
    const double obj =
        (double)budgetForSpan(caps, b, W) / floor + spanCostRatio * W;
    if (obj < best) {
      best = obj;
      span = W;
    }
  }
  return span;
}

u32 sspmtParams::resolveLayerBudget(u64 n) const {
  if (layerBudget) {
    return layerBudget;
  }
  const uint32_t numSlots = heNumSlots;
  const uint64_t m = oc::roundUpTo(bandExpansion * n, numSlots);
  const uint32_t b = (uint32_t)(m / numSlots);
  const auto caps = hallCaps(n, numSlots, m - bandWidth + 1, b,
                             layerBudgetLambda);
  if (caps.empty()) {
    return 0;
  }
  return (uint32_t)budgetForSpan(caps, b, resolveSpanBlocks(n));
}

}
