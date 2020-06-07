#pragma once

#include <arb/Acb.h>
#include <array>
#include <vector>

namespace iint {
    class Matching {
        public:
            using mpoint_t = std::array<arb::Acb,3>;
            using mpoints_t = std::vector<mpoint_t>;

            static mpoints_t points3(const arb::Acb &x1, const arb::Acb &x2, const arb::Acb &x3, double q);
    };
}
