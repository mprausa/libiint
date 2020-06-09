#pragma once

#include <iint/EllipticKernel.h>

namespace iint {
    class KappaKernel : public EllipticKernel {
        public:
            KappaKernel();

            virtual std::string str() const {
                return "kappa";
            }
    };
}
