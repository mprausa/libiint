#pragma once

#include <iint/EllipticKernel.h>

namespace iint {
    class MuKernel : public EllipticKernel {
        protected:
            int _n;
        public:
            MuKernel(int n);

            virtual std::string str() const {
                return "mu("+std::to_string(_n)+")";
            }
    };
}

