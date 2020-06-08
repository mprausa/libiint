#pragma once

#include <iint/Kernel.h>

namespace iint {
    class GPLKernel : public Kernel {
        protected:
            arb::Acb _a;
        public:
            GPLKernel(const arb::Acb &a) : _a(a) {
                _nocache = true;
            }

            virtual std::string str() const {
                return "omega("+_a.str()+")";
            }
        protected:
            virtual arb::Acb _calc(const arb::Acb &x, int k);
            virtual int _init(const arb::Acb &x);
    };
}
