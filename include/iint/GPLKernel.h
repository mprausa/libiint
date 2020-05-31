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

            virtual int start(const arb::Acb &x) const;

            virtual std::string str() const {
                return "w("+_a.str()+")";
            }
        protected:
            virtual arb::Acb calc(const arb::Acb &x, int k);
    };
}
