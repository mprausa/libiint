#pragma once

#include <arb/Acb.h>
#include <typeinfo>

namespace iint {
    class Kernel {
        protected:
            bool _nocache = false;
        public:
            arb::Acb operator() (const arb::Acb &x, int k);
            virtual int start(const arb::Acb &x) const = 0;
            virtual std::string str() const {
                return typeid(*this).name();
            }
        protected:
            virtual arb::Acb calc(const arb::Acb &x, int k) = 0;
    };

    inline std::ostream &operator<<(std::ostream &os, const Kernel &krn) {
        os << krn.str();
        return os;
    }
}
