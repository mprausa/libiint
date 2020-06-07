#pragma once

#include <arb/Acb.h>
#include <typeinfo>
#include <unordered_map>
#include <complex>

namespace iint {
    class Kernel {
        protected:
            bool _nocache = false;

            struct point_data {
                int k0;
                std::vector<arb::Acb> cache;
            };
            std::unordered_map<arb::Acb,point_data> _points;
            std::vector<std::complex<double>> _singularities;
        public:
            arb::Acb operator() (const arb::Acb &x, int k);
            int init(const arb::Acb &x);

            virtual std::string str() const {
                return typeid(*this).name();
            }
        protected:
            virtual int _init(const arb::Acb &x) = 0;
            virtual arb::Acb _calc(const arb::Acb &x, int k) = 0;
    };

    inline std::ostream &operator<<(std::ostream &os, const Kernel &krn) {
        os << krn.str();
        return os;
    }
}
