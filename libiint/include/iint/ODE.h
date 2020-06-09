#pragma once

#include <arb/Acb.h>
#include <unordered_map>

namespace iint {
    class ODE {
        protected:
            arb::Acb zero{0};
            std::vector<std::vector<int>> _p;
            int _mhat;

            struct point_data {
                long prec;
                int k0;
                int r;
                std::vector<arb::Acb> a;
                std::vector<std::vector<arb::Acb>> dp;
            };

            std::unordered_map<arb::Acb,point_data> _points;
        public:
            ODE(const std::vector<std::vector<int>> &p);

            void init(const arb::Acb &x, int k0, int r, const std::vector<arb::Acb> &a);
            int start(const arb::Acb &x);
            const arb::Acb &operator() (const arb::Acb &x, int k);
        private:
            arb::Acb C(int j, int k, const point_data &data) const;
    };
}
