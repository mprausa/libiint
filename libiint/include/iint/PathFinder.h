#pragma once

#include <vector>
#include <functional>
#include <complex>
#include <arb/Acb.h>

namespace iint {
    class PathFinder {
        protected:
            std::vector<arb::Acb> _sings;
        public:
            using mpoint_t = std::array<arb::Acb,3>;
            using mpoints_t = std::vector<mpoint_t>;

            using path_fun_t = std::function<arb::Acb(double)>;
            PathFinder(const std::vector<arb::Acb> &sings) : _sings(sings) {}

            mpoints_t operator() (const path_fun_t &fun, double q, const arb::Acb &a = arb::Acb::infty) const;

            static mpoints_t euclidean(long prec, double q, const arb::Acb &a = arb::Acb::infty);
            static mpoints_t physical(long prec, double q, const arb::Acb &a = arb::Acb::infty);
        private:
            arb::Acb radius(const arb::Acb &p) const;
            double scan(const path_fun_t &fun, double t, const arb::Acb &d) const;
            mpoints_t convert(const path_fun_t &fun, const std::vector<double> &path) const;
    };
}
