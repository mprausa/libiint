#include <iint/Matching.h>

namespace iint {
    Matching::mpoints_t Matching::points3(const arb::Acb &x1, const arb::Acb &x2, const arb::Acb &x3, double q) {
        std::vector<std::array<double,3>> y;

        for(double y0=0,ym=q;;) {
            double y1 = ym/(1.-q);
            if (y1 >= .5) y1 = .5;
            y.push_back({y0,y1,ym});
            if (y1 == .5) break;
            y0 = y1;
            ym = y0*(1.+q);
            if (ym >= .5) ym = y0/2+.25;
        }

        for (size_t n=y.size(); n>0; --n) {
            y.push_back({1.-y[n-1][1],1.-y[n-1][0],1.-y[n-1][2]});
        }

        for (size_t n=0,sz=y.size(); n<sz; ++n) {
            y.push_back({1.+y[n][0],1.+y[n][1],1.+y[n][2]});
        }

        auto b = -x1;
        auto c = -(x1 - 2*x2 + x3)/(2*x2 - 2*x3);
        auto d = -(x1*x2 - 2*x1*x3 + x2*x3)/(2*x2 - 2*x3);

        mpoints_t x;

        auto moebius = [&b,&c,&d](double y) {
            return (-b + d*y)/(1 - c*y);
        };

        for (auto &v : y) {
            x.push_back({moebius(v[0]),moebius(v[1]),moebius(v[2])});
        }

        return x;
    }
}

