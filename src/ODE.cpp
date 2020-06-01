#include <iint/ODE.h>

namespace {
        static arb::Acb pochhammer(const arb::Acb &x, int k) {
            arb::Acb v=1;
            for (int n=0; n<k; ++n) {
                v *= x+n;
            }
            return v;
        };
}

namespace iint {
    ODE::ODE(const std::vector<std::vector<int>> &p) : _p(p) {
        _mhat = 0;
        for (auto v : _p) {
            if (int(v.size())-1 > _mhat) _mhat = int(v.size())-1;
        }
    }

    void ODE::init(const arb::Acb &x, int k0, int r, const std::vector<arb::Acb> &a) {
        if (_points.count(x)) return;

        const int nhat = _p.size()-1;
        const long prec = x.default_prec();
        std::vector<std::vector<arb::Acb>> dp;

        for (auto &p : _p) {
            dp.push_back({});

            for (int n=0; n<=nhat+_mhat; ++n) {
                arb::Acb xp(1,prec);
                arb::Acb v;
                for (int m=n,mhat = p.size(); m<mhat; ++m) {
                    if (p[m]) v += pochhammer(arb::Acb(-m,prec),n) * ((n%2)?-1:1) * xp * p[m];
                    xp *= x;
                }
                dp.back().push_back(v);
            }
        }

        _points[x] = {prec,k0,r,a,dp};
    }

    const arb::Acb &ODE::operator() (const arb::Acb &x, int k) {
        assert(_points.count(x));
        auto &data = _points[x];

        if (k < data.k0) return zero;
        if (k-data.k0 < data.a.size()) return data.a[k-data.k0];
        while(k-data.k0 > data.a.size()) (*this)(x,data.k0+data.a.size());

        assert(k == data.k0+data.a.size());

        const int nhat = _p.size()-1;

        arb::Acb v;
        for (int j=std::max(-_mhat,-(k-data.k0)/2 + nhat - data.r); j<=nhat-data.r-1; ++j) {
            v -= data.a[k-data.k0-2*nhat+2*data.r+2*j]*C(j,k-2*nhat+2*data.r,data);
        }

        v /= C(nhat-data.r,k-2*nhat+2*data.r,data);

        data.a.push_back(v);
        return data.a.back();
    }

    arb::Acb ODE::C(int j, int k, const point_data &data) const {
        const int nhat = _p.size()-1;

        arb::Acb v;
        for (int n=std::max(0,j), n1=std::min(nhat,_mhat+j); n <= n1; ++n) {
            v += ((n%2)?-1:1) * pochhammer(-arb::Acb(k,data.prec)/2-j,n)/arb::Acb::fac(n-j,data.prec) * data.dp[n][n-j];
        }

        return v;
    }
}
