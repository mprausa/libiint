#include <iint/TRat.h>
#include <algorithm>

namespace iint {
    TRat::TRat(const std::vector<int> &numer, const std::vector<int> &denom)
        : _numer(numer), _denom(denom), _ode({{-1,9},{0,-1,18,-1}}) {
    }

    int TRat::init(const arb::Acb &x)  {
        if (_points.count(x)) return -_points[x].denmin;

        long prec = x.default_prec();
        point_data data;

        arb::Acb thesqrt = (1 - 18*x + x*x).sqrt().conj();
        arb::Acb t = (1 - 9*x - thesqrt)/(2*x);

        assert(!x.contains_zero());         //TODO
        assert(!thesqrt.contains_zero());   //TODO

        _ode.init(x,0,0,{thesqrt/(2*x)});

        auto do_shift = [&t,&prec](const std::vector<int> &poly, std::vector<arb::Acb> &res) {
            int N = poly.size()-1;
            res.resize(N+1);

            for (int i=0; i<=N; ++i) {
                res[i] = arb::Acb(0,prec);
                for (int n=i; n<=N; ++n) {
                    res[i] += arb::Acb(poly[n],prec) * arb::Acb::binom(n,i,prec) * t.pow(n-i);
                }
            }
        };

        data.x = x;
        do_shift(_numer,data.numer);
        do_shift(_denom,data.denom);

        data.nummin = std::distance(data.numer.begin(),
            std::find_if(data.numer.begin(),data.numer.end(),[](const arb::Acb &v) {return !v.contains_zero();}));

        data.denmin = std::distance(data.denom.begin(),
            std::find_if(data.denom.begin(),data.denom.end(),[](const arb::Acb &v) {return !v.contains_zero();}));

        _points[x] = data;

        assert(data.denmin == 0);   //TODO
        return -2*data.denmin;
    }

    const arb::Acb &TRat::operator() (const arb::Acb &x, int n) {
        assert(_points.count(x));
        auto &data = _points[x];

        if (n < -2*data.denmin) return _zero;
        if (n+2*data.denmin < data.cache.size()) return data.cache[n+2*data.denmin];
        while (n+2*data.denmin > data.cache.size()) (*this)(x,int(data.cache.size()) - 2*data.denmin);

        assert(n+2*data.denmin == data.cache.size());

        long prec = data.x.default_prec();

        arb::Acb res(0,prec);

        for (int k=0; k<=n; ++k) {
            res += arb::Acb::fac(k,prec) * rat_expansion(data,k) * bell(data,n,k);
        }

        res /= arb::Acb::fac(n,prec);

        data.cache.push_back(res);
        return data.cache.back();
    }

    arb::Acb TRat::t_expansion(const arb::Acb &x, int n) {  // expansion of t in Sqrt[x]
        long prec = x.default_prec();
        auto res = -_ode(x,n);

        assert(!x.contains_zero());     //TODO

        if (n%2) return res;

        n /= 2;

        if (n >= 0) {
            res += ((n%2) ? -1 : 1) * x.pow(-n-1)/2;
        }
        if (n == 0) {
            res -= arb::Acb(4.5,prec);
        }

        return res;
    }

    arb::Acb TRat::rat_expansion(point_data &data, int n) { // expansion of numer/denom in t
        if (n < -data.denmin) return arb::Acb(0,data.x.default_prec());
        if (n+data.denmin < data.rat.size()) return data.rat[n+data.denmin];
        while (n+data.denmin > data.rat.size()) rat_expansion(data,data.rat.size()-data.denmin);

        assert(data.rat.size() == n+data.denmin);

        arb::Acb v(0,data.x.default_prec());

        if (n+data.denmin < data.numer.size()) v = data.numer[n+data.denmin];

        for (int m=1,max=std::min(int(data.denom.size())-1-data.denmin,n+data.denmin); m<=max; ++m) {
            v -= rat_expansion(data,n-m) * data.denom[m+data.denmin];
        }

        v /= data.denom[data.denmin];

        data.rat.push_back(v);

        return v;
    }

    const arb::Acb &TRat::bell(point_data &data, int n, int k) {
        if (n < data.bell.size() && k < data.bell[n].size()) return data.bell[n][k];

        long prec = data.x.default_prec();
        if (k > n) return _zero;

        assert(n <= data.bell.size());

        if (n == data.bell.size()) data.bell.push_back({});

        assert(k == data.bell[n].size());

        arb::Acb res(0,prec);

        if (n == 0 && k == 0) {
            res = arb::Acb(1,prec);
        } else if (k > 0) {
            for (int i=1; i<=n-k+1; ++i) {
                res += i * arb::Acb::fac(n-1,prec) / arb::Acb::fac(n-i,prec) * bell(data,n-i,k-1) * t_expansion(data.x,i);
            }
        }

        data.bell[n].push_back(res);
        return data.bell[n][k];
    }

    void TRat::print(std::ostream &os) const {
        os << "(";
        bool first = true;
        for (int n=0,sz=_numer.size(); n<sz; ++n) {
            auto &v = _numer[n];
            if (!v) continue;
            if (!first) os << "+";
            os << "(" << v << ")";
            if (n == 1) {
                os << "*t";
            } else if (n > 1) {
                os << "*t^" << n;
            }
            first = false;
        }
        os << ")/(";
        first=true;
        for (int n=0,sz=_denom.size(); n<sz; ++n) {
            auto &v = _denom[n];
            if (!v) continue;
            if (!first) os << "+";
            os << "(" << v << ")";
            if (n == 1) {
                os << "*t";
            } else if (n > 1) {
                os << "*t^" << n;
            }
            first = false;
        }
        os << ")";
    }


}

