/*
 *  src/iint/TRat.cpp
 *
 *  Copyright (C) 2020 Mario Prausa
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iint/TRat.h>
#include <iint/utilities.h>
#include <algorithm>

namespace iint {
    TRat::TRat(const std::vector<int> &numer, const std::vector<int> &denom)
        : _numer(numer), _denom(denom), _ode({{-1,9},{0,-1,18,-1}}) {
    }

    int TRat::init(const arb::Acb &x)  {
        if (_points.count(x)) return _points[x].n0;

        long prec = x.default_prec();
        point_data data;

        auto x1 = xeps(x);

        auto thesqrt = (1 - 18*x1 + x1*x1).sqrt();
        auto t = x.contains_zero() ? arb::Acb(0,prec) : (1 - 9*x - thesqrt)/(2*x);
        int f = 2;

        if (x.contains_zero()) {
            _ode.init(x,-2,1,{arb::Acb(.5,prec)});
        } else if (x == 9-4*arb::Acb(5,prec).sqrt()) {
            auto c = -arb::Acb::I*arb::Acb(20,prec).pow(arb::Acb(.25,prec))/x;
            _ode.init(x,1,1,{c});
            f = 1;
        } else {
            _ode.init(x,0,0,{thesqrt/(2*x)});
        }

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
        data.reciprocal = false;
        do_shift(_numer,data.numer);
        do_shift(_denom,data.denom);

        data.nummin = std::distance(data.numer.begin(),
            std::find_if(data.numer.begin(),data.numer.end(),[](const arb::Acb &v) {return !v.contains_zero();}));

        data.denmin = std::distance(data.denom.begin(),
            std::find_if(data.denom.begin(),data.denom.end(),[](const arb::Acb &v) {return !v.contains_zero();}));

        assert(data.nummin == 0 || data.denmin == 0);

        data.n0 = f*(data.nummin - data.denmin);

        if (data.denmin > 0) {  // let's work with the reciprocal
            data.reciprocal = true;
            std::swap(data.numer,data.denom);
            std::swap(data.nummin,data.denmin);
        }

        _points[x] = data;

        return data.n0;
    }

    const arb::Acb &TRat::operator() (const arb::Acb &x, int n) {
        assert(_points.count(x));
        auto &data = _points[x];

        if (n < data.n0) return _zero;
        if (!data.reciprocal) return xrat_expansion(data,n);

        if (n-data.n0 < data.cache.size()) return data.cache[n-data.n0];
        while (n-data.n0 > data.cache.size()) (*this)(x,int(data.cache.size())+data.n0);
        assert(n-data.n0 == data.cache.size());

        arb::Acb res;

        if (n == data.n0) {
            res = 1/xrat_expansion(data,-data.n0);
        } else {
            long prec = data.x.default_prec();
            res = arb::Acb(0,prec);

            for (int k=1; k<=n-data.n0; ++k) {
                res -= xrat_expansion(data,k-data.n0) * data.cache[n-k-data.n0];
            }

            res *= data.cache.front();
        }

        data.cache.push_back(res);
        return data.cache.back();
    }

    const arb::Acb &TRat::xrat_expansion(point_data &data, int n) {
        if (n < data.xrat.size()) return data.xrat[n];
        while (n > data.xrat.size()) xrat_expansion(data,int(data.xrat.size()));

        assert(n == data.xrat.size());

        long prec = data.x.default_prec();

        arb::Acb res(0,prec);

        for (int k=0; k<=n; ++k) {
            res += arb::Acb::fac(k,prec) * rat_expansion(data,k) * bell(data,n,k);
        }

        res /= arb::Acb::fac(n,prec);

        data.xrat.push_back(res);
        return data.xrat.back();
    }

    arb::Acb TRat::t_expansion(const arb::Acb &x, int n) {  // expansion of t in Sqrt[x]
        long prec = x.default_prec();

        if (x.contains_zero()) {
            if (n<2) return arb::Acb(0,prec);
            return -_ode(x,n);
        }

        auto res = -_ode(x,n);

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
        if (n < data.nummin-data.denmin) return arb::Acb(0,data.x.default_prec());
        if (n-data.nummin+data.denmin < data.rat.size()) return data.rat[n-data.nummin+data.denmin];
        while (n-data.nummin+data.denmin > data.rat.size()) rat_expansion(data,data.rat.size()+data.nummin-data.denmin);

        assert(data.rat.size() == n-data.nummin+data.denmin);

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

