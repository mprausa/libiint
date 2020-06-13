#include <iint/IInt.h>

namespace iint {
    bool verbose = false;

    std::unordered_map<std::pair<IInt::kernels_t,arb::Acb>,std::shared_ptr<IInt>,IInt::args_hasher> IInt::_iints;
    arb::Acb IInt::_zero(0);
    arb::Acb IInt::_one(1);

    IInt::IInt(const kernels_t &kernels, const arb::Acb &x0) : _kernels(kernels), _x0(x0) {
        _constants[x0] = 0;

        if (!kernels.empty()) {
            _subiint = fetch(kernels_t(kernels.begin()+1,kernels.end()),x0);
        }
    }

    void IInt::match(const arb::Acb &x1, const arb::Acb &x2, const arb::Acb &x) {
        if (_kernels.empty()) return;
        if (_constants.count(x2)) return;

        assert(_constants.count(x1));

        _subiint->match(x1,x2,x);

        _constants[x2] = 0;
        _constants[x2] = (*this)(x1,x-x1) - (*this)(x2,x-x2);

        if (verbose) std::cout << "<" << *this << "> new constant for " << x2 << " => " << _constants[x2] << std::endl;
    }

    int IInt::start(const arb::Acb &x) {
        if (_kernels.empty()) return 0;

        int n0 = _subiint ? _subiint->start(x) : 0;
        int k0 = _kernels.front()->init(x);

        int n1 = n0+k0+2;
        return n1>=0 ? 0 : n1;
    }

    int IInt::maxlog(const arb::Acb &x) {
        int m = 0;
        for (auto &k : _kernels) {
            if (k->init(x) <= -2) ++m;
        }
        return m;
    }

    const arb::Acb &IInt::operator() (const arb::Acb &x, int n, int m) {
        if (_kernels.empty()) {
            return (n==0 && m==0) ? _one : _zero;
        }

        if (n == 0 && m == 0) {
            auto it = _constants.find(x);
            assert(it != _constants.end());
            return it->second;
        }

        auto it = _cache.find({x,n,m});
        if (it != _cache.end()) return it->second;

        arb::Acb &res = _cache[{x,n,m}];

        if (n == 0 && m>0) {
            int k0 = _kernels.front()->init(x);
            int n0 = _subiint->start(x);

            res = 0;
            for (int k=k0; k<=-2-n0; ++k) {
                res += (*_kernels.front())(x,k) * (*_subiint)(x,-2-k,m-1);
            }

            res /= arb::Acb(m,x.default_prec());
        } else if (n != 0) {
            int mhat = _subiint->maxlog(x);
            int k0 = _kernels.front()->init(x);
            int n0 = _subiint->start(x);
            long prec = x.default_prec();

            res = 0;
            for (int j=m; j<=mhat; ++j) {
                for (int k=k0; k<=n-n0-2; ++k) {
                    res -= arb::Acb::fac(j,prec) * (arb::Acb(-n,prec)/2).pow(m-j-1) / arb::Acb::fac(m,prec) *
                           (*_kernels.front())(x,k) * (*_subiint)(x,n-k-2,j);
                }
            }
        }

        return res;
    }

    arb::Acb IInt::operator() (const arb::Acb &x, const arb::Acb &delta) {
        int n0 = start(x);
        int mhat = maxlog(x);
        auto sqrtdelta = delta.sqrt();
        auto logdelta = mhat>0 ? delta.log() : arb::Acb::nan;

        arb::Acb res;

        //std::cout << "x = " << x << " delta = " << delta << " sqrtdelta = " << sqrtdelta << std::endl;

        for (int n=n0,cnt=-1; cnt < 4; ++n) {
            arb::Acb add;

            for (int m=0; m<=mhat; ++m) {
                arb::Acb v = (*this)(x,n,m);

                if (v.contains_zero()) continue;
                if (n>0) v *= sqrtdelta.pow(n);
                if (m>0) v *= logdelta.pow(m);

                add += v;
            }

            auto res1 = res + add;

            if (res1 == res) {
                if (cnt >= 0) ++cnt;
            } else {
                cnt = 0;
            }
            res = res1;

            //std::cout << *this << " x = " << x << " delta = " << delta << " res = " << res << " cnt = " << cnt << std::endl;
        }

        return res;
    }

    arb::Acb IInt::operator() (const arb::Acb &x) {
        auto mindiff = arb::Acb::infty;
        auto x0 = arb::Acb::nan;

        for (auto &c : _constants) {
            auto diff = (c.first - x).abs();
            if (diff < mindiff) {
                mindiff = diff;
                x0 = c.first;
            }
        }

        assert(!x0.is_nan());

        return (*this)(x0,x-x0);
    }

    std::shared_ptr<IInt> IInt::fetch(const kernels_t &kernels, const arb::Acb &x0) {
        auto &iint = _iints[{kernels,x0}];
        if (!iint) iint = std::make_shared<IInt>(kernels,x0);
        return iint;
    }
}

