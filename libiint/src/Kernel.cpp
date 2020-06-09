#include <iint/Kernel.h>

namespace iint {
    arb::Acb Kernel::operator() (const arb::Acb &x, int k) {
        if (_nocache) return _calc(x,k);
        assert(_points.count(x));
        auto &data = _points[x];

        if (k < data.k0) return arb::Acb(0,x.default_prec());
        if (k-data.k0 < data.cache.size()) return data.cache[k-data.k0];
        while (k-data.k0 > data.cache.size()) (*this)(x,int(data.cache.size())+data.k0);
        assert(k-data.k0 == data.cache.size());

        auto res = _calc(x,k);
        data.cache.push_back(res);
        return res;
    }

    int Kernel::init(const arb::Acb &x) {
        if (_nocache) return _init(x);

        if (_points.count(x)) return _points[x].k0;

        int k0 = _init(x);
        _points[x] = {k0,{}};
        return k0;
    }
}

