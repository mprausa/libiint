#include <iint/TauKernel.h>

namespace iint {
    TauKernel::TauKernel() : _reciprocal({0,-20,0,1},{2000,1700,560,85,5}) {}

    int TauKernel::_init(const arb::Acb &x) {
        return -_reciprocal.init(x);
    }

    arb::Acb TauKernel::_calc(const arb::Acb &x, int k) {
        auto &data = _points[x];

        assert(k >= data.k0);
        if (k == data.k0) return arb::Acb::I * arb::Acb::Pi(x.default_prec()) / _reciprocal(x,-data.k0);

        arb::Acb res;
        for (int n=1; n<=k-data.k0; ++n) {
            res -= _reciprocal(x,n-data.k0) * data.cache[k-n-data.k0];
        }

        res /= _reciprocal(x,-data.k0);

        //std::cout << "TauKernel -- k = " << k << " res = " << res << std::endl;
        return res;
    }
}
