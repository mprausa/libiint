#include <iint/GPLKernel.h>

namespace iint {
    arb::Acb GPLKernel::_calc(const arb::Acb &x, int k) {
        if (x == _a) {
            return (k==-2) ? 1 : 0;
        }

        if (k%2) return 0;
        return -(_a-x).pow(-k/2-1);
    }

    int GPLKernel::_init(const arb::Acb &x) {
        return (x == _a) ? -2 : 0;
    }
}
