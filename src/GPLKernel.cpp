#include <iint/GPLKernel.h>

namespace iint {
    int GPLKernel::start(const arb::Acb &x) const {
        return (x == _a) ? -2 : 0;
    }

    arb::Acb GPLKernel::calc(const arb::Acb &x, int k) {
        if (x == _a) {
            return (k==-2) ? 1 : 0;
        }

        if (k%2) return 0;
        return -(_a-x).pow(-k/2-1);
    }
}
