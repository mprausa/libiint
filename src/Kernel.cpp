#include <iint/Kernel.h>

namespace iint {
    arb::Acb Kernel::operator() (const arb::Acb &x, int k) {
        assert(_nocache);   //TODO
        return calc(x,k);
    }
}

