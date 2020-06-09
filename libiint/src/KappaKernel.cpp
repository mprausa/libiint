#include <iint/KappaKernel.h>

namespace iint {
    KappaKernel::KappaKernel()
        : EllipticKernel(
            {-2000,-200,215,35},
            {0,-160000,-192000,-92800,-19840,0,992,232,24,1}) {}
}
