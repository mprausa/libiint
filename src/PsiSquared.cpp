#include <iint/PsiSquared.h>
#include <acb_elliptic.h>

namespace {
    arb::Acb ellipticK(const arb::Acb &x) {
        arb::Acb res;
        acb_elliptic_k(res.get(),x.get(),x.default_prec());
        res.update_default_prec(x.default_prec());
        return res;
    }

    arb::Acb ellipticE(const arb::Acb &x) {
        arb::Acb res;
        acb_elliptic_e(res.get(),x.get(),x.default_prec());
        res.update_default_prec(x.default_prec());
        return res;
    }

}

namespace iint {
    PsiSquared::PsiSquared() : _phi({
        {2,58,58,2},
        {-1,57,-10,-38,-9,1},
        {0,-3,87,-198,150,-39,3},
        {0,0,-1,21,-58,58,-21,1}}) {}

    int PsiSquared::init(const arb::Acb &x) {
        const long prec = x.default_prec();

        if (x == 1) {
            assert(false);  //TODO
        } else {
            arb::Acb zero(0,prec);
            arb::Acb sqrt2 = arb::Acb(2,prec).sqrt();
            arb::Acb thesqrt = (1 - 18*x + x*x).sqrt().conj();
            arb::Acb t = (1 - 9*x - thesqrt)/(2*x);
            arb::Acb z = (t*(4 + t).pow(5))/((4 + 6*t + t*t).pow(2)*(20 + 8*t + t*t));
            arb::Acb sqrtz = z.sqrt();
            arb::Acb lam = 2*sqrtz/(1+sqrtz);
            arb::Acb ellK = ellipticK(lam);
            arb::Acb ellE = ellipticE(lam);
            arb::Acb psi = 2*sqrt2 * ellK/(1+sqrtz).sqrt();   // Sqrt[2]*Pi*Hypergeometric2F1[1/4,3/4,1,z]
            arb::Acb dpsi = (-ellE - (sqrtz - 1)*ellK)/(sqrt2*(sqrtz - 1)*(1+sqrtz).sqrt()*z);
            arb::Acb d2psi = ((-4+8*z)*ellE + (sqrtz-1)*(5*z-4)*ellK)/(4*sqrt2*(1-sqrtz).pow(2) * (1+sqrtz).sqrt().pow(3) * z*z);

            arb::Acb c0 = (1-x)*(3+3*x+2*thesqrt)/(1+18*x+x*x) * psi*psi;

            arb::Acb c1 = (psi*(-16*x*(1-thesqrt + x*(-11 + 2*thesqrt + x*(-61 + 71*x + 39*thesqrt)))*psi -
                          (20480*x.pow(5)*(1 - thesqrt + x*(-20 + 11*thesqrt + x*(70 - 11*thesqrt +
                           x*(-20 + x + thesqrt))))*dpsi)/((-1 + x)*(-1 + thesqrt + x*(12 + 5*x - 3*thesqrt)).pow(2))))/
                           (4*x*thesqrt*(-1 + thesqrt + x*(12 + 5*x - 3*thesqrt)).pow(2));

            arb::Acb c2 = (1024*((-1 + x).pow(3)*(-1 + 9*x + thesqrt).pow(2)*(-1 + thesqrt + x*(12 + 5*x - 3*thesqrt)).pow(4)*
                          (11*(-1 + thesqrt) + x*(330 - 231*thesqrt + x*(-2577 + 938*thesqrt + x*(2460 + 702*thesqrt +
                          x*(3987 + 1131*thesqrt + x*(1050 - 119*x + 9*thesqrt))))))*psi.pow(2) + 320*(-1 + x)*x.pow(4)*thesqrt*
                          (-1 + x + thesqrt).pow(4)*(-1 + 9*x + thesqrt)*(-1 + thesqrt + x*(12 + 5*x - 3*thesqrt)).pow(2)*
                          (1 - thesqrt + x*(-11 + 2*thesqrt + x*(-61 + 71*x + 39*thesqrt)))*psi*dpsi - 640*(-1 + x)*x.pow(3)*
                          (-1 + x + thesqrt).pow(3)*(-1 + thesqrt + x*(12 + 5*x - 3*thesqrt)).pow(2)*(2 - 2*thesqrt + x*(-90 +
                          72*thesqrt + x*(1407 - 839*thesqrt + x*(-8852 + 3461*thesqrt + x*(19045 - 3616*thesqrt +
                          x*(-2*(3333 + 59*thesqrt) + x*(481 - 1543*thesqrt + x*(-232 + 25*x + 25*thesqrt))))))))*psi*dpsi +
                          12800*x.pow(8)*thesqrt*(-1 + x + thesqrt).pow(8)*(-1 + 9*x + thesqrt)*dpsi.pow(2) + 12800*x.pow(8)*
                          thesqrt*(-1 + x + thesqrt).pow(8)*(-1 + 9*x + thesqrt)*psi*d2psi))/((-1 + x).pow(3)*thesqrt.pow(3)*
                          (-1 + 9*x + thesqrt).pow(2)*(2 - 10*x.pow(2) - 2*thesqrt + 6*x*(-4 + thesqrt)).pow(7));

            _phi.init(x,0,0,{c0,zero,c1,zero,c2});

#if 0
            for (int k=0; k<=20; ++k) {
                std::cout << "a[" << k << "] = " << _phi(x,k) << std::endl;
            }
#else
            {
                arb::Acb res;
                arb::Acb del(.25,prec);
                arb::Acb sqrtdel = del.sqrt();

                for (int k=0; k<=200; ++k) {
                    res += _phi(x,k) * sqrtdel.pow(k);
                }

                std::cout << "res = " << res << std::endl;
            }
#endif

            return 0;
        }
    }
}
