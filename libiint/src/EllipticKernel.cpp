#include <iint/EllipticKernel.h>
#include <iint/utilities.h>
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
    EllipticKernel::EllipticKernel(const std::vector<int> &numer, const std::vector<int> &denom)
        : _trat(numer,denom), _phi({
            {2,58,58,2},
            {-1,57,-10,-38,-9,1},
            {0,-3,87,-198,150,-39,3},
            {0,0,-1,21,-58,58,-21,1}}) {}

    arb::Acb EllipticKernel::_calc(const arb::Acb &x, int k) {
        int n0 = _trat.init(x);
        int k0 = _phi.start(x);

        if (k < n0+k0) return arb::Acb(0,x.default_prec());

        arb::Acb res;
        for (int n=n0; n<=k-k0; ++n) {
            res += _trat(x,n) * _phi(x,k-n);
        }

        return res*arb::Acb::I*arb::Acb::Pi(x.default_prec());
    }

    int EllipticKernel::_init(const arb::Acb &x) {
        const long prec = x.default_prec();
        arb::Acb zero(0,prec);

        int n0 = _trat.init(x);

        if (x.contains_zero()) {
            auto pi2 = arb::Acb::Pi(prec).pow(2);
            _phi.init(x,0,2,{10*pi2,zero,20*pi2,zero,220*pi2});
            return n0;
        } else if (x == 9-4*arb::Acb(5,prec).sqrt()) {
            arb::Acb sqrt2 = arb::Acb(2,prec).sqrt();
            arb::Acb sqrt5 = arb::Acb(5,prec).sqrt();
            arb::Acb z = (9 + 4*sqrt5)/18;
            arb::Acb sqrtz = z.sqrt();
            arb::Acb lam = 2*sqrtz/(1+sqrtz);
            arb::Acb ellK = ellipticK(lam);
            arb::Acb ellE = ellipticE(lam);
            arb::Acb psi = 2*sqrt2 * ellK/(1+sqrtz).sqrt();   // Sqrt[2]*Pi*Hypergeometric2F1[1/4,3/4,1,z]
            arb::Acb dpsi = (-ellE - (sqrtz - 1)*ellK)/(sqrt2*(sqrtz - 1)*(1+sqrtz).sqrt()*z);  // D[psi,z]
            arb::Acb d2psi = ((-4+8*z)*ellE + (sqrtz-1)*(5*z-4)*ellK)/(4*sqrt2*(1-sqrtz).pow(2) * (1+sqrtz).sqrt().pow(3) * z*z);   // D[psi,{z,2}]

            auto c0 = 2*arb::Acb(5,prec).sqrt()/3 * psi.pow(2);
            auto c1 = -arb::Acb::I*sqrt2*sqrt5.sqrt()*psi*(36*(2 + sqrt5)*psi - (5 + 2*sqrt5)*dpsi)/81;
            auto c2 = (-6642*(9 + 4*sqrt5)*psi.pow(2) - 5*(9 + 4*sqrt5)*dpsi.pow(2) +
                      psi*(792*(20 + 9*sqrt5)*dpsi - 5*(9 + 4*sqrt5)*d2psi))/8748;

            _phi.init(x,0,1,{c0,c1,c2});

            return n0;
        } else if (x == 1) {
            auto g = arb::Acb::Gamma(.25,prec).pow(4);
            auto pi = arb::Acb::Pi(prec);

            auto c0 = (-15 + 5*arb::Acb::I)*g/(8*pi);
            auto c1 = -25*arb::Acb::I*pi;
            auto c2 = (75 + 25*arb::Acb::I)*pi.pow(3)/g + (45 - 15*arb::Acb::I)*g/(32*pi);

            _phi.init(x,3,3,{c0,c1,c2});

            return n0+3;
        } else {
            auto x1 = xeps(x);
            arb::Acb sqrt2 = arb::Acb(2,prec).sqrt();
            arb::Acb thesqrt = (1 - 18*x1 + x1*x1).sqrt();

            arb::Acb t = (1 - 9*x1 - thesqrt)/(2*x1);
            arb::Acb z = (t*(4 + t).pow(5))/((4 + 6*t + t*t).pow(2)*(20 + 8*t + t*t));
            arb::Acb sqrtz = z.sqrt();

            arb::Acb lam1 = 2*sqrtz/(1+sqrtz);
            arb::Acb ellK1 = ellipticK(lam1);
            arb::Acb ellE1 = ellipticE(lam1);

            arb::Acb psi,dpsi,d2psi;

            if (x.real() > -9+4*arb::Acb(5,prec).sqrt() && x.abs() < .1 && (z.real() <= 1 || z.imag() > 0)) {
                psi = 2*sqrt2 * ellK1/(1+sqrtz).sqrt();   // Sqrt[2]*Pi*Hypergeometric2F1[1/4,3/4,1,z]
                dpsi = (-ellE1 - (sqrtz - 1)*ellK1)/(sqrt2*(sqrtz - 1)*(1+sqrtz).sqrt()*z);
                d2psi = ((-4+8*z)*ellE1 + (sqrtz-1)*(5*z-4)*ellK1)/(4*sqrt2*(1-sqrtz).pow(2) * (1+sqrtz).sqrt().pow(3) * z*z);
            } else {
                arb::Acb sqrtz2 = (1-z).sqrt();
                arb::Acb lam2 = 2*sqrtz2/(1+sqrtz2);
                arb::Acb ellK2 = ellipticK(lam2);
                arb::Acb ellE2 = ellipticE(lam2);

                psi = 2*sqrt2*ellK1/(1+sqrtz).sqrt() + (4*arb::Acb::I*ellK2)/(1 + sqrtz2).sqrt();   // Sqrt[2]*Pi*Hypergeometric2F1[1/4,3/4,1,z] + 2*I*Pi*Hypergeometric2F1[1/4,3/4,1,1-z]
                dpsi = -arb::Acb::I*(ellE2 + ellK2*(-1 + sqrtz2))/((-1 + sqrtz2)*(1 + sqrtz2).sqrt()*(-1 + z)) -
                       (ellE1 + ellK1*(-1 + sqrtz))/(sqrt2*(-1 + sqrtz)*(1 + sqrtz).sqrt()*z);
                d2psi = (sqrt2*(ellK1*(-1 + sqrtz)*(-4 + 5*z) + ellE1*(-4 + 8*z)))/
                        (8*(-1 + sqrtz).pow(2)*(1 + sqrtz).sqrt().pow(3)*z.pow(2)) +
                        ((2*arb::Acb::I)*(ellK2*(-1 + sqrtz2 + 5*z - 5*sqrtz2*z) +
                        4*ellE2*(-1 + sqrtz2*(-1 + sqrtz2*(2 + sqrtz2) + z))))/
                        (8*(-1 + sqrtz2).pow(2)*(1 + sqrtz2).sqrt().pow(3)*(-1 + z).pow(2));

            }

            arb::Acb c0 = (1-x1)*(3+3*x1+2*thesqrt)/(1+18*x1+x1*x1) * psi*psi;

            arb::Acb c1 = (psi*(-16*x1*(1-thesqrt + x1*(-11 + 2*thesqrt + x1*(-61 + 71*x1 + 39*thesqrt)))*psi -
                          (20480*x1.pow(5)*(1 - thesqrt + x1*(-20 + 11*thesqrt + x1*(70 - 11*thesqrt +
                           x1*(-20 + x1 + thesqrt))))*dpsi)/((-1 + x1)*(-1 + thesqrt + x1*(12 + 5*x1 - 3*thesqrt)).pow(2))))/
                           (4*x1*thesqrt*(-1 + thesqrt + x1*(12 + 5*x1 - 3*thesqrt)).pow(2));

            arb::Acb c2 = (1024*((-1 + x1).pow(3)*(-1 + 9*x1 + thesqrt).pow(2)*(-1 + thesqrt + x1*(12 + 5*x1 - 3*thesqrt)).pow(4)*
                          (11*(-1 + thesqrt) + x1*(330 - 231*thesqrt + x1*(-2577 + 938*thesqrt + x1*(2460 + 702*thesqrt +
                          x1*(3987 + 1131*thesqrt + x1*(1050 - 119*x1 + 9*thesqrt))))))*psi.pow(2) + 320*(-1 + x1)*x1.pow(4)*thesqrt*
                          (-1 + x1 + thesqrt).pow(4)*(-1 + 9*x1 + thesqrt)*(-1 + thesqrt + x1*(12 + 5*x1 - 3*thesqrt)).pow(2)*
                          (1 - thesqrt + x1*(-11 + 2*thesqrt + x1*(-61 + 71*x1 + 39*thesqrt)))*psi*dpsi - 640*(-1 + x1)*x1.pow(3)*
                          (-1 + x1 + thesqrt).pow(3)*(-1 + thesqrt + x1*(12 + 5*x1 - 3*thesqrt)).pow(2)*(2 - 2*thesqrt + x1*(-90 +
                          72*thesqrt + x1*(1407 - 839*thesqrt + x1*(-8852 + 3461*thesqrt + x1*(19045 - 3616*thesqrt +
                          x1*(-2*(3333 + 59*thesqrt) + x1*(481 - 1543*thesqrt + x1*(-232 + 25*x1 + 25*thesqrt))))))))*psi*dpsi +
                          12800*x1.pow(8)*thesqrt*(-1 + x1 + thesqrt).pow(8)*(-1 + 9*x1 + thesqrt)*dpsi.pow(2) + 12800*x1.pow(8)*
                          thesqrt*(-1 + x1 + thesqrt).pow(8)*(-1 + 9*x1 + thesqrt)*psi*d2psi))/((-1 + x1).pow(3)*thesqrt.pow(3)*
                          (-1 + 9*x1 + thesqrt).pow(2)*(2 - 10*x1.pow(2) - 2*thesqrt + 6*x1*(-4 + thesqrt)).pow(7));

            _phi.init(x,0,0,{c0,zero,c1,zero,c2});
            return n0;
        }
    }
}
