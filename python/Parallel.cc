#include "Parallel.hh"
#include "Faddeeva.hh"

#include <cmath>
#include <limits>
#include <iostream>
#include <stdexcept>

using namespace std;

const double eps = numeric_limits<double>::epsilon ();
const double sqrt_pi = sqrt (M_PI);

Parallel::Parallel
(const double ocbz, const double Omega, const double va, const double qshear,
 const int maxiter, const double tol):
    ocbz (ocbz), Omega (Omega), va (va), qshear (qshear),
    maxiter (maxiter), tol (tol),
    oc2 (ocbz*ocbz),
    spin (ocbz + 2.*Omega),
    og (sqrt (spin*(spin - qshear*Omega))),
    Tide (2.*qshear*Omega*Omega*ocbz/spin)
{
}

double Parallel::warm
(const double beta, const double kz, double gamma)
{
    // Alfvenic and thermal frequencies
    const double kva = kz*va;
    const double kvt = kva*sqrt (beta);

    // Non-trivial element of anisotropy tensor in xy-basis
    const double Qyy = og/spin;

    for (int it = 0; it < maxiter; ++it)
    {
        // Reduced frequency
        const double Re_xi = og/kvt;
        const double Im_xi = gamma/kvt;
        complex<double> xi (Re_xi, Im_xi);

        // Plasma dispersion function
        const complex<double> wofz (Faddeeva::w (xi));
        const double Re_Z = -sqrt_pi*imag (wofz);
        const double Im_Z =  sqrt_pi*real (wofz);

        // Dimensionless growth rate
        const double fac = gamma/kvt;

        // Response tensor
        const double Re_L =  fac*Im_Z;
        const double Im_L = -fac*Re_Z;

        // Dimensionless conductivity tensor K = Q.Lambda.Q in xy-basis
        const double Kxx =      Re_L; const double Kxy =     Im_L*Qyy;
        const double Kyx = -Qyy*Im_L; const double Kyy = Qyy*Re_L*Qyy;

        // Dispersion tensor
        const double Dxx = kva*kva + oc2*Kxx;
        const double Dyy = kva*kva + oc2*Kyy - Tide;
        const double Dxy = oc2*Kxy - ocbz*gamma;
        const double Dyx = oc2*Kyx + ocbz*gamma;

        // Determinant
        const double det = Dxx*Dyy - Dxy*Dyx;

        // Derivative of plasma dispersion function
        // with respect to its argument
        const double Re_xiZ = Re_xi*Re_Z - Im_xi*Im_Z;
        const double Im_xiZ = Re_xi*Im_Z + Im_xi*Re_Z;
        const double Re_dZ = -2.*(1. + Re_xiZ);
        const double Im_dZ = -2.*      Im_xiZ ;

        // Derivative of reponse tensor
        const double Re_dL = (Re_L + fac*fac*Re_dZ)/gamma;
        const double Im_dL = (Im_L + fac*fac*Im_dZ)/gamma;

        // Derivative of dimensionless conductivity tensor
        const double dKxx =      Re_dL; const double dKxy =     Im_dL*Qyy;
        const double dKyx = -Qyy*Im_dL; const double dKyy = Qyy*Re_dL*Qyy;

        // Derivative of dispersion tensor
        const double dDxx = oc2*dKxx;
        const double dDyy = oc2*dKyy;
        const double dDxy = oc2*dKxy - ocbz;
        const double dDyx = oc2*dKyx + ocbz;

        // Derivative of determinant
        const double ddet = Dxx*dDyy + Dyy*dDxx - Dxy*dDyx - Dyx*dDxy;

        // Standard Newton raphson:
        if (abs (ddet) < eps)
        //if (abs (gamma*ddet) < eps*abs (det))
        {
            cout << "Warning: Derivative is too small." << endl;
            return gamma;
        }

        const double diff = det/ddet;
        const double next = gamma - diff;

        if (abs (diff) < tol)
        //if (abs (diff) < tol*abs(gamma))
        {
            return next;
        }

        gamma = next;
    }
    throw runtime_error ("Error: Exceeded maximum number of iterations.");
}

double Parallel::cold (const double kz)
{
    const double kappa2 = 2.*(2. - qshear)*Omega*Omega;
    const double qOmega = qshear*Omega;

    const double kva = kz*va;
    const double kva2 = kva*kva;

    const double a = ocbz*ocbz;
    const double b = kva2*kva2 + kva2*(2.*ocbz - qOmega)*ocbz
        + kappa2*ocbz*ocbz;
    const double c = kva2*(spin - qOmega)*(kva2*spin - 2.*qOmega*Omega*ocbz);

    return sqrt ((sqrt (b*b - 4.*a*c) - b)/(2.*a));
}
