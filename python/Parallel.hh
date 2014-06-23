#include <complex>

class Parallel
{
    public:
        Parallel (const double, const double, const double,
                const double, const int, const double);
        double warm (const double, const double, double);
        double cold (const double);
        const double ocbz, Omega, va, qshear;
        const int maxiter;
        const double tol;
    private:
        const double oc2, spin, og, Tide;
};
