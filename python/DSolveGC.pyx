# distutils: language = c++
# distutils: sources = Faddeeva.cc

DEF TEST = 0

from libc.math cimport sqrt, cos, sin, exp
from libc.math cimport M_PI as pi

cdef extern from "Faddeeva.hh" namespace "Faddeeva":
    double erfcx (double)

import numpy as np
from numpy cimport ndarray

cdef class DSolveGC:

    cdef double beta, psi, va, Omega, qshear, tol
    cdef int maxiter

    def __init__ (self, double beta, double psi,
            double va = 1., double Omega = 1., double qshear = 1.5,
            int maxiter = 50, double tol = 1.48e-8):

        self.beta = beta
        self.psi = psi
        self.va = va
        self.Omega = Omega
        self.qshear = qshear
        self.maxiter = maxiter
        self.tol = tol

    def __call__ (self, double kx, double kz, double gamma):

        drel = DispersionRelation (kx, kz, self.beta, self.psi,
                self.va, self.Omega, self.qshear)

        cdef double f, df, diff, gamma1
        cdef int it

        for it in range (self.maxiter):
            f, df = drel (gamma)
            IF TEST:
                print "Test passed."
                raise SystemExit
            diff = f/df
            gamma1 = gamma - diff
            if abs (diff) < self.tol:
                return gamma1
            gamma = gamma1

        msg = "Failed to converge after %d iterations, value is %s"
        raise RuntimeError (msg % (self.maxiter, gamma))

cdef class DispersionRelation:

    cdef double kx, kz, beta, psi, va, Omega, qshear
    cdef double by, bz
    cdef double vt, kvt
    cdef double kx2, kz2
    cdef double va2, vt2
    cdef double by2, bz2

    def __init__ (self, double kx, double kz, double beta, double psi,
            double va = 1., double Omega = 1., double qshear = 1.5):

        IF TEST:
            beta = 182./17.
            psi = -46./88.
            kx = -43./97.
            kz = 37./215.

        # Arguments
        self.kx = kx
        self.kz = kz
        self.beta = beta
        self.psi = psi
        self.va = va
        self.Omega = Omega
        self.qshear = qshear

        # Derived constants
        self.bz = cos (psi)
        self.by = -sin (psi)

        self.vt = sqrt (beta)
        self.kvt = kz*self.bz*self.vt

        self.kx2 = kx*kx
        self.kz2 = kz*kz

        self.va2 = va*va
        self.vt2 = self.vt*self.vt

        self.by2 = self.by*self.by
        self.bz2 = self.bz*self.bz

    def __call__ (self, double gamma):

        IF TEST:
            gamma = 37./117.

        # Plasma dispersion function W (zeta)
        cdef double xi = gamma/self.kvt
        cdef double W = 1. - sqrt (pi)*xi*erfcx (xi)
        # Derivative with respect to gamma, *not* xi
        cdef double dW = (2*xi*W + (W-1.)/xi)/self.kvt

        if TEST:
            assert "%.12f" % W == "0.371635652233"
            assert "%.12f" % dW == "-0.999537112231"

        cdef double xx = \
                self.kz2*self.va2 \
                + gamma*gamma \
                + (1.-W)*self.by2*self.kz2*self.vt2
        cdef double xy = \
                self.kx*self.kz*self.by*self.va2 \
                - 2*gamma*self.Omega*(self.bz + self.by2*W/self.bz) \
                + self.vt2*(1.-W)*self.by*self.kx*self.kz
        cdef double yx = \
                self.kx*self.kz*self.by*self.va2 \
                + 2*gamma*self.Omega*(self.bz + self.by2*W/self.bz) \
                + self.vt2*(1.-W)*self.by*self.kx*self.kz
        cdef double yy = \
                (self.kz2*self.bz2 + self.kx2)*self.va2 \
                + gamma*gamma \
                - 2*self.qshear*self.Omega*self.Omega \
                + self.vt2*(1.-W)*self.kx2 \
                + 8*self.by2*W*xi*xi*self.Omega*self.Omega

        IF TEST:
            assert "%.12f" % xx == "0.179281184665"
            assert "%.12f" % xy == "-0.909936408580"
            assert "%.12f" % yx == "0.321329137162"
            assert "%.12f" % yy == "-1.047927427130"

        cdef double dxx = \
                2*gamma - dW*self.by2*self.kz2*self.vt2
        cdef double dxy = \
                -2*self.Omega*(self.bz + self.by2*(W + gamma*dW)/self.bz) \
                - self.vt2*dW*self.by*self.kx*self.kz
        cdef double dyx = \
                +2*self.Omega*(self.bz + self.by2*(W + gamma*dW)/self.bz) \
                - self.vt2*dW*self.by*self.kx*self.kz
        cdef double dyy = \
                2*gamma - self.vt2*dW*self.kx2 \
                + 8*self.by2*(dW*xi*xi + 2*xi*W/self.kvt) \
                *self.Omega*self.Omega

        IF TEST:
            assert "%.12f" % dxx == "0.711469245220"
            assert "%.12f" % dxy == "-2.172439810162"
            assert "%.12f" % dyx == "1.357313045215"
            assert "%.12f" % dyy == "3.866967770051"

        cdef double f = xx*yy - xy*yx
        cdef double df = dxx*yy + xx*dyy - dxy*yx - xy*dyx

        IF TEST:
            assert "%.12f" % f == "0.104515410463"
            assert "%.12f" % df == "1.880843194667"

        return (f, df)
