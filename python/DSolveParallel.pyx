# distutils: language = c++
# distutils: sources = Parallel.cc Faddeeva.cc

import numpy as np
from numpy cimport ndarray
from libc.math cimport sqrt

cdef extern from "Parallel.hh":
    cdef cppclass Parallel:
        Parallel (double, double, double, double, int, double)
        double warm (double, double, double) except +
        double cold (double)
        double ocbz, Omega, va, qshear
        int maxiter
        double tol

cdef class DSolveParallel:

    cdef Parallel *parallel

    def __cinit__ (self, double ocbz,
            double Omega=1., double va=1., double qshear=1.5,
            int maxiter=50, double tol=1.48e-8):
        self.parallel = new Parallel (ocbz, Omega, va, qshear, maxiter, tol)

    cpdef warm (self, double beta,
            ndarray[double,ndim=1] kz, ndarray[double,ndim=1] guess):

        cdef int ikz
        cdef int nkz = guess.shape[0]
        cdef ndarray[double,ndim=1] gamma = np.empty ([nkz], dtype=np.double)
        for ikz in range (nkz):
            try:
                gamma[ikz] = self.parallel.warm (beta, kz[ikz], guess[ikz])
            except RuntimeError as e:
                message = e.message + " ikz = %d, kz = %g"
                raise RuntimeError (message % (ikz, kz[ikz]))
        return gamma

    cpdef cold (self, ndarray[double, ndim=1] kz):

        cdef int ikz
        cdef int nkz = kz.shape[0]
        cdef ndarray[double,ndim=1] gamma = np.empty ([nkz], dtype=np.double)
        for ikz in range (nkz):
            gamma[ikz] = self.parallel.cold (kz[ikz])
        return gamma

    property ocbz:
        def __get__ (self): return self.parallel.ocbz
    property Omega:
        def __get__ (self): return self.parallel.Omega
    property va:
        def __get__ (self): return self.parallel.va
    property qshear:
        def __get__ (self): return self.parallel.qshear
    property maxiter:
        def __get__ (self): return self.parallel.maxiter
    property tol:
        def __get__ (self): return self.parallel.tol

    property kmax:
        def __get__ (self):
            cdef double spin = self.ocbz + 2.*self.Omega
            return sqrt (2.*self.qshear*self.Omega*self.Omega*self.ocbz/spin)

    def __dealloc__ (self):
        del self.parallel
