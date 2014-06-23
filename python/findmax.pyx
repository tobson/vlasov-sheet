from numpy cimport ndarray

cpdef findmax (ndarray[double,ndim=1] kz, ndarray[double,ndim=1] gamma):

    cdef int i, imax
    cdef double old, new

    imax = 1
    old = gamma[imax]
    for i in range (2, kz.shape[0]-1):
        new = gamma[i]
        if new > old:
            imax = i
            old = new

    cdef double kz1 = kz[imax - 1]
    cdef double kz2 = kz[imax    ]
    cdef double kz3 = kz[imax + 1]

    cdef double gamma1 = gamma[imax - 1]
    cdef double gamma2 = gamma[imax    ]
    cdef double gamma3 = gamma[imax + 1]

    cdef double kz_opt, gamma_max

    kz_opt = gamma1*(kz2 - kz3)*(kz2 + kz3) \
           + gamma2*(kz3 - kz1)*(kz3 + kz1) \
           + gamma3*(kz1 - kz2)*(kz1 + kz2)

    kz_opt /= 2.*(gamma1*(kz2 - kz3) \
                + gamma2*(kz3 - kz1) \
                + gamma3*(kz1 - kz2))

    gamma_max = \
            gamma1*(kz_opt - kz2)*(kz_opt - kz3)/((kz1 - kz2)*(kz1 - kz3)) \
          + gamma2*(kz_opt - kz3)*(kz_opt - kz1)/((kz2 - kz3)*(kz2 - kz1)) \
          + gamma3*(kz_opt - kz1)*(kz_opt - kz2)/((kz3 - kz1)*(kz3 - kz2))

    return kz_opt, gamma_max
