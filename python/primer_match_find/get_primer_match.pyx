# %%cython
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3
from collections import deque
from Bio.Seq import reverse_complement

cpdef getmatch(str p, str seq, int mm, str c):
    """
    p = primer/probe sequence
    s = chromosome sequence
    mm = max allowed mismatches
    """

    cdef:
        Py_ssize_t      l = len(p), ls = len(seq)
        Py_ssize_t      i, j, b, count
        str             s, seqr
    out = deque()
    # print(c)
    for i in range(ls-l+1):

        s = seq[i:i+l]
        b = 0
        for j in range(l):
            if p[j] != s[j]: b+=1
            if b > mm: break
        if b <= mm:
            out.append([p, s, c, str(i+1), '+', str(b)])
    seqr = reverse_complement(seq)
    for i in range(ls-l+1):
        s = seqr[i:i+l]
        b = 0
        for j in range(l):
            if p[j] != s[j]: b+=1
            if b > mm: break
        if b <= mm:
            out.append([p, s, c, str(ls-i), '-', str(b)])
    return(out)
