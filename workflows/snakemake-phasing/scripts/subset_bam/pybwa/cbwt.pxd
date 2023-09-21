from libc.stdint cimport *

cdef extern from "bwa/bwt.h":
    ctypedef uint64_t bwtint_t

    ctypedef struct bwt_t:
        bwtint_t primary        # S^{-1}(0), or the primary index of BWT
        bwtint_t L2[5]          # C(), cumulative count
        bwtint_t seq_len        # sequence length
        bwtint_t bwt_size       # size of bwt, about seq_len / 4
        uint32_t *bwt           # BWT
        uint32_t cnt_table[256] # Occurence array, separated to two parts
        # Suffix array
        int sa_intv
        bwtint_t n_sa
        bwtint_t *sa
