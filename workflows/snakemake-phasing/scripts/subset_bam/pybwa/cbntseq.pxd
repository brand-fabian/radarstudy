from libc.stdint cimport *


cdef extern from "bwa/bntseq.h":
    ctypedef struct bntann1_t:
        int64_t offset
        int32_t len
        int32_t n_ambs
        uint32_t gi
        int32_t is_alt
        char *name
        char *anno

    ctypedef struct bntamb1_t:
        int64_t offset
        int32_t len
        char amb

    ctypedef struct bntseq_t:
        int64_t l_pac
        int32_t n_seqs
        uint32_t seed
        bntann1_t *anns
        int32_t n_holes
        bntamb1_t *ambs
        void *fp_pac

    unsigned char nst_nt4_table[256]