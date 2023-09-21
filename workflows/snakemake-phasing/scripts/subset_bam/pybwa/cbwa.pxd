from libc.stdint cimport *
from cbntseq cimport bntseq_t
from cbwt cimport bwt_t

cdef extern from "bwa/bwa.h":
    ctypedef struct bwaidx_t:
        bwt_t *bwt
        bntseq_t *bns
        uint8_t *pac

        int is_shm
        int64_t l_mem
        uint8_t *mem

    ctypedef struct bseq1_t:
        int l_seq, id
        char *name
        char *comment
        char *seq
        char *qual
        char *sam

    bseq1_t *bseq_read(int chunk_size, int *n_, void *ks1_, void *ks2_)

    int bwa_idx_build(const char *fa, const char *prefix, int algo_type, int block_size) nogil
    bwaidx_t *bwa_idx_load(const char *hint, int which)
    void bwa_idx_destroy(bwaidx_t *idx)

    void bwa_print_sam_hdr(const bntseq_t *bns, const char *hdr_line)