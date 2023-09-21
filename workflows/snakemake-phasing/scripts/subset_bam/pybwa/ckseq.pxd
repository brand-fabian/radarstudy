cdef extern from "bwa/kseq.h":
    """
    # include <zlib.h>
    KSEQ_DECLARE(gzFile)
    """
    ctypedef void* gzFile
    ctypedef struct kseq_t:
        pass

    kseq_t *kseq_init(gzFile fd)
    void kseq_destroy(kseq_t *ks)
    int kseq_read(kseq_t *seq)