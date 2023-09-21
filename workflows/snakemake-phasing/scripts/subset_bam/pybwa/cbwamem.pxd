from libc.stdint cimport *
from cbwt cimport bwt_t
from cbntseq cimport bntseq_t
from cbwa cimport bwaidx_t, bseq1_t


cdef extern from "bwa/bwamem.h":
    ctypedef struct mem_opt_t:
        int a, b                 # match score and mismatch penalty
        int o_del, e_del
        int o_ins, e_ins
        int pen_unpaired         # phread-scaled penalty for unpaired reads
        int pen_clip5, pen_clip3 # clipping penalty. This score is not deducted from the DP score
        int w                    # band width
        int zdrop                # Z-dropoff
        uint64_t max_mem_intv
        int T                    # output score threshold; only affecting output
        int flag                 # see MEM_F* macros
        int min_seed_len         # minimum seed length
        int min_chain_weight 
        int max_chain_extend 
        float split_factor       # split into a seed if MEM is longer than min_seed_len*split_factor
        int split_width          # split into a seed if its occurence is smaller than this value
        int max_occ              # skip a seed if its occurence is larger than this value
        int max_chain_gap        # do not chain seed if it is max_chain_gap-bp away from the closest seed
        int n_threads            # number of threads
        int chunk_size           # process chunk_size-bp sequences in a batch
        float mask_level         # regard a hit as redundant if the overlap with another better hit is over mask_level times the min length of the two hits
        float drop_ratio         # drop a chain if its seed coverage is below drop_ratio times the seed coverage of a better chain overlapping with the small chain
        float XA_drop_ratio      # when counting hits for the XA tag, ignore alignments with score < XA_drop_ratio * max_score  only effective for the XA tag
        float mask_level_redun 
        float mapQ_coef_len 
        int mapQ_coef_fac 
        int max_ins              # when estimating insert size distribution, skip pairs with insert longer than this value
        int max_matesw           # perform maximally max_matesw rounds of mate-SW for each end
        int max_XA_hits, max_XA_hits_alt   # if there are max_hits or fewer, output them all
        int8_t mat[25]           # scoring matrix; mat[0] == 0 if unset

    ctypedef struct mem_alnreg_t:
        int64_t rb, re # [rb,re): reference sequence in the alignment
        int qb, qe     # [qb,qe): query sequence in the alignment
        int rid        # reference seq ID
        int score      # best local SW score
        int truesc     # actual score corresponding to the aligned region; possibly smaller than $score
        int sub        # 2nd best SW score
        int alt_sc
        int csub       # SW score of a tandem hit
        int sub_n      # approximate number of suboptimal hits
        int w          # actual band width used in extension
        int seedcov    # length of regions coverged by seeds
        int secondary  # index of the parent hit shadowing the current hit; <0 if primary
        int secondary_all
        int seedlen0   # length of the starting seed
        int n_comp
        bint is_alt # number of sub-alignments chained together
        float frac_rep
        uint64_t hash

    ctypedef struct mem_alnreg_v:
        size_t n, m
        mem_alnreg_t *a

    ctypedef struct mem_pestat_t:
        int low, high   # lower and upper bounds within which a read pair is considered to be properly paired
        bint failed     # non-zero if the orientation is not supported by sufficient data
        double avg, std # mean and stddev of the insert size distribution

    ctypedef struct mem_aln_t:
        int64_t pos
        int rid
        int flag
        bint is_rev
        bint is_alt
        uint8_t mapq
        uint32_t NM
        int n_cigar
        uint32_t *cigar
        char *XA
        int score, sub, alt_sc

    mem_opt_t* mem_opt_init() nogil

    mem_alnreg_v mem_align1(
        const mem_opt_t *opt,
        const bwt_t *bwt,
        const bntseq_t *bns,
        const uint8_t *pac,
        int l_seq,
        const char *seq
    ) nogil

    mem_aln_t mem_reg2aln(
        const mem_opt_t *opt,
        const bntseq_t *bns,
        const uint8_t *pac,
        int l_seq,
        const char *seq,
        const mem_alnreg_t *ar
    ) nogil

    void mem_process_seqs(
        const mem_opt_t *opt,
        const bwt_t *bwt,
        const bntseq_t *bns,
        const uint8_t *pac,
        int64_t n_processed,
        int n,
        bseq1_t *seqs,
        const mem_pestat_t *pes0
    ) nogil