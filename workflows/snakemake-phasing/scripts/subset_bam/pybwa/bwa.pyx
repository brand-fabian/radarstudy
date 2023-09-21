import typing
cimport cbwamem
cimport cbwa
cimport ckseq
cimport cbntseq
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from cpython cimport array as c_array
from libc.stdlib cimport free, malloc, realloc
from libc.string cimport memset, memcpy
from libc.stdint cimport *
from pysam.libcalignmentfile cimport AlignedSegment


__package_name__ = "pybwa"
__version__ = "0.1"
__bwa_version__ = "0.7.17"
__bwa_pg__ = f"@PG\tID:{__package_name__}\tPN:bwa\tVN:{__version__}(bwa {__bwa_version__})\tCL:*"


cpdef enum MEM_F:
    PE             = 0x2
    NOPAIRING      = 0x4
    ALL            = 0x8
    NO_MULTI       = 0x10
    NO_RESCUE      = 0x20
    REF_HDR        = 0x100
    SOFTCLIP       = 0x200
    SMARTPE        = 0x400
    PRIMARY5       = 0x800
    KEEP_SUPP_MAPQ = 0x1000
    XB             = 0x2000


cpdef enum BWTALGO:
    AUTO  = 0
    RB2   = 1
    BWTSW = 2
    IS    = 3


cpdef enum BWA_IDX:
    BWT   = 0x1
    BNS   = 0x2
    PAC   = 0x4
    I_ALL = 0x7


cdef extern from "<zlib.h>":
    ctypedef void *gzFile
    ctypedef int64_t z_off_t

    int gzclose(gzFile fp)
    gzFile gzopen(char *path, char *mode)
    int gzread(gzFile fp, void *buf, unsigned int n)
    char *gzerror(gzFile fp, int *errnum)


cdef inline char *copy_string(const char* s, size_t l, bint dupempty):
    cdef char *copy = NULL
    cdef char *terminator = '\0'

    if l > 0 or dupempty:
        copy = <char*> malloc(l + 1)
        if not copy:
            raise MemoryError()

        memcpy(copy, s, l)
        copy[l] = terminator[0]
    return copy


cdef inline char *get_qual_string(arr, size_t l):
    cdef char *copy = NULL
    cdef char *terminator = '\0'

    if l > 0:
        copy = <char*> malloc(l + 1)
        if not copy:
            raise MemoryError()

        for i, c in enumerate(arr):
            copy[i] = c + 0x21
        
        copy[l] = terminator[0]
    return copy


cdef inline str str_or_none(const char* s):
    if s is not NULL:
        try:
            return (<bytes>s).decode('ascii')
        except UnicodeDecodeError:
            return None
    else:
        return None
    

cdef void destroy_seqs(cbwa.bseq1_t *seqs, uint64_t n_seqs):
    if n_seqs > 0 and seqs is not NULL:
        for i in range(n_seqs):
            if seqs[i].name is not NULL:
                free(seqs[i].name)
            if seqs[i].comment is not NULL:
                free(seqs[i].comment)
            if seqs[i].seq is not NULL:
                free(seqs[i].seq)
            if seqs[i].qual is not NULL:
                free(seqs[i].qual)
            if seqs[i].sam is not NULL:
                free(seqs[i].sam)
        free(seqs)


cdef inline void two_bit_to_seq(const char *s, size_t l, char *t):
    """Decode two bit format to a char sequence.
    
    Warning
    -------
    Does not bounds-check the char* target given. Make sure enough
    space is allocated before calling this function.
    """
    cdef size_t pos = 0
    cdef char *terminator = '\0'
    if l > 0:
        while pos < l:
            if s[pos] > 4:
                t[pos] = s[pos]
            else:
                t[pos] = "ACGTN"[<int>s[pos]]
            pos += 1
    t[pos] = terminator[0]


cdef inline char *seq_to_two_bit(const char *s, size_t l):
    cdef size_t pos = 0
    cdef char *terminator = '\0'
    cdef char *res = NULL
    
    if l > 0:
        res = <char*>malloc(l + 1)
        if res == NULL:
            raise MemoryError()
        
        while pos < l:
            res[pos] = cbntseq.nst_nt4_table[<int>s[pos]]
            pos += 1
        res[l] = terminator[0]

    return res


cdef class BwamemOptions:
    """Bwamem algorithm options.

    Wrapper around mem_opt_t*. Can be used to access and change the runtime
    options for the BWA algorithm.

    Notes
    -----
    Some parameters are renamed from the original struct names to more
    descriptive ones, usually found in the respective comments.
    
    Attribute description are copied from original bwa source code at 
    https://github.com/lh3/bwa/blob/3ddd7b87d41f89a404d95f083fb37c369f70d783/bwamem.h#L52.
    """
    cdef cbwamem.mem_opt_t *_opt

    def __cinit__(self):
        self._opt = cbwamem.mem_opt_init()
        if self._opt is NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._opt is not NULL:
            free(self._opt)

    cdef cbwamem.mem_opt_t *raw(self):
        return self._opt

    @staticmethod
    cdef BwamemOptions wrap(cbwamem.mem_opt_t *opt):
        cdef BwamemOptions v = BwamemOptions()
        if v._opt is not NULL:
            free(v._opt)
        v._opt = opt
        return v

    @property
    def match_score(self) -> int:
        """int: Match score (Original: a)."""
        return self._opt.a

    @match_score.setter
    def match_score(self, int val):
        self._opt.a = val

    @property
    def mismatch_penalty(self) -> int:
        """int: Mismatch penalty (Original: b)."""
        return self._opt.b

    @mismatch_penalty.setter
    def mismatch_penalty(self, int val):
        self._opt.b = val

    @property
    def o_del(self) -> int:
        """int: o_del"""
        return self._opt.o_del

    @o_del.setter
    def o_del(self, int val):
        self._opt.o_del = val

    @property
    def e_del(self) -> int:
        """int: e_del"""
        return self._opt.e_del

    @e_del.setter
    def e_del(self, int val):
        self._opt.e_del = val

    @property
    def o_ins(self) -> int:
        """int: o_ins"""
        return self._opt.o_ins

    @o_ins.setter
    def o_ins(self, int val):
        self._opt.o_ins = val

    @property
    def e_ins(self) -> int:
        """int: e_ins"""
        return self._opt.e_ins

    @e_ins.setter
    def e_ins(self, int val):
        self._opt.e_ins = val

    @property
    def pen_unpaired(self) -> int:
        """int: Phred-scaled penalty for unpaired reads."""
        return self._opt.pen_unpaired

    @pen_unpaired.setter
    def pen_unpaired(self, int val):
        self._opt.pen_unpaired = val

    @property
    def pen_clip5(self) -> int:
        """int: Clipping penalty."""
        return self._opt.pen_clip5

    @pen_clip5.setter
    def pen_clip5(self, int val):
        self._opt.pen_clip5 = val

    @property
    def pen_clip3(self) -> int:
        """int: Clipping penalty."""
        return self._opt.pen_clip3

    @pen_clip3.setter
    def pen_clip3(self, int val):
        self._opt.pen_clip3 = val

    @property
    def band_width(self) -> int:
        """int: band width (Original: w)"""
        return self._opt.w

    @band_width.setter
    def band_width(self, int val):
        self._opt.w = val

    @property
    def z_drop(self) -> int:
        """int: Z-dropoff"""
        return self._opt.zdrop

    @z_drop.setter
    def z_drop(self, int val):
        self._opt.zdrop = val

    @property
    def max_mem_intv(self) -> int:
        """int: max_mem_intv"""
        return self._opt.max_mem_intv

    @max_mem_intv.setter
    def max_mem_intv(self, int val):
        if val < 0:
            raise ValueError()
        self._opt.max_mem_intv = val

    @property
    def threshold(self) -> int:
        """int: Output score threshold."""
        return self._opt.T

    @threshold.setter
    def threshold(self, int val):
        self._opt.T = val

    @property
    def flag(self) -> int:
        """int: Flag values. Retrieve individual flag using bitwise-and with MEM_F.*"""
        return self._opt.flag

    def set_flag(self, MEM_F val):
        """Set a flag.
        
        Parameters
        ----------
        val : MEM_F
              Set this flag.
        """
        self._opt.flag |= val

    def unset_flag(self, MEM_F val):
        """Unset a flag.
        
        Parameters
        ----------
        val : MEM_F
              Unset this flag.
        """
        self._opt.flag &= ~val

    def has_flag(self, MEM_F val):
        """Check if a flag is set.
        
        Parameters
        ----------
        val : MEM_F
              A valid flag value (MEM_F enum value).
              
        Returns
        -------
        bool
            True, if the flag is set.
        """
        return (self._opt.flag & val) > 0

    @flag.setter
    def flag(self, int val):
        self._opt.flag = val

    @property
    def min_seed_len(self) -> int:
        """int: Minimum seed length."""
        return self._opt.min_seed_len

    @min_seed_len.setter
    def min_seed_len(self, int val):
        self._opt.min_seed_len = val

    @property
    def min_chain_weight(self) -> int:
        """int: min_chain_weight"""
        return self._opt.min_chain_weight

    @min_chain_weight.setter
    def min_chain_weight(self, int val):
        self._opt.min_chain_weight = val

    @property
    def max_chain_extend(self) -> int:
        """int: max_chain_extend"""
        return self._opt.max_chain_extend

    @max_chain_extend.setter
    def max_chain_extend(self, int val):
        self._opt.max_chain_extend = val

    @property
    def split_factor(self) -> float:
        """float: Split into a seed if MEM is longer than this value."""
        return self._opt.split_factor

    @split_factor.setter
    def split_factor(self, float val):
        self._opt.split_factor = val

    @property
    def split_width(self) -> int:
        """int: Split into a seed if its occurence is smaller than this vaue."""
        return self._opt.split_width

    @split_width.setter
    def split_width(self, int val):
        self._opt.split_width = val

    @property
    def max_occ(self) -> int:
        """int: Skip a seed if its occurence is larger than this value."""
        return self._opt.max_occ

    @max_occ.setter
    def max_occ(self, int val):
        self._opt.max_occ = val

    @property
    def max_chain_gap(self) -> int:
        """int: Do not chain seed if it is max-chain-gap-bp away from the closest seed."""
        return self._opt.max_chain_gap

    @max_chain_gap.setter
    def max_chain_gap(self, int val):
        self._opt.max_chain_gap = val
    
    @property
    def n_threads(self) -> int:
        """int: Number of threads for alignment."""
        return self._opt.n_threads

    @n_threads.setter
    def n_threads(self, int val):
        self._opt.n_threads = val
    
    @property
    def chunk_size(self) -> int:
        """int: Process this many base pairs per batch."""
        return self._opt.chunk_size

    @chunk_size.setter
    def chunk_size(self, int val):
        self._opt.chunk_size = val

    @property
    def mask_level(self) -> float:
        """float: Regard hits as redundant, if the overlap with another better hit is over mask_level times the min length of the two hit."""
        return self._opt.mask_level

    @mask_level.setter
    def mask_level(self, float val):
        self._opt.mask_level = val

    @property
    def drop_ratio(self) -> float:
        """float: Drop a chain, if its seed coverage is below drop_ratio * seed coverage of a better chaiin overlapping with the smal chain."""
        return self._opt.drop_ratio

    @drop_ratio.setter
    def drop_ratio(self, float val):
        self._opt.drop_ratio = val

    @property
    def XA_drop_ratio(self) -> float:
        """float: When counting hits for the XA tag, ignore alignments with score < XA_drop_ratio * max_score"""
        return self._opt.XA_drop_ratio

    @XA_drop_ratio.setter
    def XA_drop_ratio(self, float val):
        self._opt.XA_drop_ratio = val
    
    @property
    def mask_level_redun(self) -> float:
        """float: mask_level_redun"""
        return self._opt.mask_level_redun

    @mask_level_redun.setter
    def mask_level_redun(self, float val):
        self._opt.mask_level_redun = val

    @property
    def mapQ_coef_len(self) -> float:
        """float: mapQ_coef_len"""
        return self._opt.mapQ_coef_len

    @mapQ_coef_len.setter
    def mapQ_coef_len(self, float val):
        self._opt.mapQ_coef_len = val

    @property
    def mapQ_coef_fac(self) -> int:
        """int: mapQ_coef_fac"""
        return self._opt.mapQ_coef_fac

    @mapQ_coef_fac.setter
    def mapQ_coef_fac(self, int val):
        self._opt.mapQ_coef_fac = val

    @property
    def max_ins(self) -> int:
        """int: When estimating insert size distribution, skip pairs with insert longer than this."""
        return self._opt.max_ins

    @max_ins.setter
    def max_ins(self, int val):
        self._opt.max_ins = val

    @property
    def max_matesw(self) -> int:
        """int: Perform this many rounds of mate-SW for each end at max."""
        return self._opt.max_matesw

    @max_matesw.setter
    def max_matesw(self, int val):
        self._opt.max_matesw = val

    @property
    def max_XA_hits(self) -> int:
        """int: Output all XA hits if they are less than this max."""
        return self._opt.max_XA_hits

    @max_XA_hits.setter
    def max_XA_hits(self, int val):
        self._opt.max_XA_hits = val

    @property
    def max_XA_hits_alt(self) -> int:
        """int: Output all XA alt hits if they are less than this max."""
        return self._opt.max_XA_hits_alt

    @max_XA_hits_alt.setter
    def max_XA_hits_alt(self, int val):
        self._opt.max_XA_hits_alt = val

    @property
    def mat(self) -> typing.List[int]:
        """List[int]: (read only) Scoring matrix"""
        return [x for x in self._opt.mat[:25]]

    def __str__(self) -> str:
        return f"BwamemOptions(match_score={self.match_score},mismatch_penalty={self.mismatch_penalty},o_del={self.o_del}," \
            + f"e_del={self.e_del},o_ins={self.o_ins},e_ins={self.e_ins},pen_unpaired={self.pen_unpaired}," \
            + f"pen_clip5={self.pen_clip5},pen_clip3={self.pen_clip3},band_width={self.band_width}," \
            + f"z_drop={self.z_drop},max_mem_intv={self.max_mem_intv},threshold={self.threshold}," \
            + f"flag={self.flag},min_seed_len={self.min_seed_len},min_chain_weight={self.min_chain_weight}," \
            + f"max_chain_extend={self.max_chain_extend},split_factor={self.split_factor}," \
            + f"split_width={self.split_width},max_occ={self.max_occ},max_chain_gap={self.max_chain_gap}," \
            + f"n_threads={self.n_threads},chunk_size={self.chunk_size},mask_level={self.mask_level}," \
            + f"drop_ratio={self.drop_ratio},XA_drop_ratio={self.XA_drop_ratio}," \
            + f"mask_level_redun={self.mask_level_redun},mapq_coef_len={self.mapQ_coef_len}," \
            + f"mapq_coef_fac={self.mapQ_coef_fac},max_ins={self.max_ins},max_matesw={self.max_matesw}," \
            + f"max_XA_hits={self.max_XA_hits},max_XA_hits_alt={self.max_XA_hits_alt},mat={self.mat})"


cdef class Bseq:
    """Read-only sequence object.

    A basic sequence object used by bwa internally to represent reads from
    a fastq file. If the read was aligned to a reference, the sam property
    is set.

    Notes
    -----
    This object can not easily, and should not be modified from python.
    """
    cdef cbwa.bseq1_t *_seq

    def __cinit__(self):
        self._seq = NULL

    def __dealloc__(self):
        """Delete the instance.

        Notes
        -----
        The memory of these objects is not handled by this type, but by the
        list type that generates them.
        """
        pass

    @staticmethod
    cdef Bseq wrap(cbwa.bseq1_t *seq):
        cdef Bseq v = Bseq()
        v._seq = seq
        return v

    cdef cbwa.bseq1_t *raw(self):
        return self._seq

    def is_valid(self) -> bool:
        """Check if this instance is valid.
        
        Checks whether the underlying pointer that this object references
        is valid and can be read.
        
        Returns
        -------
        bool
            True if this object does not reference a null pointer.
        """
        if self._seq is not NULL:
            return True
        else:
            raise False

    @property
    def l_seq(self) -> int:
        """int: Length of the sequence in bases."""
        return self._seq.l_seq if self.is_valid() else -1

    @property
    def id(self) -> int:
        """int: Id of the read in the sequence list."""
        return self._seq.id if self.is_valid() else -1

    @property
    def name(self) -> str:
        """str: Read name. Part of the fastq file containing instrument and flow
           cell name and position of the read on the flow cell."""
        return str_or_none(self._seq.name) if self.is_valid() else None

    @property
    def comment(self) -> str:
        """str: Fastq comment. Typically contains barcode, filter etc."""
        return str_or_none(self._seq.comment) if self.is_valid() else None

    @property
    def seq(self) -> str:
        """str: Sequence bases."""
        cdef char *s = NULL
        if self.is_valid():
            s = <char*>malloc(self._seq.l_seq)
            if s is NULL:
                raise MemoryError()
            two_bit_to_seq(self._seq.seq, self._seq.l_seq, s)
            free(self._seq.seq)
            self._seq.seq = s
        return str_or_none(self._seq.seq) if self.is_valid() else None

    @property
    def qual(self) -> str:
        """str: ASCII encoded phred-scaled base quality scores."""
        return str_or_none(self._seq.qual) if self.is_valid() else None

    @property
    def sam(self) -> str:
        """str: If set, this contains the alignment of the read as sam string.
          Adheres to the sam specification for formatting."""
        return str_or_none(self._seq.sam) if self.is_valid() else None

    def __str__(self) -> str:
        def _append_attr(obj, s, attr):
            if getattr(obj, attr) is not None:
                s += f"{attr}=\"{getattr(obj, attr)}\""
            else:
                s += f"{attr}=None"
            return s

        s = f"{self.__class__.__name__}=("
        if not self.is_valid():
            s += "[invalid]"
        else:
            s += f"l_seq={self.l_seq},id={self.id},"
            s = _append_attr(self, s, "name") + ","
            s = _append_attr(self, s, "comment") + ","
            s = _append_attr(self, s, "seq") + ","
            s = _append_attr(self, s, "qual") + ","
            s = _append_attr(self, s, "sam")
        s += ")"
        
        return s


cdef class SequenceListForwardIterator:
    """Sequence List forward iterator
    
    Forward iterator through a list of sequence (bseq1_t)
    objects.
    """
    cdef cbwa.bseq1_t *_seq_list   # List of sequence vectors
    cdef uint64_t _curr_idx        # Index into _seq_lists
    cdef uint32_t _n_seqs          # Number of sequence vectors
    cdef uint32_t _step            # Step size
    cdef uint32_t *_ref_count      # Ref count to the seq_list

    def __cinit__(self):
        self._seq_list = NULL
        self._ref_count = NULL
        self._curr_idx = 0
        self._n_seqs = 0
        self._step = 1

    def __dealloc__(self):
        self._ref_count[0] -= 1
        if self._ref_count[0] == 0:
            free(self._ref_count)
            destroy_seqs(self._seq_list, self._n_seqs)

    def __next__(self):
        """Get the next Bseq value in the list.
        
        Returns
        -------
        Bseq
            A wrapped bseq1_t sequence instance.
        """
        if self._seq_list is NULL:
            raise ValueError("Invalid Iterator")
        self._curr_idx += self._step
        if self._curr_idx > self._n_seqs:
            raise StopIteration
        return Bseq.wrap(&self._seq_list[self._curr_idx-self._step])


cdef class SequenceListView:
    """Partial view to a sequence list.
    
    A view to a part of a sequence list object.
    """
    cdef cbwa.bseq1_t *_seq_list
    cdef uint64_t _start
    cdef uint64_t _stop
    cdef uint64_t _step
    cdef uint64_t _n_seqs
    cdef bint _is_paired_end
    cdef uint32_t *_ref_count

    def __cinit__(self):
        self._seq_list = NULL
        self._ref_count = NULL
        self._start = 0
        self._stop = 0
        self._step = 1
        self._n_seqs = 0

    def __dealloc__(self):
        self._ref_count[0] -= 1
        if self._ref_count[0] == 0:
            free(self._ref_count)
            destroy_seqs(self._seq_list, self._n_seqs)
            
    @property
    def is_paired_end(self) -> bool:
        """bool: Contains interleaved paired end reads, if true."""
        return self._is_paired_end
            
    def __len__(self) -> int:
        return self._n_seqs

    def __iter__(self):
        """Create an iterator for this view.
        
        Returns
        -------
        SequenceListForwardIterator
            An iterator.
        """
        cdef SequenceListForwardIterator it = SequenceListForwardIterator()
        it._seq_list = self._seq_list
        it._curr_idx = self._start
        it._n_seqs = self._stop
        it._ref_count = self._ref_count
        it._step = self._step
        self._ref_count[0] += 1
        return it

    def __getitem__(self, acc) -> typing.Union[Bseq,SequenceListView]:
        """Get a member or part of the list.
        
        Parameters
        ----------
        acc : Union[int, slice]
              An integer or slice object. Int gives the 0-based index of the
              sequence to be returned. Slice objects will return a new View on
              the data for the defined interval.
        
        Returns
        -------
        Union[Bseq, SequenceListView]
            If an integer was used to access data it returns the Bseq object
            for the corresponding position in the list. Otherwise a view on
            the interval defined by the slice object is returned.
            
        Raises
        ------
        TypeError
            Raises a TypeError if something else than an integer or slice object
            is passed to the function.
        IndexError
            The integer index is out of range.
        """
        cdef SequenceListView slv
        cdef uint64_t start, stop, step
        cdef uint64_t idx = 0
        cdef uint64_t length = self._stop - self._start

        if isinstance(acc, int):
            if acc < 0:
                idx = self._start + length + (self._step * acc)
            else:
                idx = self._start + (self._step * acc)
            if idx < self._start or idx >= self._stop:
                raise IndexError(f"Index out of range [0,{(length+1)//self._step}).")
            return Bseq.wrap(&self._seq_list[idx])
        elif isinstance(acc, slice):
            if acc.start is not None and acc.start > 0:
                start = self._start + (self._step * acc.start)
            elif acc.start is not None and acc.start < 0:
                start = self._start + length + (self._step * acc.start)
            else:
                # acc.start is None or 0
                start = self._start
            
            if acc.stop is None or self._step * acc.stop > self._stop:
                stop = self._stop
            if acc.stop < 0:
                stop = self._start + length + (self._step * acc.stop)
            else:
                stop = self._start + (self._step * acc.stop)

            if acc.step is not None and acc.step > 0:
                step = acc.step
            else:
                step = 1

            slv = SequenceListView()
            slv._seq_list = self._seq_list
            slv._ref_count = self._ref_count
            self._ref_count[0] += 1
            slv._start = start
            slv._stop = stop
            slv._step = self._step * step
            slv._n_seqs = (stop - start) // (self._step * step)
            slv._is_paired_end = self._is_paired_end
            return slv
        else:
            raise TypeError("Type not supported by getitem.")


cdef class _BseqRead:
    """Helper class to fill bseq1_t structs.

    Helper class to fill the bseq1_t struct with the contents of the strings
    set in the instance constructor.
    """
    cdef bytes name
    cdef bytes comment
    cdef bytes sequence
    cdef bytes qual

    def __cinit__(
        self,
        const unsigned char[:] name,
        const unsigned char[:] comment,
        const unsigned char[:] sequence,
        const unsigned char[:] qual
    ):
        self.name = bytes(name)    
        self.comment = bytes(comment)
        self.sequence = bytes(sequence)
        self.qual = bytes(qual)

    def __dealloc__(self):
        pass

    cdef void to_bseq1(self, cbwa.bseq1_t *s):
        s.name = copy_string(self.name, len(self.name), True)
        s.comment = copy_string(self.comment, len(self.comment), False)
        s.seq = copy_string(self.sequence, len(self.sequence), True)
        s.qual = copy_string(self.qual, len(self.qual), False)
        s.sam = copy_string("", 0, False)
        s.l_seq = len(self.sequence)


cdef class SequenceList:
    """A memory view on (un-)aligned sequence objects.
    
    Internally, bwa represents the reads from a given set of fastq files as
    `bseq1_t` struct objects. The SequenceList wraps a `bseq1_t*` (read list)
    and provides convenient python based accessor methods.
    
    Accessor methods return Bseq objects or Iterators on the base memory
    managed by the SequenceList.
    
    Warning
    -------
    If this object and all derived iterators (SequenceListView,
    SequenceListForwardIterator) have gone out of scope all Bseq objects
    returned by this type will be inaccessible as well (the underlying
    memory was deallocated).
    
    Make sure not to access Bseq objects from a SequenceList after it has
    gone out of scope.
    
    Notes
    -----
    Under no circumstances you should (try to) change the underlying sequence
    objects.
    
    The List features a basic ref-count for the base `bseq1_t*` object
    that holds the memory. If the list goes out of scope (is deallocated) and
    all iterators (SequenceListView, SequenceListForwardIterator) the memory
    of the list is deallocated.
    """
    cdef cbwa.bseq1_t *_seq_list
    cdef uint64_t _n_seqs
    cdef bint _is_paired_end
    cdef uint32_t *_ref_count

    def __cinit__(self):
        self._seq_list = NULL
        self._n_seqs = 0
        self._is_paired_end = False

        self._ref_count = <uint32_t*>malloc(sizeof(uint32_t))
        self._ref_count[0] = 1

    def __dealloc__(self):
        self._ref_count[0] -= 1
        if self._ref_count[0] == 0:
            if self._ref_count is not NULL:
                free(self._ref_count)
            destroy_seqs(self._seq_list, self._n_seqs)
            
    cdef cbwa.bseq1_t *raw(self):
        return self._seq_list

    @staticmethod
    def fromfiles(const char *reads1, const char *reads2 = NULL,
                  uint64_t chunk_size = 10000000):
        """Create a sequence list from files.

        Read some fastq files and populate a sequence list with the content
        from these files. If `reads2` is set, we assume that the files contain
        paired end reads.

        Parameters
        ----------
        reads1, reads2 : bytes
                         (const char*) str filepaths to the input files. Must
                         be gzipped.
        chunk_size : uint64_t
                     Number of bases to read per iteration.

        Returns
        -------
        SequenceList
            SequenceList object containing unaligned sequences from the fastq
            files.

        Raises
        ------
        ValueError
            Raises a value error if any of the files passed could not be opened
            for reading.
        """
        cdef SequenceList sl = SequenceList()

        cdef gzFile _fp1
        cdef gzFile _fp2
        cdef ckseq.kseq_t *_ks1 = NULL
        cdef ckseq.kseq_t *_ks2 = NULL
        cdef int n_seqs
        cdef cbwa.bseq1_t *seqs
        cdef int n_processed = 0
        cdef cbwamem.mem_pestat_t pes[4]

        try:
            memset(pes, 0, 4 * sizeof(cbwamem.mem_pestat_t))

            _fp1 = gzopen(reads1, "r")
            if _fp1 is NULL:
                raise ValueError(f"Could not open {reads1}.")
            _ks1 = ckseq.kseq_init(_fp1)

            if reads2 is not NULL:
                _fp2 = gzopen(reads2, "r")
                if _fp2 is NULL:
                    raise ValueError(f"Could not open {reads2}.")
                _ks2 = ckseq.kseq_init(_fp2)
                sl._is_paired_end = True

            seqs = cbwa.bseq_read(chunk_size, &n_seqs, _ks1, _ks2)
            if n_seqs > 0:
                sl.extend(seqs, n_seqs)
            print(f"Read {n_seqs} sequences.")

            while n_seqs > 0:
                seqs = cbwa.bseq_read(chunk_size, &n_seqs, _ks1, _ks2)
                if n_seqs > 0:
                    sl.extend(seqs, n_seqs)
                    print(f"Read {n_seqs} sequences.")
        finally:
            if _ks1 is not NULL:
                ckseq.kseq_destroy(_ks1)
            if _fp1 is not NULL:
                gzclose(_fp1)
            if reads2 is not NULL:
                if _ks2 is not NULL:
                    ckseq.kseq_destroy(_ks2)
                if _fp2 is not NULL:
                    gzclose(_fp2)
        return sl

    @staticmethod
    def frombytes(reads1: typing.Iterable[bytes], reads2: typing.Iterable[bytes] = None,
                  uint64_t pre_alloc = 1, uint64_t max_len = 0xFFFFFFFFFFFFFFFF):
        """Create a sequence list from string iterables.

        Read the string iterables and parse them into a SequenceList
        representation that can be used for alignment. If `reads2` is passed,
        we assume that the iterables contain paired end reads.

        Parameters
        ----------
        reads1, reads2 : Iterable[bytes]
                         Some string iterables. Length must be divisble by 4
                         (see fastq specification). If both are passed, the
                         lengths must be equal. If you pass strings from
                         python, make sure to encode them to bytes objects.
        pre_alloc : uint64_t
                    Number of reads to allocate space for on the first pass.
                    If you know the number of reads passed, you can pass it
                    here to save on some memory (re-)allocation.
        max_len : uint64_t
                  Maximum number of sequences to read.
        
        Returns
        -------
        SequenceList
            A sequence list with the data parsed from both input streams.

        Raises
        ------
        ValueError
            The input iterables do not match the specification.
        """
        cdef SequenceList sl
        cdef uint64_t len1 = 0, len2 = 0, l_s = 0, idx = 0
        cdef cbwa.bseq1_t *seqs = \
            <cbwa.bseq1_t*>malloc(pre_alloc * sizeof(cbwa.bseq1_t))
        l_s = pre_alloc

        for r in reads1:
            if r[0] == ord(b'@'):
                # Write previous sequence
                if len1 > 1:
                    if idx + 1 >= l_s:
                        l_s = l_s << 1
                        seqs = <cbwa.bseq1_t*>realloc(seqs, l_s * sizeof(cbwa.bseq1_t))

                    _BseqRead(
                        name1,
                        comment1,
                        seq1,
                        qual1
                    ).to_bseq1(&seqs[idx])
                    seqs[idx].id = idx
                    idx += 1
                    if reads2 is not None:
                        _BseqRead(
                            name2,
                            comment2,
                            seq2,
                            qual2
                        ).to_bseq1(&seqs[idx])
                        seqs[idx].id = idx
                        idx += 1

            read1 = r
            len1 += 1
            if reads2 is not None:
                try:
                    read2 = next(reads2)
                    len2 += 1
                except StopIteration:
                    raise ValueError("Read iterables have different lengths.")
            else:
                read2 = None

            if idx >= max_len:
                break

            if len1 % 4 == 0:
                # Quality Scores
                qual1 = read1
                if reads2 is not None:
                    qual2 = read2
            elif len1 % 4 == 1:
                # Sequence ID
                name1, comment1 = read1.split(b" ")
                name1 = name1[1:]
                if read2 is not None:
                    name2, comment2 = read2.split(b" ")
                    name2 = name2[1:]
            elif len1 % 4 == 2:
                # Sequence
                seq1 = read1
                if reads2 is not None:
                    seq2 = read2
            elif len1 % 4 == 3:
                # +
                pass

        if len1 != len2:
            free(seqs)
            raise ValueError("Read iterables have different lengths.")
        if len1 % 4 != 0 or len2 % 4 != 0:
            free(seqs)
            raise ValueError("Malformed input.")

        sl = SequenceList()
        sl._n_seqs = idx
        sl._seq_list = seqs
        sl._is_paired_end = reads2 is not None
        return sl
    
    @staticmethod
    def frombseqs(sequences: typing.Iterable[Bseq], bint is_paired_end = False,
                  uint64_t pre_alloc = 1, uint64_t max_len = 0xFFFFFFFFFFFFFFFF):
        """Create a sequence list from bseq iterables.
        
        Create a new SequenceList object based on the sequences given in
        the bseq iterable. The content of the iterable _must_ be Bseq objects
        created by this module. Each bseq object is copied and the sam
        field is discarded.
        
        If the data contains paired end reads, they are expected to be
        interleaved, s.t. all even indices and zero contain forward reads
        and all odd indices contain reverse reads. This is also how BWA
        handles them. Additionally, one must pass the flag is_paired_end to
        clearly signal the intent.
        
        Parameters
        ----------
        sequences : Iterable[Bseq]
                    Any python iterable (i.e. objects with the __next__ method)
                    returning Bseq objects created by this lib.
        is_paired_end : bool
                        Sets the is_paired_end flag. Default: False
        pre_alloc : uint64_t
                    Number of sequences objects to pre allocate. Improves
                    performance if the number of sequences is known in advance.
        max_len : uint64_t
                  Number of sequences to copy. If set, only the first n
                  sequences will be copied.
        
        Returns
        -------
        SequenceList
            A new SequenceList object containing data from the input iterable.
            
        Raises
        ------
        ValueError
            The iterable did not return proper Bseq objects.
        """
        cdef SequenceList sl
        cdef uint64_t l_s = 0, idx = 0
        cdef cbwa.bseq1_t *seqs = <cbwa.bseq1_t*>malloc(pre_alloc * sizeof(cbwa.bseq1_t))
        l_s = pre_alloc
        
        for b in sequences:
            if isinstance(b, Bseq):
                if idx >= l_s:
                    l_s = l_s << 1
                    seqs = <cbwa.bseq1_t*>realloc(seqs, l_s*sizeof(cbwa.bseq1_t))
                # (Deep)-Copy the original sequence object
                seqs[idx].l_seq = b.l_seq
                seqs[idx].id = idx
                seqs[idx].name = copy_string(
                    bytes(b.name.encode()),
                    len(bytes(b.name.encode())),
                    True
                )
                seqs[idx].comment = copy_string(
                    bytes(b.comment.encode()),
                    len(bytes(b.comment.encode())),
                    False
                )
                seqs[idx].seq = copy_string(
                    b.seq.encode(),
                    len(b.seq.encode()),
                    True
                )
                seqs[idx].qual = copy_string(
                    bytes(b.qual.encode()),
                    len(bytes(b.qual.encode())),
                    False
                )
                seqs[idx].sam = copy_string("", 0, False)
                idx += 1
            else:
                raise ValueError("Must provide Bseq iterable")
        
        sl = SequenceList()
        sl._n_seqs = idx
        sl._is_paired_end = is_paired_end
        sl._seq_list = seqs
        return sl

    @staticmethod
    def from_alignedsegments(sequences: typing.Iterable[AlignedSegment], bint is_paired_end = False,
                             uint64_t pre_alloc = 1, uint64_t max_len = 0xFFFFFFFFFFFFFFFF):
        """Create a SequenceList from an iterable of pysam.AlignedSegment objects.

        Create a new SequenceList object by copying data from the pysam AlignedSegment
        objects passed to this function. For each aligned segment, we copy all contents
        of the object, namely the `query_name`, `query_sequence` and `query_qualities`
        properties to a Bseq object. Alignment information is dropped when converting
        to the Bseq object.

        Warning
        -------
        When working with paired end data (e.g. `is_paired_end = True`), the caller is
        responsible for passing in interleaved reads. This method will not sort the input
        iterable by `query_name`, or make any other effort to put read pairs next to
        each other, as it is expected by the `Bwa.align` method.

        Parameters
        ----------
        sequences: Iterable[AlignedSegment]
                   An iterable of pysam.AlignedSegment objects.
        is_paired_end: bool
                       Set the `is_paired_end` flag on the resulting SequenceList object.
                       If true, make sure the input iterable provides fastq pairs in order.
        pre_alloc: uint64_t
                   Size of the array to allocate on first attempt. If the number of sequences
                   is known in advance, this can be used to save some allocations.
        max_len: uint64_t
                 Number of elements to copy from the input iterable to the SequenceList.

        Returns
        -------
        SequenceList
            A new SequenceList object containing data from the input iterable.

        Raises
        ------
        ValueError
            The input iterable contains other objects than `pysam.AlignedSegment` instances.
        MemoryError
            The runtime could not allocate enough memory to store the Bseq objects.
        """
        cdef SequenceList sl
        cdef uint64_t l_s = 0, idx = 0
        cdef cbwa.bseq1_t *seqs = <cbwa.bseq1_t*>malloc(pre_alloc * sizeof(cbwa.bseq1_t))
        if seqs is NULL:
            raise MemoryError()

        l_s = pre_alloc

        for read in sequences:
            if isinstance(read, AlignedSegment):
                if idx >= l_s:
                    l_s = l_s << 1
                    seqs = <cbwa.bseq1_t*>realloc(seqs, l_s * sizeof(cbwa.bseq1_t))
                    if seqs is NULL:
                        raise MemoryError()

                if idx >= max_len:
                    break

                if read.query_length == 0:
                    # Ignore empty reads
                    continue

                seqs[idx].l_seq = read.query_length
                seqs[idx].id = idx
                seqs[idx].name = copy_string(
                    bytes(read.query_name.encode()),
                    len(bytes(read.query_name.encode())),
                    False,
                )
                ##
                # Here, we utilize the fact that AlignedSegment.get_forward_qualities
                # from pysam is idempotent w.r.t. the original orientation of the read.
                # If the read was not reverse complemented (i.e. flag 0x10 is not set),
                # the function will return the same read / qual sequence. If the read
                # was reverse complemented (flag & 0x10 != 0) then get_forward_sequence
                # will return the reverse complement sequence and get_forward_qualities
                # will return the reversed qual array.
                seqs[idx].comment = copy_string("", 0, False)
                seqs[idx].seq = seq_to_two_bit(
                    bytes(read.get_forward_sequence().encode()),
                    read.query_length,
                )
                seqs[idx].qual = get_qual_string(
                    read.get_forward_qualities(),
                    read.query_length
                )
                seqs[idx].sam = copy_string("", 0, False)
                idx += 1
            else:
                raise ValueError("Must provide a AlignedSegment iterable")
        
        sl = SequenceList()
        sl._n_seqs = idx
        sl._is_paired_end = is_paired_end
        sl._seq_list = seqs
        return sl

    cdef void extend(self, cbwa.bseq1_t *seqs, uint32_t n_seq):
        cdef uint64_t prev_size = self._n_seqs

        if self._seq_list is NULL:
            self._seq_list = <cbwa.bseq1_t*>malloc(n_seq * sizeof(cbwa.bseq1_t))
            if self._seq_list is NULL:
                raise MemoryError()
            memcpy(self._seq_list, seqs, n_seq * sizeof(cbwa.bseq1_t))
            self._n_seqs = n_seq
        else:
            self._n_seqs += n_seq
            new_ptr =  <cbwa.bseq1_t*>realloc(self._seq_list, self._n_seqs * sizeof(cbwa.bseq1_t))
            if new_ptr is NULL:
                raise MemoryError()
            self._seq_list = new_ptr
            memcpy(&self._seq_list[prev_size], seqs, n_seq * sizeof(cbwa.bseq1_t))

    @property
    def is_paired_end(self):
        """bool: List contains paired end reads, if true.
        
        Denotes whether the list contains paired end reads or not.
        If true, paired end reads are interleaved (i.e. even
        indices contain forward reads, odd ones reverse reads).
        """
        return self._is_paired_end

    def __iter__(self) -> SequenceListForwardIterator:
        """Create an iterator for this list.
        
        Returns
        -------
        SequenceListForwardIterator
            An iterator
        """
        cdef SequenceListForwardIterator it = SequenceListForwardIterator()
        it._seq_list = self._seq_list
        it._n_seqs = self._n_seqs
        it._ref_count = self._ref_count
        self._ref_count[0] += 1
        return it

    def __len__(self) -> int:
        return self._n_seqs

    def __getitem__(self, acc):
        """Get a member or part of the list.
        
        Parameters
        ----------
        acc : Union[int, slice]
              An integer or slice object. Int gives the 0-based index of the
              sequence to be returned. Slice objects will return a new View on
              the data for the defined interval.
        
        Returns
        -------
        Union[Bseq, SequenceListView]
            If an integer was used to access data it returns the Bseq object
            for the corresponding position in the list. Otherwise a view on
            the interval defined by the slice object is returned.
            
        Raises
        ------
        TypeError
            Raises a TypeError if something else than an integer or slice object
            is passed to the function.
        IndexError
            The integer index is out of range.
        """
        cdef SequenceListView slv
        cdef uint64_t start = 0, stop = 0, step = 0
        cdef uint32_t i = 0

        if isinstance(acc, int):
            if acc < 0:
                i = self._n_seqs + acc
            else:
                i = acc
            if i < 0 or i >= self._n_seqs:
                raise IndexError(f"Index out of range [0,{self._total_seqs})")
            else:
                return Bseq.wrap(&self._seq_list[i])
        elif isinstance(acc, slice):
            if acc.start is not None and acc.start >= 0:
                start = acc.start
            elif acc.start is not None and acc.start < 0:
                start = self._n_seqs + acc.start
            else:
                start = 0

            if acc.stop is None or acc.stop > self._n_seqs:
                stop = self._n_seqs
            elif acc.stop < 0:
                stop = self._n_seqs + acc.stop
            else:
                stop = acc.stop

            if acc.step is not None and acc.step > 0:
                step = acc.step
            else:
                step = 1

            # Allocate SequenceListView
            slv = SequenceListView()
            slv._seq_list = self._seq_list
            slv._ref_count = self._ref_count
            self._ref_count[0] += 1
            slv._start = start
            slv._stop = stop
            slv._step = step
            slv._n_seqs = (stop - start) // step
            slv._is_paired_end = self._is_paired_end
            return slv
        else:
            raise TypeError("Type not supported by getitem.")

    def __str__(self) -> str:
        return f"{self.__class__.__name__}=(" \
            + f"n_seqs={self._n_seqs},is_paired_end={self._is_paired_end})"


cdef class Bwa:
    """BWA Alignment wrapper.
    
    Thin wrapper around some of the BWA API's. Currently, only MEM alignment
    and building/loading indices is supported. By default, no options are set.
    You can use the options members to fine-tune the BWA algorithm settings.
    
    Parameters
    ----------
    index : bytes
            (const char*) Path to a BWA index prefix.
    bwtalgo : BWA_IDX
              Algorithm to load the BWA index with. By default, all bwa index
              versions are used.
    """
    cdef BwamemOptions _opt
    cdef cbwa.bwaidx_t *_idx

    def __cinit__(
        self,
        const char *index,
        int bwtalgo = BWA_IDX.I_ALL,
    ):
        self._opt = BwamemOptions()
        self._idx = cbwa.bwa_idx_load(index, bwtalgo)

    def align(self, SequenceList sl, uint64_t chunk_size = 1000000):
        """Align a set of (paired end) fastq's against a reference.

        Given the location of two fastq files, we return an alignment of all
        sequences against the reference genome the BWA instance was constructed
        with.

        Parameters
        ----------
        sl : SequenceList
             A sequence list object containing all reads that should be aligned.
        chunk_size : uint64_t, optional
                     Number of reads per iteration of the read-align loop.

        Returns
        -------
        SequenceList
            The input sequence list object with the `sam` member set for all
            `bseq1_t` structs that are part of the list.
        """
        cdef uint64_t n_processed = 0
        cdef cbwamem.mem_pestat_t pes[4]
        memset(pes, 0, 4 * sizeof(cbwamem.mem_pestat_t))
        
        if sl.is_paired_end and not self._opt.has_flag(MEM_F.PE):
            self._opt.set_flag(MEM_F.PE)

        cbwamem.mem_process_seqs(
            self._opt.raw(),
            self._idx.bwt,
            self._idx.bns,
            self._idx.pac,
            n_processed,
            sl._n_seqs,
            sl.raw(),
            pes
        )
        return sl

    def __dealloc__(self):
        cbwa.bwa_idx_destroy(self._idx)

    @property
    def options(self) -> BwamemOptions:
        """BwamemOptions: Read/Write mem algorithm options."""
        return self._opt

    @staticmethod
    def build_index(const char *fa, const char *prefix, int algo_type = BWTALGO.AUTO,
                    int block_size = 10000000):
        """Build a bwa index for a fasta file.

        Build the BWT index for the given fastq file and save the files to the
        location at prefix.

        The prefix can then be used to create a bwa object using the default
        constructor.

        Parameters
        ----------
        fa : bytes
             (const char*) Location of the input fasta file.
        prefix : bytes
                 (const char*) Prefix of the output files.
        algo_type : BWTALGO
                    BWTALGO member. Algorithm used to build the BWT index.
        block_size : int
                     BWT block size parameter.
        """
        cdef bint ret = False
        ret = cbwa.bwa_idx_build(fa, prefix, algo_type, block_size)
        return ret
    
    def get_sam_header(self) -> str:
        """Re-Implementation of `bwa_print_sam_hdr`.

        Print the SAM header lines for the current index.
        
        Returns
        -------
        str
            The SAM header for the currently loaded reference index. An additional
            program line is added.

        Notes
        -----
        The bwa_print_sam_hdr function only prints the header lines to stdout,
        and not to a char*. Therefore we reimplement it.
        """
        cdef cbntseq.bntseq_t *bns = self._idx.bns
        cdef uint32_t i
        header = str()

        for i in range(bns.n_seqs):
            header += f"@SQ\tSN:{(<bytes>bns.anns[i].name).decode('ascii')}\tLN:{bns.anns[i].len}"
            if bns.anns[i].is_alt:
                header += "\tAH:*\n"
            else:
                header += "\n"

        header += __bwa_pg__ + "\n"
        return header