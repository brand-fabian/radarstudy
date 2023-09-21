from multiprocessing.sharedctypes import Value
import bwa
import pysam
import os
from meta import ReferenceMatcher
import pandas
import typing
import logging
from itertools import chain


logger = logging.getLogger("subset-bam")


class GetReadsOptions(typing.NamedTuple):
    flanking_region: int = 5000
    target_reference_name: str = "grch37"
    locus_col: str = "locus"
    locus_col_remap: str = "grch38.locus"
    threads: int = 4
    fetch_unmapped_reads: bool = False
    ignore_unmatched_reads: bool = False
DefaultGetReadsOptions = GetReadsOptions()


class Locus(typing.NamedTuple):
    contig: str
    pos: int

    @classmethod
    def fromstring(cls, x: str):
        if x is None:
            return None
        contig, pos = x.split(":")
        return cls(
            contig=contig,
            pos=int(pos),
        )


def _get_interleaved_reads(
    bam: pysam.AlignmentFile,
    contig: str,
    start: int,
    end: int,
    fetch_unmapped_reads = False,
    ignore_unmatched_reads = False,
):
    """Get interleaved read pairs from a bam file region.
    
    Parameters
    ----------
    bam : pysam.AlignmentFile
          PySAM BamReader instance.
    contig : str
             Contig to fetch from the bam/sam reader.
    start : int
            Start of the chromosomal region to fetch.
    end : int
          End of the chromosomal region to fetch.
    fetch_unmapped_reads : boolean
                           (Default = False) If true, attempt to fetch
                           the mate reads from the bam file for all singular
                           reads after finishing reading the region. Very slow.
                           If `ignore_unmatched_reads` is False, then the singleton
                           reads will still be emitted by the generator.
    ignore_unmatched_reads : boolean
                             (Default = False) Do not emit singleton reads,
                             and do not fetch their mate reads.
    
    Returns
    -------
    Generator[pysam.AlignedSegment]
        Iterator of interleaved read pairs.
    """
    read_pairs = {}
    matched_pairs = 0
    try:
        for read in bam.fetch(contig, start, end):
            if read.query_name not in read_pairs:
                read_pairs[read.query_name] = []
            read_pairs[read.query_name].append(read)
            if len(read_pairs[read.query_name]) == 2:
                yield read_pairs[read.query_name][0]
                yield read_pairs[read.query_name][1]
                del read_pairs[read.query_name]
                matched_pairs += 1
        logger.debug("got {} matched and {} unmatched read pairs".format(
            matched_pairs, len(read_pairs),
        ))
        if not ignore_unmatched_reads:
            for pair in read_pairs.values():
                yield pair[0]
                if fetch_unmapped_reads:
                    yield bam.mate(pair[0])
    except ValueError as e:
        logger.error("could not read {}:{}-{}: {}".format(
            contig, start, end, str(e)
        ))


def _get_alignedsegment(bseq, header):
    # Read supplementary alignments
    if '\n' in bseq.sam[:-1]:
        return [
            pysam.AlignedSegment.fromstring(s, header)
            for s in bseq.sam[:-1].split('\n')
        ]
    else:
        return [
            pysam.AlignedSegment.fromstring(bseq.sam, header)
        ]


class Mapper:
    def __init__(self, bam: str, reference: str, threads: int):
        self._bwa = bwa.Bwa(reference.encode())
        self._bwa.options.n_threads = threads
        self._bwa.options.set_flag(bwa.MEM_F.NO_MULTI)
        self._bam = os.path.abspath(bam)

    def remap(self, contig: str, start: int, end: int, header = None,
              fetch_unmapped_reads: bool = False, ignore_umatched_reads: bool = False):
        if header is None:
            header = self.header
        bam = pysam.AlignmentFile(self._bam)
        try:
            sequence_list = bwa.SequenceList.from_alignedsegments(
                _get_interleaved_reads(bam, contig, start, end,
                                       fetch_unmapped_reads=fetch_unmapped_reads,
                                       ignore_unmatched_reads=ignore_umatched_reads),
                is_paired_end=(fetch_unmapped_reads or ignore_umatched_reads)
            )
            self._bwa.align(sequence_list)
            return list(chain.from_iterable(map(
                lambda bseq: _get_alignedsegment(bseq, header),
                sequence_list,
            )))
        finally:
            bam.close()

    @property
    def header(self):
        return pysam.AlignmentHeader.from_text(self._bwa.get_sam_header())


def get_reads(
    sites: pandas.DataFrame,
    bam: str,
    reference: str,
    reference_matcher: ReferenceMatcher,
    options: GetReadsOptions = DefaultGetReadsOptions,
) -> typing.Tuple[pysam.AlignmentHeader, typing.List[pysam.AlignedSegment]]:
    """Read sites from the underlying bam file accounting for references.
    
    Fetch reads from a bam file and optionally remap them to a different
    genome build.

    Parameters
    ----------
    sites : pandas.DataFrame
            Sites table with `options.locus_col` or `options.locus_col_remap`
            columns.
    bam : str
          Path to a bam file
    reference : str
                Path to BWA-indexed .fasta reference file for the genome build
                to remap to, if there is a mismatch.
    reference_matcher : ReferenceMatcher
                        ReferenceMatcher object to get the genome reference from
                        the input bam file.
    options : GetReadsOptions
              Generic Options.
    
    Returns
    -------
    pysam.AlignmentHeader
        Header of the returned AlignedSegmment objects.
    List[pysam.AlignedSegment]
        List of aligned segment mapped to the original genome build (`reference`).
    """
    reference_genome = reference_matcher.get(bam)
    logger.debug("detected reference {} for file {}".format(reference_genome, bam))
    if reference_genome != options.target_reference_name:
        logger.info("applying remapping to {} for file {}".format(options.target_reference_name, bam))
        loci = list(filter(
            lambda x: x is not None,
            map(Locus.fromstring, sites[options.locus_col_remap]))
        )
        bwa = Mapper(bam, reference, options.threads)
        logger.info("remapping {} loci".format(len(loci)))
        reads = []
        for l in loci:
            reads.extend(bwa.remap(
                l.contig,
                max(0, l.pos - options.flanking_region),
                l.pos + options.flanking_region,
                fetch_unmapped_reads=options.fetch_unmapped_reads,
                ignore_umatched_reads=options.ignore_unmatched_reads,
            ))
        return bwa.header, reads
    else:
        logger.info("using bam {} as is".format(bam))
        loci = list(filter(
            lambda x: x is not None,
            map(Locus.fromstring, sites[options.locus_col]))
        )
        reads = []
        bam_f = pysam.AlignmentFile(os.path.abspath(bam))
        try:
            for l in loci:
                reads.extend(bam_f.fetch(
                    l.contig,
                    max(0, l.pos - options.flanking_region),
                    l.pos + options.flanking_region,
                ))
            return bam_f.header, reads
        finally:
            bam_f.close()


def fix_readgroup(
    header: pysam.AlignmentHeader,
    reads: typing.List[pysam.AlignedSegment],
    sample: str,
    query_name_start: int = 0,
    overwrite_query_name: bool = False,
) -> typing.Tuple[pysam.AlignmentHeader, typing.List[pysam.AlignedSegment]]:
    """Fix read group annotations in a set of reads.

    For the given list of reads, fix the read group annotations, by
    overwriting it with the given sample id.

    Parameters
    ----------
    header : pysam.AlignmentHeader
             Header to add/change 'RG' tags in
    reads : List[pysam.AlignedSegment]
            List of reads that will be annotated with the given read group.
    sample : str
             Sample ID to use for the read groups.
    
    Returns
    -------
    pysam.AlignmentHeader
        The new alignment file header.
    List[pysam.AlignedSegment]
        The list of changed read objects.
    """
    header_d = header.to_dict()
    if 'RG' in header_d:
        header_d['RG']['SM'] = sample
        header_d['RG']['ID'] = sample
    else:
        header_d['RG'] = [{
            'ID': sample,
            'PL': 'Illumina',
            'SM': sample,
            'LB': sample,
        }]

    query_name = query_name_start
    for read in reads:
        read.set_tag('RG', sample, value_type='Z', replace=True)
        if overwrite_query_name:
            read.query_name = str(query_name)
        query_name += 1

    return (
        pysam.AlignmentHeader.from_dict(header_d),
        reads,
    )


def merge_readgroup(*headers: typing.List[pysam.AlignmentHeader]):
    """Merge read groups from a set of headers."""
    rgs = chain.from_iterable([ h.to_dict()['RG'] for h in headers ])
    rg_samples = set()
    unique_rgs = []
    for rg in rgs:
        if rg['ID'] not in rg_samples:
            rg_samples.add(rg['ID'])
            unique_rgs.append(rg)

    hdr = headers[0].to_dict()
    hdr['RG'] = unique_rgs
    return pysam.AlignmentHeader.from_dict(hdr)