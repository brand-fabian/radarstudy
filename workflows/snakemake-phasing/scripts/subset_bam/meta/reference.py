import pysam
import typing
import os
import logging

logger = logging.getLogger(__name__)

class ReferenceMatcher:
    """Get reference from bam files"""
    class _Reference:
        def __init__(self, fasta: str):
            self._fasta = os.path.abspath(fasta)
            self._contigs = {}
        
        @property
        def contigs(self) -> typing.Dict[str, int]:
            if len(self._contigs) == 0:
                fasta = pysam.FastaFile(self._fasta)
                try:
                    self._contigs = {
                        r: fasta.get_reference_length(r)
                        for r in fasta.references
                    }
                finally:
                    fasta.close()
            return self._contigs

    def __init__(self):
        self._references: typing.Dict[str, ReferenceMatcher._Reference] = {}

    def add(self, name: str, fasta: str):
        self._references[name] = ReferenceMatcher._Reference(fasta)

    def get(self, bam: str):
        bam_f = pysam.AlignmentFile(os.path.abspath(bam))
        try:
            bam_references: typing.Dict[str, int] = {
                r[0]: r[1]
                for r in zip(bam_f.references, bam_f.header.lengths)
            }
            for key, ref in self._references.items():
                if ref.contigs == bam_references:
                    return key
        finally:
            bam_f.close()
        logger.warning("could not assert reference for {}".format(bam))
        return None