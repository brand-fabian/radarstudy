import bwa
import pysam
import sys

if len(sys.argv) < 2:
    print("multi.py <path to mm10>")
    sys.exit(1)

mm10 = sys.argv[1]

b = bwa.Bwa(b"data/thy1.fa")
sl = bwa.SequenceList.fromfiles(
    b"data/R1.fastq.gz",
    b"data/R2.fastq.gz"
)
v_aln = b.align(sl)
bam = pysam.AlignmentFile("test.bam", "wb", text=b.get_sam_header())

vec_seqs = []
i =  0
for read, mate in zip(sl[::2], sl[1::2]):
    i += 1
    read_s = pysam.AlignedSegment.fromstring(read.sam, bam.header)
    mate_s = pysam.AlignedSegment.fromstring(mate.sam, bam.header)
    if not read_s.is_unmapped or not mate_s.is_unmapped:
        if read.name != mate.name:
            print("Name mismatch:", read.name, mate.name)
        vec_seqs.extend([read, mate])
    if i % 10000 == 0:
        print(i)

sl2 = bwa.SequenceList.frombseqs(vec_seqs, is_paired_end=True)
print(len(sl2))
mouse = bwa.Bwa(mm10.encode())
m_aln = mouse.align(sl2)

for i, x in enumerate(m_aln):
    print(i, x)
    if i > 10:
        break