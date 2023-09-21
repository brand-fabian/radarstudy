import bwa

b = bwa.Bwa(b"data/thy1.fa")
sl = bwa.SequenceList.fromfiles(
    b"data/R1.fastq.gz",
    b"data/R2.fastq.gz"
)
aln = b.align(sl)
for i, x in enumerate(sl):
    print(i, x)
    if i >= 10:
        break