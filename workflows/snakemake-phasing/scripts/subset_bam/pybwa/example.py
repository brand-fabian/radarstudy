import bwa

b = bwa.Bwa(b"/mnt/e/Projects/pybwa/data/mm10.fa")
for x in b.align(
    b"/mnt/e/Projects/pybwa/data/GFP_S28_R1_001.fastq.gz",
    b"/mnt/e/Projects/pybwa/data/GFP_S28_R2_001.fastq.gz"
):
    print(str(x[1]), str(x[0]))