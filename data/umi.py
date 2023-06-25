from Bio import SeqIO
import sys, gzip

with gzip.open(sys.argv[1],"rt") as t:
    with gzip.open(sys.argv[2], "wt") as tt:
        for record in SeqIO.parse(t,"fastq"):
            a=record.id.split("_")
            record.description=""
            record.id=a[0]+":UMI_"+a[1]
            SeqIO.write(record, tt, "fastq")
