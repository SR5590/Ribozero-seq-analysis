from Bio import SeqIO
import sys, gzip
from Bio.Seq import Seq
with gzip.open(sys.argv[1],"rt") as t:
    with gzip.open(sys.argv[2], "wt") as tt:
        with gzip.open(sys.argv[3], "wt") as ttt:
            for record in SeqIO.parse(t,"fastq"):
                umi=record.id.split(":")[-1].replace('+','')
                record.description=""
                record.id=':'.join(record.id.split(':')[:-1])
                SeqIO.write(record, tt, "fastq")
                record.letter_annotations={}
                record.seq=Seq(umi)
                record.letter_annotations["phred_quality"] = [40] * len(record)
                SeqIO.write(record,ttt,"fastq")
