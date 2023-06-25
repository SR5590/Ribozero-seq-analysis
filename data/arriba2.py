import glob
import sys
import subprocess
id=sys.argv[1]
a=glob.glob("../04_Arriba/"+id+ "/*Fusion*")
bam="../03_star/"+id+"/UmiReads/"+id+"Aligned.sortedByCoord.out.bam"
for file in a:
    subprocess.call(["draw_fusions.R",
        "--annotation=/mnt/g27prist/CMTD/Stephan/bcbio_installation/genomes/Hsapiens/GRCh37/rnaseq/ref-transcripts.gtf",
        "--fusions="+file, "--output=" +file+".pdf", "--alignments="+bam,
        "--cytobands=/mnt/g27prist/CMTD/Stephan/bcbio_installation/genomes/Hsapiens/GRCh37/rnaseq/fusion-blacklist/arriba-cytobands.tsv",
        "--minConfidenceForCircosPlot=high",
        "--proteinDomains=/mnt/g27prist/CMTD/Stephan/bcbio_installation/genomes/Hsapiens/GRCh37/rnaseq/fusion-blacklist/arriba-protein-domains.gff3"])
    subprocess.call(["pdftoppm", file+".pdf" ,file,"-png"])
    #subprocess.call(["touch", ""])

subprocess.call(["touch","../04_Arriba/"+id+"/"+id+".done2"])
