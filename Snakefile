import os
import pandas as pd
from collections import Counter
import glob
SAMPLES = []
DIRECTION=["1","2"]
GROUPS=set()
sample_sheet=pd.read_csv("samples.tsv", sep="\t")
count={}
syn=[]

sample_sheet["working"]=sample_sheet["ID"].astype(str)+sample_sheet["Type"]

def check_symlink(file1, file2):
    try:
        os.symlink(file1, file2)
    except FileExistsError:
        print("Link for " + file1 + " is already present in 01_raw")

WORKING_SAMPLES=list(sample_sheet["working"])

sample_sheet.to_csv("sample_working.tsv",index=False, sep="\t")


sample_sheet
for b in set(WORKING_SAMPLES):
    os.makedirs("../01_raw/" +b+ "/fastqc", exist_ok=True)
    check_symlink(sample_sheet[sample_sheet["working"]==b]["forward reads"].values[0], "../01_raw/"+b +"/"+ b +"_1P.fastq.gz")
    check_symlink(sample_sheet[sample_sheet["working"]==b]["reverse reads"].values[0], "../01_raw/"+b +"/"+ b +"_2P.fastq.gz")



rule all:
    input:
        expand('../01_raw/{sample}/fastqc/{sample}_{direction}P_fastqc.html',direction=DIRECTION, sample=WORKING_SAMPLES),
#        expand("../03_umi/{sample}/{sample}_umi_{direction}P.fastq.gz",sample=WORKING_SAMPLES, direction=DIRECTION),
#        expand("../02_trim/{sample}/{sample}_umi_trimmed_{direction}P.fastq.gz",sample=WORKING_SAMPLES, direction=DIRECTION),
#        expand("../03_star/{sample}/{sample}Aligned.sortedByCoord.out.bam",sample=WORKING_SAMPLES, direction=DIRECTION),
#        expand("../03_star/{sample}/{sample}Aligned.sortedByCoord.out.bam.bai",sample=WORKING_SAMPLES, direction=DIRECTION),
        expand("../04_Arriba/{sample}/{sample}.fusion.tsv",sample=WORKING_SAMPLES, direction=DIRECTION),
        expand("../04_Arriba/{sample}/{sample}.done",sample=WORKING_SAMPLES, direction=DIRECTION),
        expand("../04_Arriba/{sample}/{sample}.done2",sample=WORKING_SAMPLES, direction=DIRECTION),
        expand("../05_FeatureCount/{sample}_FeatureCountTable",sample=WORKING_SAMPLES),
        expand('../06_quality/{sample}.pdf',sample=WORKING_SAMPLES),
        expand("../03_star/{sample}/{sample}_CallCons.bam",sample=WORKING_SAMPLES),
        expand("../03_star/{sample}/{sample}_GroupRead.bam",sample=WORKING_SAMPLES)






rule fastqc1:
    input:
        r = '../01_raw/{sample}/{sample}_{direction}P.fastq.gz'
    threads: 1
    priority: 50
    output: '../01_raw/{sample}/fastqc/{sample}_{direction}P_fastqc.html'
    conda: "envs/fastqc.yaml"
    shell: 'fastqc -o ../01_raw/{wildcards.sample}/fastqc -t {threads} --extract {input.r}'




rule trimming:
    input:
        umi1="../01_raw/{sample}/{sample}_1P.fastq.gz",
        umi2="../01_raw/{sample}/{sample}_2P.fastq.gz"

    output:
        p1="../02_trim/{sample}/{sample}_umi_trimmed_1P.fastq.gz",
        u1="../02_trim/{sample}/{sample}_umi_trimmed_1U.fastq.gz",
        p2="../02_trim/{sample}/{sample}_umi_trimmed_2P.fastq.gz",
        u2="../02_trim/{sample}/{sample}_umi_trimmed_2U.fastq.gz",
    conda:
        "envs/trimmomatic.yaml"
    threads: 8
    priority: 50
    shell:
        'trimmomatic PE -threads {threads} \
        {input.umi1} {input.umi2} \
        {output.p1} {output.u1} \
        {output.p2} {output.u2} \
        ILLUMINACLIP:data/adapters.fa:2:30:10:2:true \
        MINLEN:36'

rule umi_prep1:
    input:
        p1="../02_trim/{sample}/{sample}_umi_trimmed_1P.fastq.gz",
    output:
        r1="../02_trim/{sample}/{sample}_umiClean_trimmed_1P.fastq.gz",
        r1_umi="../02_trim/{sample}/{sample}_umis_1P.fastq.gz",
    conda: "envs/umi2.yaml",
    threads: 1
    priority: 50
    shell:"python data/umi_tag.py {input.p1} {output.r1} {output.r1_umi}"


rule umi_prep2:
    input:
        p2="../02_trim/{sample}/{sample}_umi_trimmed_2P.fastq.gz",
    output:
        r2="../02_trim/{sample}/{sample}_umiClean_trimmed_2P.fastq.gz",
        r2_umi="../02_trim/{sample}/{sample}_umis_2P.fastq.gz",
    conda: "envs/umi2.yaml",
    threads: 1
    priority: 50
    shell:"python data/umi_tag.py {input.p2} {output.r2} {output.r2_umi}"



rule StarMapping:
    input:
        r1 = "../02_trim/{sample}/{sample}_umiClean_trimmed_1P.fastq.gz",
        r2 = "../02_trim/{sample}/{sample}_umiClean_trimmed_2P.fastq.gz",
        Genom="/mnt/ceph/ressources/01_refGenomes/Hsapiens/GRCh37/RefSeq/fasta/star_2.7.9",
            #sampleID="{sample}",
    output:
          "../03_star/{sample}/{sample}Aligned.sortedByCoord.out.bam",
    conda:
          "envs/star.yaml"
    threads: 8
    priority: 50
    shell:
            'STAR \
            --genomeDir {input.Genom} \
            --readFilesIn {input.r1} {input.r2} \
            --runThreadN {threads} \
            --outFileNamePrefix ../03_star/{wildcards.sample}/{wildcards.sample} \
            --readFilesCommand zcat \
            --alignIntronMax 500000 \
            --alignMatesGapMax 500000 \
            --outBAMcompression 0 \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMprimaryFlag OneBestScore \
            --outFilterMultimapNmax 100 \
            --outFilterMismatchNoverLmax 0.05 \
            --chimSegmentMin 15 \
            --chimOutType WithinBAM \
            --chimScoreMin 1 \
            --chimScoreJunctionNonGTAG 0 \
            --chimJunctionOverhangMin 15 \
            --chimSegmentReadGapMax 3 \
            --alignSJstitchMismatchNmax 5 -1 5 5'




rule Samtools:
    input:
        bam="../03_star/{sample}/{sample}Aligned.sortedByCoord.out.bam",
        r_umi="../02_trim/{sample}/{sample}_umis_1P.fastq.gz",

    output:
        index="../03_star/{sample}/{sample}Aligned.sortedByCoord.out.bam.bai",
        anno="../03_star/{sample}/{sample}_anno.bam",
        sorted="../03_star/{sample}/{sample}_fgbioSorted.bam",
        CallCons="../03_star/{sample}/{sample}_CallCons.bam",
        GroupRead="../03_star/{sample}/{sample}_GroupRead.bam",
        r1_out="../03_star/{sample}/UmiReads/{sample}_umiReads_1P.fastq.gz",
        r2_out="../03_star/{sample}/UmiReads/{sample}_umiReads_2P.fastq.gz",

    threads:8
    conda: "envs/samtools.yaml"
    shell:"samtools index {input.bam} ;\
           fgbio -XX:-UseGCOverheadLimit -Xms750m -Xmx32g --tmp-dir=../{wildcards.sample}/tmp/  AnnotateBamWithUmis -i {input.bam} -f {input.r_umi}  -o {output.anno} ;\
           fgbio -XX:-UseGCOverheadLimit -Xms750m -Xmx32g --tmp-dir=../{wildcards.sample}/tmp/ SortBam -i {output.anno} -s Queryname -o {output.sorted} ;\
           fgbio -XX:-UseGCOverheadLimit -Xms750m -Xmx32g --tmp-dir=../{wildcards.sample}/tmp/ SetMateInformation -i {output.sorted} -x true | \
           fgbio -XX:-UseGCOverheadLimit -Xms750m -Xmx32g --tmp-dir=../{wildcards.sample}/tmp/  GroupReadsByUmi --edits=1 --min-map-q=1 -t RX -s adjacency -o {output.GroupRead}; \
           fgbio -XX:-UseGCOverheadLimit -Xms750m -Xmx32g  --tmp-dir=../{wildcards.sample}/tmp/ CallMolecularConsensusReads --min-input-base-quality=2 --min-reads=1 --max-reads=1000000 \
           --output-per-base-tags=false --sort-order=:none: -i {output.GroupRead} -o {output.CallCons} ; \
           samtools fastq -1 {output.r1_out} -2 {output.r2_out} {output.CallCons} "





rule StarMapping2:
    input:
        r1="../03_star/{sample}/UmiReads/{sample}_umiReads_1P.fastq.gz",
        r2="../03_star/{sample}/UmiReads/{sample}_umiReads_2P.fastq.gz",
        Genom="/mnt/ceph/ressources/01_refGenomes/Hsapiens/GRCh37/RefSeq/fasta/star_2.7.9",
            #sampleID="{sample}",
    output:
          "../03_star/{sample}/UmiReads/{sample}Aligned.sortedByCoord.out.bam",
    conda:
          "envs/star.yaml"
    threads: 8
    priority: 50
    shell:
            'STAR \
            --genomeDir {input.Genom} \
            --readFilesIn {input.r1} {input.r2} \
            --runThreadN {threads} \
            --outFileNamePrefix ../03_star/{wildcards.sample}/UmiReads/{wildcards.sample} \
            --readFilesCommand zcat \
            --alignIntronMax 500000 \
            --alignMatesGapMax 500000 \
            --outBAMcompression 0 \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMprimaryFlag OneBestScore \
            --outFilterMultimapNmax 100 \
            --outFilterMismatchNoverLmax 0.05 \
            --chimSegmentMin 15 \
            --chimOutType WithinBAM \
            --chimScoreMin 1 \
            --chimScoreJunctionNonGTAG 0 \
            --chimJunctionOverhangMin 15 \
            --chimSegmentReadGapMax 3 \
            --alignSJstitchMismatchNmax 5 -1 5 5'


rule IndexSam:
    input:
          bam="../03_star/{sample}/UmiReads/{sample}Aligned.sortedByCoord.out.bam",
    output:
          bai="../03_star/{sample}/UmiReads/{sample}Aligned.sortedByCoord.out.bam.bai",
    conda: "envs/samtools.yaml"
    threads:8
    shell: 'samtools index {input.bam}'




rule FeatureCount:
    input:
         bam="../03_star/{sample}/UmiReads/{sample}Aligned.sortedByCoord.out.bam",
         bai="../03_star/{sample}/UmiReads/{sample}Aligned.sortedByCoord.out.bam.bai"
    output:
        out="../05_FeatureCount/{sample}_FeatureCountTable"
    conda:
        "envs/subread.yaml"
    threads:8
    priority: 50
    shell:
         "featureCounts -a /mnt/ceph/ressources/01_refGenomes/Hsapiens/GRCh37/RefSeq/GRCh37_latest_genomic_replace.gtf \
         -o {output.out} {input.bam}"


rule run_Arriba:
    input:
        bai="../03_star/{sample}/UmiReads/{sample}Aligned.sortedByCoord.out.bam.bai",
        bam="../03_star/{sample}/UmiReads/{sample}Aligned.sortedByCoord.out.bam"
    conda:
        "envs/arriba.yml"
    output:
        out="../04_Arriba/{sample}/{sample}.fusion.tsv"
    threads:4
    shell:
        "arriba -x {input.bam} -o {output.out} -a /mnt/ceph/ressources/01_refGenomes/Hsapiens/GRCh37/RefSeq/arriba_v2.3.0/hs37d5.fa \
        -g /mnt/ceph/ressources/01_refGenomes/Hsapiens/GRCh37/RefSeq/arriba_v2.3.0/gencode.v19.annotation.gtf \
        -b /mnt/ceph/ressources/01_refGenomes/Hsapiens/GRCh37/RefSeq/arriba_v2.3.0/blacklist_hg19_hs37d5_GRCh37_v2.3.0.tsv.gz -T -P -f duplicates;\
        "

rule Arriba_post:
    input:
        tab="../04_Arriba/{sample}/{sample}.fusion.tsv",

    conda: "envs/arriba.yml"
    output:
        done="../04_Arriba/{sample}/{sample}.done",
        tab2="../04_Arriba/{sample}/{sample}_filtered.fusion.tsv"

    threads:1
    shell: "python data/arriba.py {input.tab} {output.tab2} ../04_Arriba/{wildcards.sample}/{wildcards.sample}"

rule Arriba_plot:
    input:
        done="../04_Arriba/{sample}/{sample}.done",
        bai="../03_star/{sample}/UmiReads/{sample}Aligned.sortedByCoord.out.bam.bai"
    conda:
        "envs/arriba.yml"
    output:
        done2="../04_Arriba/{sample}/{sample}.done2"
    threads:1
    shell: "python data/arriba2.py {wildcards.sample}"


rule plotGepado:
    input:
        '../01_raw/{sample}/fastqc/{sample}_1P_fastqc.html',
        '../01_raw/{sample}/fastqc/{sample}_2P_fastqc.html',
        "../05_FeatureCount/{sample}_FeatureCountTable"

    output:
        out='../06_quality/{sample}.pdf'
    conda:
        'envs/plot.yml'
    shell:
        'python data/plt.py {wildcards.sample} '
