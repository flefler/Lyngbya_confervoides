# QA/QC for raw genomic reads
Last updated Oct 30 2025 by Forrest W. Lefler

# Illumina reads
Runs fastP and decontaminates human reads with BowTie2
fastP removes poly G ```-g``` and ploy X ```-x``` tails (common with novaseq), detects for adaptor sequences, and removes reads shorter than 150bp
```
mamba activate QAQC_Assemble
for SAMPLE in M182 M309 M348 M349 M399; do

    N=${SAMPLE}

    F=00_RAWREADS/${SAMPLE}/*1.fq.gz
    R=00_RAWREADS/${SAMPLE}/*2.fq.gz
    REF=/homo_mus_phix

    mkdir -p 00_READS/${SAMPLE}
    mkdir -p 01_QC/FASTP

    CMD="fastq_preprocessor_short.py \
    --forward_reads ${F} --reverse_reads ${R} \
    -n ${N} -o TMP/ \
    -x ${REF} -p 16 -m 150 --fastp_options '-D -x -g --detect_adapter_for_pe -c -j -y -l 150 -j 01_QC/FASTP/${SAMPLE}.json' &&\
    mv TMP/${SAMPLE}/intermediate/2__bowtie2/*1.fastq.gz 00_READS/${SAMPLE}/${SAMPLE}_1.fastq.gz &&\
    mv TMP/${SAMPLE}/intermediate/2__bowtie2/*2.fastq.gz 00_READS/${SAMPLE}/${SAMPLE}_2.fastq.gz &&\
    rm -rf TMP/${SAMPLE}"

done
```

# ONT reads
just run fastP long
Removes reads shorter than 500bp, and keeps only >Q20 reads
```
mamba activate QAQC_Assemble
mkdir -p 01_QC/FASTPLong
for SAMPLE in M137 M349 M348 ; do

    N=${SAMPLE}_long

    mkdir -p 00_LONGREADS/${SAMPLE}

    READS=00_RAWLONGREADS/${SAMPLE}/${SAMPLE}.fastq.gz

    CMD="fastplong -i ${READS} -o 00_LONGREADS/${SAMPLE}/${SAMPLE}.fastq.gz -x -w 16 -l 500 -m 20 -M 20 -5 -3 -y \
    -j 01_QC/FASTPLong/${SAMPLE}_long.json -h 01_QC/FASTPLong/${SAMPLE}_long.html"

done
```
