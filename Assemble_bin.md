# Assembly Script
Last updated Oct 14 2025 by Forrest W. Lefler

# Assemble ONT reads with myloasm
We will use the mmlong2-lite package for this.
It will assemble with myloasm, polish with medaka, and check for misassemblies with rust-anvio-mis.
It does more, but we will stop it here because we will polish with short reads, if we have, and check for MGE contigs with geNomad before binning. Two steps not included in mmlong2-lite.
```
mamba activate mmlong
for SAMPLE in M137 M349 M348 ; do

    N=mmlong_${SAMPLE}
    mkdir -p 02_mmlong/${SAMPLE}

    CMD="cd 02_mmlong/${SAMPLE} &&\
    mmlong2-lite -np /blue/hlaughinghouse/flefler/CYANO_GENOMES/00_LONGREADS/${SAMPLE}/*fastq \
    -p 32 -myl -med --run_until curation"

done
```

##  Polish contigs with PYPOLCA
```
mamba activate Longreads
for SAMPLE in M348 M349 ; do

    N=${SAMPLE}_polish

    F=00_READS/${SAMPLE}/*1.fastq.gz
    R=00_READS/${SAMPLE}/*2.fastq.gz

    CONTIGS=02_mmlong/${SAMPLE}/mmlong2/tmp/curation/asm_curated.fasta

    CMD="pypolca run -t 16 -f \
        -1 ${F} -2 ${R} \
        -a ${CONTIGS} \
        -o 02_mmlong/${SAMPLE}/PYPOLCA/ &&\
        cp 02_mmlong/${SAMPLE}/PYPOLCA/pypolca_corrected.fasta 02_mmlong/${SAMPLE}/polished_contigs.fasta"
one
```

## Filter illumina reads which do NOT map to the myloasm assembly
We will reassemble those
```
mamba activate Longreads
for SAMPLE in M348 M349 ; do

    mkdir -p 02_mmlong/${SAMPLE}/UNMAPPED/FILTERING/READS

    N=${SAMPLE}filter

    F=00_READS/${SAMPLE}/*1.fastq.gz
    R=00_READS/${SAMPLE}/*2.fastq.gz

    CMD="coverm make -t 32 \
        -r 02_mmlong/${SAMPLE}/polished_contigs.fasta \
        -1 ${F} \
        -2 ${R} \
        -o 02_mmlong/${SAMPLE}/UNMAPPED/FILTERING &&\
        coverm filter -t 32 --proper-pairs-only \
        -b 02_mmlong/${SAMPLE}/UNMAPPED/FILTERING/*.bam -o 02_mmlong/${SAMPLE}/UNMAPPED/FILTERING/inverse_filtered.bam --inverse &&\
        samtools sort -@ 32 -n 02_mmlong/${SAMPLE}/UNMAPPED/FILTERING/inverse_filtered.bam -o 02_mmlong/${SAMPLE}/UNMAPPED/FILTERING/inverse_filtered_sorted.bam &&\
        bedtools bamtofastq -i 02_mmlong/${SAMPLE}/UNMAPPED/FILTERING/inverse_filtered_sorted.bam \
        -fq 02_mmlong/${SAMPLE}/UNMAPPED/FILTERING/READS/unmapped_R1.fastq \
        -fq2 02_mmlong/${SAMPLE}/UNMAPPED/FILTERING/READS/unmapped_R2.fastq"

done
```

## metaSPADES
Assemble the illumina reads which dont map to the myloasm assembly, remove contigs less than 2000 bp, cat with myloasm assembly
```
mamba activate Longreads
for SAMPLE in M348 M349 ; do

    N=${SAMPLE}_metaSPAdes

    F=02_mmlong/${SAMPLE}/UNMAPPED/FILTERING/READS/unmapped_R1.fastq
    R=02_mmlong/${SAMPLE}/UNMAPPED/FILTERING/READS/unmapped_R2.fastq
    LONG=00_LONGREADS/${SAMPLE}/*.fastq

    CMD="spades.py --meta \
        -1 ${F} -2 ${R} \
        --nanopore ${LONG} \
        -k 21,33,55,77,99,127 \
        -o 02_mmlong/${SAMPLE}/SPADES \
        -m 250 -t 16 &&\
        seqkit seq -m 2000 02_mmlong/${SAMPLE}/SPADES/scaffolds.fasta -o 02_mmlong/${SAMPLE}/SPADES/cleaned_scaffolds.fasta &&\
        zcat 02_mmlong/${SAMPLE}/SPADES/cleaned_scaffolds.fasta.gz >> 02_mmlong/${SAMPLE}/polished_contigs.fasta"

done
```

# Assemble Illumina only
```
mamba activate Longreads
for SAMPLE in M309 M399 M182 ; do

    N=${SAMPLE}_metaSPAdes

    F=00_READS/${SAMPLE}/*_1.fastq.gz
    R=00_READS/${SAMPLE}/*_2.fastq.gz

    CMD="spades.py --meta \
        -1 ${F} -2 ${R} \
        -k 21,33,55,77,99,127 \
        -o 02_SPADES/${SAMPLE} \
        -m 250 -t 16 &&\
        seqkit seq -m 2000 02_mmlong/${SAMPLE}/SPADES/scaffolds.fasta -o 02_mmlong/${SAMPLE}/SPADES/cleaned_scaffolds.fasta &&\
        zcat 02_mmlong/${SAMPLE}/SPADES/cleaned_scaffolds.fasta.gz >> 02_mmlong/${SAMPLE}/polished_contigs.fasta"

done
```

# geNomad
Identify viral and plasmid contigs with geNomad
Helps reduce the amount of data we need to bin, speeds things up
```
mamba activate genomad
for SAMPLE in M309 M399 M182 M137 M348 M349 ; do
    N=${SAMPLE}_geNomad

    CMD="genomad end-to-end --cleanup --disable-find-proviruses --quiet --enable-score-calibration --threads 16 \
    --lenient-taxonomy --full-ictv-lineage --composition metagenome --min-score 0.7 \
    02_mmlong/${SAMPLE}/finalcontigs.fasta 02_mmlong/${SAMPLE}/04_geNomad/ /orange/hlaughinghouse/flefler/databases/genomad_db &&\
    seqkit grep -v -f <(cat 02_mmlong/${SAMPLE}/04_geNomad/*_summary/*_{virus,plasmid}_summary.tsv \
    | cut -f 1 | grep -v "seq_name" | sort | uniq) 02_mmlong/${SAMPLE}/finalcontigs.fasta \
    -o 02_mmlong/${SAMPLE}/Prok_final.fasta"

done
```

# BIN
## Classify contigs METABULI
This really isnt that useful to try and interperate, but we need it for VAMB.
```
mamba activate metabuli
for SAMPLE in M309 M399 M182 M137 M348 M349 ; do
    N=${SAMPLE}_metabuli

    mkdir -p 02_mmlong/${SAMPLE}/02_METABULI/

    CONTIGS=02_mmlong/${SAMPLE}/Prok_final.fasta

    CMD="metabuli classify \
        --seq-mode 3 \
        ${CONTIGS} \
        /orange/hlaughinghouse/flefler/databases/METABULI/gtdb \
        02_mmlong/${SAMPLE}/02_METABULI/ \
        ${SAMPLE} \
        --min-score 0.008 &&\
        taxconverter metabuli \
        -c 02_mmlong/${SAMPLE}/02_METABULI/${SAMPLE}_classifications.tsv \
        -r 02_mmlong/${SAMPLE}/02_METABULI/${SAMPLE}_report.tsv \
        -o 02_mmlong/${SAMPLE}/02_METABULI/${SAMPLE}_result.tsv"

done
```

## Make BAM files, long reads
```
mamba activate Longreads
for SAMPLE in M137 M348 M349 ; do

    mkdir -p 02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/BAM

    N=${SAMPLE}_BAM

    CONTIGS=02_mmlong/${SAMPLE}/Prok_final.fasta
    LONG=00_LONGREADS/${SAMPLE}/*fastq

    CMD="coverm make -t 16 \
    -r ${CONTIGS} \
    --single ${LONG} \
    -o 02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/BAM \
    -p minimap2-ont &&\
    mv 02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/BAM/*.bam 02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/BAM/${SAMPLE}.bam"

    sbatch -A hlaughinghouse --qos=hlaughinghouse-b -J ${N} -c 16 --mem=20G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt \
    --export=ALL --mail-type=FAIL --mail-user=flefler@ufl.edu -t 1:00:00 --wrap="${CMD}"

done
```

## Make BAM files, short reads
```
mamba activate Longreads
for SAMPLE in M309 M399 M182 M348 M349 ; do

    CONTIGS=02_mmlong/${SAMPLE}/Prok_final.fasta

    mkdir -p 02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/BAM/short

    N=${SAMPLE}_sBAM

    F=00_READS/${SAMPLE}/*1.fastq.gz
    R=00_READS/${SAMPLE}/*2.fastq.gz

    CMD="coverm make -t 16 \
    -r ${CONTIGS} -1 ${F} -2 ${R} \
    -o 02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/BAM/short &&\
    mv 02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/BAM/short/*.bam 02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/BAM/${SAMPLE}_short.bam"

    sbatch -A hlaughinghouse --qos=hlaughinghouse-b -J ${N} -c 16 --mem=20G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt \
    --export=ALL --mail-type=FAIL --mail-user=flefler@ufl.edu -t 1:00:00 --wrap="${CMD}"
done
```

## comebin
```
module load comebin
for SAMPLE in M309 M399 M182 M137 M348 M349 ; do
    mkdir -p 02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/comebin

    N=${SAMPLE}_comebin

    CONTIGS=02_mmlong/${SAMPLE}/Prok_final.fasta

    OUT_DIR=02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/comebin

    CMD="run_comebin.sh \
    -a ${CONTIGS} \
    -o ${OUT_DIR} \
    -p 02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/BAM/ \
    -t 32"

done
```
## Semibin
```
module load semibin
for SAMPLE in M309 M399 M182 M137 M348 M349 ; do

    mkdir -p 02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/

    N=${SAMPLE}_Semibin

    CONTIGS=02_mmlong/${SAMPLE}/Prok_final.fasta

    OUT_DIR=02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/Semibin

    CMD="SemiBin2 single_easy_bin \
    --random-seed 42069 \
    --engine gpu \
    --sequencing-type long_read \
    -i ${CONTIGS} \
    -b 02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/BAM/*.bam \
    -o ${OUT_DIR}"

done
```
## METADECODER
```
mamba activate metadecoder
SAMPLES=`cut -f 1 hybrid.txt | sed '1d'`
for SAMPLE in M309 M399 M182 M137 M348 M349 ; do
    mkdir -p 02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/METADECODER/FILES
    mkdir -p 02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/METADECODER/BINS

    CONTIGS=/blue/hlaughinghouse/flefler/CYANO_GENOMES/02_mmlong/${SAMPLE}/Prok_final.fasta

    HERE=/blue/hlaughinghouse/flefler/CYANO_GENOMES/02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/METADECODER

    BAM=/blue/hlaughinghouse/flefler/CYANO_GENOMES/02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/BAM/*.bam

    N=${SAMPLE}_metadecoder

    CMD="cd ${HERE} &&\
    metadecoder coverage --threads 16 -b ${BAM} -o FILES/METADECODER.COVERAGE &&\
    metadecoder seed --threads 16 -f ${CONTIGS} -o FILES/METADECODER.SEED &&\
    metadecoder cluster -f ${CONTIGS} -c FILES/METADECODER.COVERAGE -s FILES//METADECODER.SEED -o ${SAMPLE} &&\
    gzip *.fasta && mv *.fasta.gz BINS/ && mv *.{dpgmm,kmers} FILES/"

done
```

## VAMB
```
mamba activate VAMB
for SAMPLE in M309 M399 M182 M137 M348 M349 ; do

    rm -rf 02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/VAMB

    CONTIGS=02_mmlong/${SAMPLE}/Prok_final.fasta

    N=${SAMPLE}_VAMB

    OUT_DIR=02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/VAMB

    CMD="vamb bin taxvamb -m 2000 -p 32 \
    --seed 42069 \
    --minfasta 300000 \
    --outdir ${OUT_DIR} \
    --fasta ${CONTIGS} \
    --taxonomy 02_mmlong/${SAMPLE}/02_METABULI/${SAMPLE}_result.tsv \
    --bamdir 02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/BAM/"

done
```

## binette
```
mamba activate binette
for SAMPLE in M309 M399 M182 M137 M348 M349 ; do

    N=${SAMPLE}_binette

    CONTIGS=02_mmlong/${SAMPLE}/Prok_final.fasta

    CMD="binette -t 16 \
    --bin_dirs \
    02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/Semibin/output_bins/ \
    02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/VAMB/bins/ \
    02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/comebin/comebin_res/comebin_res_bins/ \
    02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/METADECODER/BINS/ \
    -c ${CONTIGS} \
    -m 90 \
    -o 02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/binette"

done
```


## Collect all MAGs in one place
```
mkdir -p 02_mmlong/99_bins
for SAMPLE in M309 M399 M182 M137 M348 M349 ; do

    outdir="02_mmlong/99_bins"

    # 1. Copy and rename all bin FASTA files from binette
    bin_dir="02_mmlong/${SAMPLE}/05_BIN/${SAMPLE}/binette/final_bins"
    for file in "$bin_dir"/*.fa; do
        # Check if any files exist to avoid errors
        [ -e "$file" ] || continue
        # Make a new basename for output file
        newbase="${SAMPLE}_$(basename "$file" .fa).fasta"
        cp "$file" "$outdir/$newbase"
    done

    # 2. Rename contig headers in each output FASTA file
    for f in "$outdir"/*.fasta; do
        prefix=$(basename "$f" .fasta)
        awk -v p="$prefix" 'BEGIN{c=0} /^>/ {c++; printf ">%s_contig_%05d\n", p, c} !/^>/ {print}' "$f" > "$outdir/${prefix}.tmp"
        mv "$outdir/${prefix}.tmp" "$f"
    done
done
```

```
