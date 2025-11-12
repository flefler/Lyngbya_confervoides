# Analyses of MAGs

## CheckM2
```
mamba activate binette
export CHECKM2DB="/orange/hlaughinghouse/flefler/databases/CheckM2_1.1.0/CheckM2_database/uniref100.KO.1.dmnd"

N=checkm2

mkdir -p 02_mmlong/07_INFO/

CMD="checkm2 predict --quiet \
--threads 32 \
--input 02_mmlong/99_bins/ \
--output-directory 02_mmlong/07_INFO/CHECKM2 \
-x .fasta"
```

## GTDBtk
```
mamba activate gtdbtk_2.4.1

N=GTDB

CMD="gtdbtk classify_wf \
--genome_dir 02_mmlong/99_bins \
--out_dir 02_mmlong/07_INFO/GTDB \
--mash_db GTDB/release226/gtdb_r226.msh \
-x fasta --force \
--cpus 32 --pplacer_cpus 32"
```

## coverM
Long only samples
```
mamba activate Longreads
mkdir -p 02_mmlong/07_INFO/coverM
for SAMPLE in M137 ; do

    N=${SAMPLE}_coverM

    LONG=00_LONGREADS/${SAMPLE}/*.fastq

    CMD="coverm genome \
    --single ${LONG} \
    --genome-fasta-files 02_mmlong/99_bins/${SAMPLE}*.fasta \
    -x .fasta --threads 8 --methods mean relative_abundance \
    -o 02_mmlong/07_INFO/coverM/${SAMPLE}_coverM.tsv"

    sbatch -A hlaughinghouse --qos=hlaughinghouse-b -J ${N} -c 8 --mem=20G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt \
    --export=ALL --mail-type=NONE --mail-user=flefler@ufl.edu -t 1:00:00 --wrap="${CMD}"
    
done
```
Hybrid samples
```
mamba activate Longreads
for SAMPLE in M182 M309 M348 M349 M399; do

    N=${SAMPLE}_coverM

    F=00_READS/${SAMPLE}/*1.fastq.gz
    R=00_READS/${SAMPLE}/*2.fastq.gz
    LONG=00_LONGREADS/${SAMPLE}/*.fastq

    CMD="coverm genome \
    --read1 ${F} \
    --read2 ${R} \
    --single ${LONG} \
    --genome-fasta-files 02_mmlong/99_bins/${SAMPLE}*.fasta \
    -x .fasta --threads 8 --methods mean relative_abundance \
    -o 02_mmlong/07_INFO/coverM/${SAMPLE}_coverM.tsv"

    sbatch -A hlaughinghouse --qos=hlaughinghouse-b -J ${N} -c 8 --mem=20G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt \
    --export=ALL --mail-type=NONE --mail-user=flefler@ufl.edu -t 1:00:00 --wrap="${CMD}"
    
done
```

## BAKTA
```
mamba activate bakta
mkdir -p 02_mmlong/07_INFO/BAKTA

for SAMPLE in 02_mmlong/99_bins/*.fasta; do

    prefix=$(basename $SAMPLE .fasta)

    N=${prefix}_bakta

    CMD="bakta --skip-plot  \
    --min-contig-length 5000 \
    --keep-contig-headers \
    --db /orange/hlaughinghouse/flefler/databases/bakta/db \
    --output 02_mmlong/07_INFO/BAKTA/${prefix} \
    --prefix ${prefix} \
    ${SAMPLE}"

    sbatch -A hlaughinghouse --qos=hlaughinghouse-b -J ${N} -c 16 --mem=100G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt \
    --export=ALL --mail-type=FAIL --mail-user=flefler@ufl.edu -t 4:00:00 --wrap="${CMD}"

done
```

