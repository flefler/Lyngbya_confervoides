# LYNGBYA script
Last updated July 9 2025 by Forrest W. Lefler

# OSCILLATORIALES

## OSCILLATORIALES Alignment
```
mamba activate gtotree2
N=GTOTREE_OSCILLATORIALES

CMD="cd OSCILLATORIALES &&\
GToTree -f OSCILLATORIALES.txt \
-m OSCILLATORIALES_mapping.txt \
-H Cyanobacteria \
-o TREE \
-B -n 2 -j 8 -k -M 16 -T IQ-TREE"

sbatch -A hlaughinghouse --qos=hlaughinghouse-b -J ${N} -c 16 --mem=100G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt \
--export=ALL --mail-type=FAIL --mail-user=flefler@ufl.edu -t 24:00:00 --wrap="${CMD}"
```

## Par genes
```
module load pargenes

N=Pargenes_OSCILLATORIALES

CMD="pargenes.py --seed 42069 \
-a LYNGBYA/OSCILLATORIALES2/TREE/run_files/individual_alignments \
-o LYNGBYA/OSCILLATORIALES2/TREE/Par_genes \
-c 32 -d aa -b 1000 --autoMRE --use-modeltest --modeltest-criteria AICc --modeltest-perjob-cores 8 --use-astral --continue"

sbatch -A hlaughinghouse --qos=hlaughinghouse-b -J ${N} -c 32 --mem=100G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt \
--export=ALL --mail-type=FAIL --mail-user=flefler@ufl.edu -t 48:00:00 --wrap="${CMD}"
```

## weighted astral
```
module load aster
astral-weighted -a OSCILLATORIALES2/OSCILLATORIALES_mapping.txt -v 1 -R -S -u 1 -i OSCILLATORIALES2/TREE/Par_genes/astral_run/gene_trees.newick -o OSCILLATORIALES2/TREE/OSCILLATORIALESAASTER.nwk
```

## SKANI
```
mamba activate skani
mkdir OSCILLATORIALES2/ANI
skani dist -t 16 --ql OSCILLATORIALES2/OSCILLATORIALES.txt --rl OSCILLATORIALES2/OSCILLATORIALES.txt -s 50 --slow -o OSCILLATORIALES2/ANI/all-to-all_results.txt
```

## ezAAI
### Extract protein sequences
```
mamba activate ezaai
mkdir -p OSCILLATORIALES2/AAI/{db,FILES}
SAMPLES=`cut -f 1 OSCILLATORIALES2/OSCILLATORIALES.txt`
for SAMPLE in $SAMPLES ; do

    filename=$(basename "$SAMPLE")
    N="${filename%%.*}"

    CMD="seqkit seq -i ${SAMPLE} -o OSCILLATORIALES2/AAI/FILES/${N}.fasta &&\
    EzAAI extract -i OSCILLATORIALES2/AAI/FILES/${N}.fasta -o OSCILLATORIALES2/AAI/db/${N}.faa"

    sbatch -A hlaughinghouse --qos=hlaughinghouse-b -J ${N} -c 1 --mem=10G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt \
    --export=ALL --mail-type=FAIL --mail-user=flefler@ufl.edu -t 1:00:00 --wrap="${CMD}"

done
```

### Run AAI calculation
```
mamba activate ezaai
N=AAI_OSCILLATORIALES
rm -rf OSCILLATORIALES2/AAI/FILES

CMD="EzAAI calculate -i OSCILLATORIALES2/AAI/db/ -j OSCILLATORIALES2/AAI/db/ -mtx -t 32 -o OSCILLATORIALES2/AAI/AAI"
sbatch -A hlaughinghouse --qos=hlaughinghouse-b -J ${N} -c 32 --mem=100G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt \
--export=ALL --mail-type=FAIL --mail-user=flefler@ufl.edu -t 73:00:00 --wrap="${CMD}"
```

# Okeania

## Okeania Alignment
```
mamba activate gtotree2
N=GTOTREE_Okeania

CMD="cd OKEANIA &&\
GToTree -f Okeania.txt \
-m Okeania_mapping.txt \
-H Cyanobacteria \
-o TREE \
-B -n 2 -j 8 -k -M 16 -T IQ-TREE"
```

## Par genes
```
mamba activate pargenes
N=Okeania_pargenes

CMD="pargenes.py --seed 42069 \
-a LYNGBYA/OKEANIA/TREE/run_files/individual_alignments \
-o LYNGBYA/OKEANIA/TREE/Par_genes \
-c 32 -d aa -b 1000 --autoMRE --use-modeltest --modeltest-criteria AICc --modeltest-perjob-cores 8 --use-astral"
```

## weighted astral
```
module load aster
astral-weighted -a OKEANIA/Okeania_mapping2.txt -v 1 -R -S -u 1 -i OKEANIA/TREE/Par_genes/astral_run/gene_trees.newick -o OKEANIA/TREE/OKEANIAASTER.nwk
```

## SKANI
```
mamba activate skani
mkdir OKEANIA/ANI
skani dist -t 16 --ql OKEANIA/Okeania.txt --rl OKEANIA/Okeania.txt -s 50 --slow -o OKEANIA/ANI/all-to-all_results.txt
```

# Capilliphycus, Limnoraphis, Lyngbya

## Alignment
```
mamba activate gtotree2
N=GTOTREE_CLL
CMD="cd CLL &&\
GToTree -f CLL.txt \
-m CLL_mapping2.txt \
-H Cyanobacteria \
-o TREE \
-B -n 2 -j 8 -k -M 16 -T IQ-TREE"
```

## Par genes
```
mamba activate pargenes
module load pargenes

N=Pargenes_Limnoraphis

CMD="pargenes.py --seed 42069 \
-a LYNGBYA/CLL/TREE/run_files/individual_alignments \
-o LYNGBYA/CLL/TREE/Par_genes \
-c 64 -d aa -p 10 -b 1000 --autoMRE --use-modeltest --modeltest-criteria AICc --modeltest-perjob-cores 8 \
--use-astral"
```

## weighted astral
```
module load aster
astral-weighted -a CLL_mapping2.txt -v 1 -R -S -u 1 -i TREE/Par_genes/astral_run/gene_trees.newick -o TREE/CLL_ASTER.nwk
```

## SKANI
```
mamba activate skani
mkdir CLL/ANI
skani dist -t 16 --ql CLL/CLL.txt --rl CLL/CLL.txt -s 50 --medium -o CLL/ANI/all-to-all_results.txt --slow
```

## ezAAI
### Extract protein sequences
```
mamba activate ezaai
mkdir -p CLL/AAI/{db,FILES}
SAMPLES=`cut -f 1 CLL/CLL.txt`
for SAMPLE in $SAMPLES ; do

    filename=$(basename "$SAMPLE")
    N="${filename%%.*}"

    CMD="seqkit seq -i ${SAMPLE} -o CLL/AAI/FILES/${N}.fasta &&\
    EzAAI extract -i CLL/AAI/FILES/${N}.fasta -o CLL/AAI/db/${N}.faa"

done
```
### Run AAI calculation
```
mamba activate ezaai
N=AAI_CLL
rm -rf CLL/AAI/FILES

CMD="EzAAI calculate -i CLL/AAI/db/ -j CLL/AAI/db/ -mtx -t 32 -o CLL/AAI/AAI"
```




