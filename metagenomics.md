# Bacterial metagenomic/WGS sequencing analysis

Yiming Qu

Apr 2nd, 2025

## Table of Contents

1. [Install all dependencies upfront](#install-all-dependencies-upfront)
2. [Raw sequence quality control](#raw-sequence-quality-control)
3. [Removal of host contamination (for metagenomic sequencing only)](#removal-of-host-contamination-for-metagenomic-sequencing-only)
4. [Genome assembly and annotation](#genome-assembly-and-annotation)
5. [Map reads to the assembled genome](#map-reads-to-the-assembled-genome)
6. [Quality assessment of assembled genomes.](#quality-assessment-of-assembled-genomes)
7. [GTDB-Tk taxonomy assignment (for isolates only)](#gtdb-tk-taxonomy-assignment-for-isolates-only)
8. [Metaphlan taxonomy profile](#metaphlan-taxonomy-profile)
9. [kraken2 taxonomy profile](#kraken2-taxonomy-profile)
10. [Metagenomic binning](#metagenomic-binning)
11. [DRAM protein functional profile](#dram-protein-functional-profile)
12. [HUMAnN functional annotation](#humann-functional-annotation)
13. [Metabolomic modeling](#metabolomic-modeling)

## Install all dependencies upfront

Copy the following content into an environment file, say `./environment.yml` and run `mamba env create -f environment.yml` to setup the environment.

```
name: <env_name>
channels:
  - bioconda
  - conda-forge
dependencies:
  - python=3.10
  - pip
  - seqkit
  - fastqc
  - fastp
  - spades
  - prokka
  - bwa-mem2
  - samtools
  - sambamba
  - quast
  - checkm-genome=1.2.3
  - pandas
  - gtdbtk
  - metaphlan
  - kraken2
  - bracken
  - metabat2
  - dram
  - pip:
      - humann
```

## Raw sequence quality control

Assuming your raw fastq files are in `$DATA_DIR/fastq` with file names `*_R1_001.fastq.gz` for read 1 and `*_R2_001.fastq.gz` for read 2. Ideally, you can directly take the output from Illumina demultiplexing, no need for renaming.

```bash
#!/bin/bash

# Define directories
DATA_DIR=<your_data_dir>

# Create output subdirectories
mkdir -p "$DATA_DIR/fastqc_before"
mkdir -p "$DATA_DIR/fastqc_after"
mkdir -p "$DATA_DIR/fastp"

seqkit stats -j 8 $DATA_DIR/fastq/*.fastq.gz > "$DATA_DIR/input.stats"

# Loop over all paired-end files and process sequentially
for READ1 in $DATA_DIR/fastq/*_R1_001.fastq.gz; do
    # Derive READ2 and sample name
    READ2="${READ1%_R1_001.fastq.gz}_R2_001.fastq.gz"
    BASENAME=$(basename "$READ1")
    SAMPLE_NAME="${BASENAME%%_*}"

    echo "Processing sample: $SAMPLE_NAME"

    # Step 1: FastQC before fastp
    echo "Running FastQC before fastp for $SAMPLE_NAME..."
    fastqc -o "$DATA_DIR/fastqc_before" "$READ1" "$READ2"

    # Step 2: Run fastp for adapter trimming and quality filtering
    echo "Running fastp for $SAMPLE_NAME..."
    fastp \
        -i "$READ1" \
        -I "$READ2" \
        -o "$DATA_DIR/fastp/${SAMPLE_NAME}_1.fq.gz" \
        -O "$DATA_DIR/fastp/${SAMPLE_NAME}_2.fq.gz" \
        --cut_tail \
        --html "$DATA_DIR/fastp/${SAMPLE_NAME}.html" \
        --json "$DATA_DIR/fastp/${SAMPLE_NAME}.json"

    # Step 3: FastQC after fastp
    echo "Running FastQC after fastp for $SAMPLE_NAME..."
    fastqc -o "$DATA_DIR/fastqc_after" \
        "$DATA_DIR/fastp/${SAMPLE_NAME}_1.fq.gz" \
        "$DATA_DIR/fastp/${SAMPLE_NAME}_2.fq.gz"

    echo "Sample $SAMPLE_NAME processing complete."
done

seqkit stats -j 8 $DATA_DIR/fastp/*.fq.gz > "$DATA_DIR/fastp.stats"

```

## Removal of host contamination (for metagenomics only)

For metagenomic sequencing on samples with, host contamination is very common (can account for 10~95% of sequencing data) and in most cases not helpful for downstream analysis. Therefore, we can remove these host DNA to get cleaner data to work with later.

```bash
#!/usr/bin/env bash

DATA_DIR=<your_data_dir>
REF_GENOME=<path_to_ref_genome_with_index>
THREADS=16

# Create output directories
mkdir -p "$DATA_DIR/map_host"
mkdir -p "$DATA_DIR/fastp_no_host"

# Index the reference genome if needed
# if [ ! -f "${REF_GENOME%.gz}.0127.bwt.2bit.64" ]; then
#     echo "Indexing the reference genome with bwa-mem2..."
#     gunzip -c "$REF_GENOME" > "${REF_GENOME%.gz}"
#     bwa-mem2 index "${REF_GENOME%.gz}"
# fi

for READ1 in "$DATA_DIR/fastp/"*_1.fq.gz; do
    READ2="${READ1%_1.fq.gz}_2.fq.gz"
    BASENAME=$(basename "$READ1")
    SAMPLE="${BASENAME%_1.fq.gz}"

    BAM_OUT="$DATA_DIR/map_host/${SAMPLE}.bam"
    UNMAPPED1="$DATA_DIR/fastp_no_host/${SAMPLE}_1.fq.gz"
    UNMAPPED2="$DATA_DIR/fastp_no_host/${SAMPLE}_2.fq.gz"

    # Align to host genome
    bwa-mem2 mem -t "$THREADS" "${REF_GENOME%.gz}" "$READ1" "$READ2" | \
        samtools view -@ "$THREADS" -bS -o "$BAM_OUT" -

    # Extract unmapped read pairs
    samtools fastq -@ "$THREADS" \
        -f 12 -F 256 "$BAM_OUT" \
        -1 "$UNMAPPED1" -2 "$UNMAPPED2" \
        -0 /dev/null -s /dev/null -n
done

# Generate stats on fastp output
seqkit stats -j 8 "$DATA_DIR/fastp_no_host/"*.fq.gz > "$DATA_DIR/fastp_no_host.stats"
```

In the following sections, I assume that we still use sequencing files in `fastp` instead of `fastp_no_host`. Change scripts accordingly if `fastp_no_host` better fits your need.

## Genome assembly and annotation

Taking the cleaned fastq files in the `fastp` folder you got from the previous step, here we run `spades` to assemble genomes and `prokka` to annotate proteins.

```bash
#!/usr/bin/env bash

DATA_DIR=<your_data_dir>
THREADS=16

mkdir -p "$DATA_DIR/spades"
mkdir -p "$DATA_DIR/assembly"
mkdir -p "$DATA_DIR/prokka"

for READ1 in "$DATA_DIR/fastp/"*_1.fq.gz; do
    READ2="${READ1%_1.fq.gz}_2.fq.gz"
    BASENAME=$(basename "$READ1")
    SAMPLE_NAME="${BASENAME%_1.fq.gz}"

    if [[ ! -f "$READ2" ]]; then
        echo "Warning: No matching R2 file found for $READ1. Skipping..."
        continue
    fi

    echo "Running SPAdes for $SAMPLE_NAME"
    SAMPLE_OUT="$DATA_DIR/spades/${SAMPLE_NAME}_assembly"
    mkdir -p "$SAMPLE_OUT"

    spades.py \
        --isolate \
        -1 "$READ1" \
        -2 "$READ2" \
        -o "$SAMPLE_OUT" \
        --threads $THREADS

    cp "$SAMPLE_OUT/contigs.fasta" "$DATA_DIR/assembly/${SAMPLE_NAME}.fna"

    echo "Running Prokka for $SAMPLE_NAME"
    PROKKA_OUT="$DATA_DIR/prokka/${SAMPLE_NAME}"
    mkdir -p "$PROKKA_OUT"

    prokka \
        --rnammer \
        --rfam \
        --cpus $THREADS \
        --outdir "$PROKKA_OUT" \
        --prefix "$SAMPLE_NAME" \
        --force \
        --centre X --compliant \
        "$DATA_DIR/assembly/${SAMPLE_NAME}.fna"
done
```

## Map reads to the assembled genome

Here we use `bwa-mem2` to map reads to assembled genomes, so we will have an idea of coverage and strain heterogeneity.

```bash
#!/usr/bin/env bash

DATA_DIR=<your_data_dir>
THREADS=16

mkdir -p "$DATA_DIR/bwa"

for FNA in "$DATA_DIR/assembly/"*.fna; do
    SAMPLE_NAME=$(basename "$FNA" .fna)
    READ1="$DATA_DIR/fastp/${SAMPLE_NAME}_1.fq.gz"
    READ2="$DATA_DIR/fastp/${SAMPLE_NAME}_2.fq.gz"
    BWA_OUT="$DATA_DIR/bwa/$SAMPLE_NAME"

    if [[ ! -f "$READ1" || ! -f "$READ2" ]]; then
        echo "Missing reads for $SAMPLE_NAME. Skipping..."
        continue
    fi

    echo "Indexing $SAMPLE_NAME assembly"
    bwa-mem2 index "$FNA"

    echo "Mapping reads for $SAMPLE_NAME"
    mkdir -p "$BWA_OUT"
    bwa-mem2 mem -t "$THREADS" "$FNA" "$READ1" "$READ2" \
        | samtools view -b -@ "$THREADS" -o "$BWA_OUT/${SAMPLE_NAME}.bam" -
    
    sambamba sort "$BWA_OUT/${SAMPLE_NAME}.bam" -o "$BWA_OUT/${SAMPLE_NAME}.sorted.bam"
done
```

## Quality assessment of assembled genomes (for isolates only)

We use `quast.py` to assess assembly quality (it also requires mapping output generated above) and `checkm` to assess phylogenetic assembly quality.

```bash
#!/usr/bin/env bash

DATA_DIR=<your_data_dir>
THREADS=16

mkdir -p "$DATA_DIR/quast"
mkdir -p "$DATA_DIR/checkm"

for FNA in "$DATA_DIR/assembly/"*.fna; do
    SAMPLE_NAME=$(basename "$FNA" .fna)
    READ1="$DATA_DIR/fastp/${SAMPLE_NAME}_1.fq.gz"
    READ2="$DATA_DIR/fastp/${SAMPLE_NAME}_2.fq.gz"
    PROKKA_OUT="$DATA_DIR/prokka/$SAMPLE_NAME"
    BWA_BAM="$DATA_DIR/bwa/$SAMPLE_NAME/${SAMPLE_NAME}.sorted.bam"
    QUAST_OUT="$DATA_DIR/quast/$SAMPLE_NAME"

    if [[ ! -f "$READ1" || ! -f "$READ2" || ! -f "$BWA_BAM" ]]; then
        echo "Skipping $SAMPLE_NAME: Missing required files"
        continue
    fi

    echo "Running QUAST for $SAMPLE_NAME"
    mkdir -p "$QUAST_OUT"
    quast.py \
        "$FNA" \
        -g "$PROKKA_OUT/${SAMPLE_NAME}.gff" \
        -1 "$READ1" \
        -2 "$READ2" \
        --bam "$BWA_BAM" \
        --threads "$THREADS" \
        --output-dir "$QUAST_OUT"
done

# Combine QUAST results
python -c "
import pandas as pd
import glob
files = glob.glob('$DATA_DIR/quast/*/transposed_report.tsv')
combined = pd.concat([pd.read_table(f) for f in files], axis=0)
combined.to_csv('$DATA_DIR/quast/combined_report.tsv', sep='\t', index=False)
"

# CheckM
CHECKM_OUT=$DATA_DIR/checkm
checkm taxonomy_wf \
    domain Bacteria \
    "$DATA_DIR/assembly" \
    "$CHECKM_OUT/tax_wf_outs" \
    -x fna \
    --tab_table \
    -f "$CHECKM_OUT/tax_wf_outs/tax_wf_result.tsv" \
    -t "$THREADS"

checkm lineage_wf \
    "$DATA_DIR/assembly" \
    "$CHECKM_OUT/lin_wf_outs" \
    -x fna \
    --tab_table \
    -f "$CHECKM_OUT/lin_wf_outs/lin_wf_result.tsv" \
    -t "$THREADS" \
    --pplacer_threads 2

checkm nx_plot -x fna "$DATA_DIR/assembly" "$CHECKM_OUT/plots"
checkm len_hist -x fna "$DATA_DIR/assembly" "$CHECKM_OUT/plots"
checkm marker_plot -x fna "$CHECKM_OUT/lin_wf_outs" "$DATA_DIR/assembly" "$CHECKM_OUT/plots"
```

## GTDB-Tk taxonomy assignment (for isolates only)

In this and the following two sections, we use three methods to assess the taxonomic profile of our sequencing data (`gtdbtk`, `metaphlan`, and `kraken2`). `gtdbtk` only makes sense for isolates or binned MAGs.

### Get GTDB-Tk database from AWS S3

Let's say we want to store the database at `~/data/gtdbtk` and we are using release 220.

```bash
mkdir -p ~/data/gtdbtk
aws s3 cp s3://yiming-qu/data/gtdbtk/release220.tar.gz ~/data/gtdbtk
aws s3 cp s3://yiming-qu/data/gtdbtk/release220_sketch.msh ~/data/gtdbtk
tar xzf ~/data/gtdbtk/release220.tar.gz
```

If you want to get the data from source, simply download the tar from GTDB-Tk documentation.

It's okay if the mash file is not available, the following `gtdbtk classify_wf` command will compute the mash and store at `--mash_db` if it does not exist.

```bash
DATA_DIR=<data_dir>
GTDBTK_DATA_PATH=$HOME/data/gtdbtk/release220
GTDBTK_MASH=$HOME/data/gtdbtk/release220_sketch.msh
gtdbtk classify_wf --genome_dir $DATA_DIR/assembly \
    --out_dir $DATA_DIR/gtdbtk_classify_wf \
    -x fna \
    --cpus $THREADS \
    --keep_intermediates \
    --mash_db $GTDBTK_MASH
```

## Metaphlan taxonomy profile

By default, when you run it for the first time, metaphlan will download the latest database if it does not exist. So we don't need separate commands to get the database.

### Run metaphlan

```bash
#!/usr/bin/env bash

DATA_DIR=<your_data_dir>
# Number of threads for parallel operations
THREADS=16

mkdir -p "$DATA_DIR/metaphlan"

for READ1 in "$DATA_DIR/fastp/"*_1.fq.gz; do
    READ2="${READ1%_1.fq.gz}_2.fq.gz"
    BASENAME=$(basename "$READ1")
    SAMPLE="${BASENAME%_1.fq.gz}"

    Run MetaPhlAn4 on raw reads (FASTQ)
    metaphlan \
        "${READ1},${READ2}" \
        --input_type fastq \
        --bowtie2out "$DATA_DIR/metaphlan/${SAMPLE}_metaphlan_bowtie2_output.bz2" \
        --nproc "$THREADS" \
        -o "$DATA_DIR/metaphlan/${SAMPLE}_metaphlan_profile.txt"
done
```

## kraken2 taxonomy profile

> To run kraken2, the entire database needs to be stored in memory. Therefore, you need at least ~90GB memory to run kraken2.

### Prepare kraken2 database

Kraken database download and processing and bracken database construction take a long time (several hours), large disk  space (hundreds of GBs to 1TB), and large memory (~90GB). Therefore, we try to reuse previously prepared database, though the database will not be perfectly up to date.

Very likely the database is already on AWS S3, but in case it's horribly outdated or you are interested in trying it out yourself, follow these steps (assuming database is stored in `~/data/kraken2_db/standard`):

1. Download and build kraken standard database

```bash
kraken2-build --standard --db ~/data/kraken2_db/standard --threads 16
```

2. Build bracken database from kraken database

```bash
bracken-build -d ~/data/kraken2_db/standard -l 150 -k 25 -t 8
```

3. Upload full database to AWS S3

````bash
aws s3 sync ~/data/kraken2_db/standard s3://yiming-qu/data/kraken2_db/standard
````

4. Clean kraken database

```bash
kraken2-build --clean --db ~/data/kraken2_db/standard
```

5. Upload cleaned database to AWS S3

```bash
aws s3 sync ~/data/kraken2_db/standard s3://yiming-qu/data/kraken2_db/standard_clean
```

Here is the structure of the database folder:

```bash
-rw-rw-r-- 1 ubuntu ubuntu  78G Dec 21 21:03 database.kraken
-rw-rw-r-- 1 ubuntu ubuntu 4.1M Dec 21 21:51 database150mers.kmer_distrib
-rw-rw-r-- 1 ubuntu ubuntu  21M Dec 21 21:51 database150mers.kraken
-rw-rw-r-- 1 ubuntu ubuntu  84G Dec 15 08:40 hash.k2d
-rw-rw-r-- 1 ubuntu ubuntu   64 Dec 15 08:40 opts.k2d
-rw-rw-r-- 1 ubuntu ubuntu 4.0M Dec 15 04:38 taxo.k2d
-rw-rw-r-- 1 ubuntu ubuntu  859 Dec 15 04:18 unmapped.txt
```

where the `database*` are bracken database and others are kraken database.

### Start from existing kraken2 database

#### Get database from AWS S3

```bash
aws s3 sync s3://yiming-qu/data/kraken2_db/standard_clean ~/data/kraken2_db/standard_clean
```

### Copy kraken database to shared memory to speed up calculation

Assuming kraken2 database is in `~/data/kraken2_db/standard`

```bash
sudo sudo mount -o remount,size=90G /dev/shm
sudo mkdir -p /dev/shm/kraken2_db
sudo mkdir -p /dev/shm/kraken2_db/standard_clean
sudo cp ~/data/kraken2_db/standard_clean/hash.k2d /dev/shm/kraken2_db/standard_clean/
sudo cp ~/data/kraken2_db/standard_clean/taxo.k2d /dev/shm/kraken2_db/standard_clean/
sudo cp ~/data/kraken2_db/standard_clean/opts.k2d /dev/shm/kraken2_db/standard_clean/
sudo cp ~/data/kraken2_db/standard_clean/unmapped.txt /dev/shm/kraken2_db/standard_clean/
```

### Run kraken2 and bracken

Assuming paired-end raw reads after quality control are in `$DATA_DIR/fastp` with names `<sample>_1.fq.gz` and `<sample>_2.fq.gz`

```bash
#!/usr/bin/env bash

DATA_DIR=<your_data_dir>
KRAKEN_DB=~/data/kraken2_db/standard_clean  # Kraken2 DB on disk
KRAKEN_DB_MEM=/dev/shm/kraken2_db/standard_clean  # Kraken2 DB in memory
# Number of threads for parallel operations
THREADS=16

mkdir -p "$DATA_DIR/kraken"

for READ1 in "$DATA_DIR/fastp/"*_1.fq.gz; do
    READ2="${READ1%_1.fq.gz}_2.fq.gz"
    BASENAME=$(basename "$READ1")
    SAMPLE="${BASENAME%_1.fq.gz}"
    CONTIGS="$DATA_DIR/assembly/$SAMPLE.fna"
    BWA_BAM="$DATA_DIR/bwa/$SAMPLE/$SAMPLE.sorted.bam"

    # Run Kraken2 on raw reads
    kraken2 \
        --db "$KRAKEN_DB_MEM" \
        --paired \
        --threads "$THREADS" \
        --output "$DATA_DIR/kraken/${SAMPLE}_kraken_output.txt" \
        --report "$DATA_DIR/kraken/${SAMPLE}_kraken_report.txt" \
        --use-names \
        --memory-mapping \
        "$READ1" "$READ2"
    
    # Run Bracken on the Kraken2 report
    bracken \
        -d "$KRAKEN_DB" \
        -i "$DATA_DIR/kraken/${SAMPLE}_kraken_report.txt" \
        -o "$DATA_DIR/kraken/${SAMPLE}_bracken_report.txt" \
        -r 150 \
        -l S
done
```

### Release database from memory

```bash
sudo rm -rf /dev/shm/kraken2_db
```

### Important output files

Coming next.

## Metagenomic binning

Metagenomic binning groups metagenome-assembled genomes based on cooccurrence such that each group ideally represents one organism. We can also use metagenomic binning on isolate WGS as a decontamination method. 

```bash
#!/usr/bin/env bash

DATA_DIR=<your_data_dir>
GTDBTK_DATA_PATH=$HOME/data/gtdbtk/release220
GTDBTK_MASH=$HOME/data/gtdbtk/release220_sketch.msh

# Number of threads for parallel operations
THREADS=16

for READ1 in "$DATA_DIR/fastp/"*_1.fq.gz; do
    READ2="${READ1%_1.fq.gz}_2.fq.gz"
    BASENAME=$(basename "$READ1")
    SAMPLE="${BASENAME%_1.fq.gz}"
    CONTIGS="$DATA_DIR/assembly/$SAMPLE.fna"
    BWA_BAM="$DATA_DIR/bwa/$SAMPLE/$SAMPLE.sorted.bam"

    echo "Running MetaBAT2 for binning: $$SAMPLE"
    METABAT_OUT="$DATA_DIR/metabat/$SAMPLE"
    mkdir -p "$METABAT_OUT"
    jgi_summarize_bam_contig_depths --outputDepth "$METABAT_OUT/depth.txt" $BWA_BAM
    metabat2 \
        -i "$CONTIGS" \
        -a "$METABAT_OUT/depth.txt" \
        -o "$METABAT_OUT/$SAMPLE-bin" \
        --unbinned
    # rename binned fasta files
    cd $METABAT_OUT
    mkdir others
    mv $SAMPLE-bin.lowDepth.fa others
    mv $SAMPLE-bin.unbinned.fa others
    mv $SAMPLE-bin.tooShort.fa others
    for file in $SAMPLE-bin.*.fa; do
        # Replace the first dot with an underscore
        new_name=$(echo "$file" | sed 's/\./-/')
        mv "$file" "$new_name"
    done
    cd -

    META_GTDBTK_OUT=$DATA_DIR/meta_gtdbtk_classify_wf/$SAMPLE
    gtdbtk classify_wf \
        --genome_dir $METABAT_OUT \
        --out_dir $META_GTDBTK_OUT \
        -x fa \
        --cpus $THREADS \
        --keep_intermediates \
        --mash_db $GTDBTK_MASH

    META_CHECKM_OUT=$DATA_DIR/meta_checkm/$SAMPLE
    checkm taxonomy_wf \
        domain Bacteria \
        $METABAT_OUT \
        $META_CHECKM_OUT/tax_wf_outs \
        -x fa \
        --tab_table \
        -f $META_CHECKM_OUT/tax_wf_outs/tax_wf_result.tsv \
        -t $THREADS
    checkm lineage_wf \
        $METABAT_OUT \
        $META_CHECKM_OUT/lin_wf_outs \
        -x fa \
        --tab_table \
        -f $META_CHECKM_OUT/lin_wf_outs/lin_wf_result.tsv \
        -t $THREADS \
        --pplacer_threads $THREADS
    checkm nx_plot -x fa $METABAT_OUT $META_CHECKM_OUT/plots
    checkm len_hist -x fa $METABAT_OUT $META_CHECKM_OUT/plots
    checkm marker_plot -x fa $META_CHECKM_OUT/lin_wf_outs $METABAT_OUT $META_CHECKM_OUT/plots
done
```

### Remove contamination for isolate sequencing

If your samples are isolates, ideally each sample should have only one contig, or all contigs should come from the same organism, but there can be contamination. Especially when sequencing depth is adequate, contamination organisms will also get assembled. Plus there can be viruses or plasmids that we don't care about but have a very high copy number.

We can essentially treat these samples as a community, bin the contigs, and take the bin that accounts for the most sequencing depth, and treat that as the cleaned isolate assembly.

This is implemented in the first section of `s3://yiming-qu/20241206_llin_isolate_assembly/analysis_coverage.ipynb`. Run it to get cleaned isolate assemblies. Practically we aggregate bins if they have the same GTDB-Tk taxonomy, and take the aggregated bin with highest sequencing depth.

### Rerun some analysis after decontamination

Since we now only work with a subset of all the contigs we got earlier, we need to rerun some analysis on the cleaned contigs. All these steps here have been introduced previously.

Theoretically for things like annotation, we can just extract from previous results for the cleaned contigs, but for simplicity and in case that some software relies on global information, we rerun the analysis.

```bash
#!/usr/bin/env bash

# Number of threads for parallel operations
THREADS=16

# Directories
DATA_DIR=<your_data_dir>  # adjust as needed
GTDBTK_DATA_PATH="$HOME/data/gtdbtk/release220"
GTDBTK_MASH="$HOME/data/gtdbtk/release220_sketch.msh"

# Create output directories
mkdir -p "$DATA_DIR/clean_prokka"
mkdir -p "$DATA_DIR/clean_bwa"
mkdir -p "$DATA_DIR/clean_quast"
mkdir -p "$DATA_DIR/clean_checkm"
mkdir -p "$DATA_DIR/clean_gtdbtk_classify_wf"

# Loop through each .fna in clean_assembly
for CLEAN_FASTA in "$DATA_DIR"/clean_assembly/*.fna; do
    if [[ ! -f "$CLEAN_FASTA" ]]; then
        echo "Warning: No .fna files found in $DATA_DIR/clean_assembly. Exiting..."
        exit 1
    fi

    BASENAME="$(basename "$CLEAN_FASTA")"
    SAMPLE_NAME="${BASENAME%.fna}"

    READ1="$DATA_DIR/fastp/${SAMPLE_NAME}_1.fq.gz"
    READ2="$DATA_DIR/fastp/${SAMPLE_NAME}_2.fq.gz"

    # Check for fastq reads
    if [[ ! -f "$READ1" ]]; then
        echo "Warning: $READ1 not found for sample $SAMPLE_NAME. Skipping..."
        continue
    fi
    if [[ ! -f "$READ2" ]]; then
        echo "Warning: $READ2 not found for sample $SAMPLE_NAME. Skipping..."
        continue
    fi

    echo "Processing sample: $SAMPLE_NAME"
    echo "Using cleaned assembly: $CLEAN_FASTA"

    # -------------------------------------------------------------------------
    # 1. Build BWA-MEM2 index
    # -------------------------------------------------------------------------
    echo "Indexing assembly for BWA-MEM2..."
    bwa-mem2 index "$CLEAN_FASTA"

    # -------------------------------------------------------------------------
    # 2. Align reads and generate sorted BAM
    # -------------------------------------------------------------------------
    CLEAN_BWA_OUT="$DATA_DIR/clean_bwa/${SAMPLE_NAME}"
    mkdir -p "$CLEAN_BWA_OUT"

    echo "Aligning reads..."
    bwa-mem2 mem -t "$THREADS" "$CLEAN_FASTA" "$READ1" "$READ2" \
        | samtools view -b -@ "$THREADS" -o "$CLEAN_BWA_OUT/${SAMPLE_NAME}.bam" -

    sambamba sort "$CLEAN_BWA_OUT/${SAMPLE_NAME}.bam" -t "$THREADS"

    # -------------------------------------------------------------------------
    # 3. Run Prokka
    # -------------------------------------------------------------------------
    CLEAN_PROKKA_OUT="$DATA_DIR/clean_prokka/${SAMPLE_NAME}"
    mkdir -p "$CLEAN_PROKKA_OUT"

    echo "Running Prokka..."
    prokka \
        --rnammer \
        --rfam \
        --cpus "$THREADS" \
        --outdir "$CLEAN_PROKKA_OUT" \
        --prefix "$SAMPLE_NAME" \
        --force \
        --centre X --compliant \
        "$CLEAN_FASTA"

    # -------------------------------------------------------------------------
    # 4. Run QUAST
    # -------------------------------------------------------------------------
    CLEAN_QUAST_OUT="$DATA_DIR/clean_quast/${SAMPLE_NAME}"
    mkdir -p "$CLEAN_QUAST_OUT"

    echo "Running QUAST..."
    quast.py \
        "$CLEAN_FASTA" \
        -g "$CLEAN_PROKKA_OUT/${SAMPLE_NAME}.gff" \
        -1 "$READ1" \
        -2 "$READ2" \
        --bam "$CLEAN_BWA_OUT/${SAMPLE_NAME}.sorted.bam" \
        --threads "$THREADS" \
        --output-dir "$CLEAN_QUAST_OUT"

    echo "QUAST completed for $SAMPLE_NAME."
done

# Combine QUAST reports
echo "Combining QUAST outputs..."
python -c "
import glob
import pandas as pd

files = glob.glob('$DATA_DIR/clean_quast/*/transposed_report.tsv')
if not files:
    print('Warning: No QUAST transposed_report.tsv files found to combine.')
else:
    combined = pd.concat([pd.read_table(f) for f in files], axis=0)
    combined.to_csv('$DATA_DIR/clean_quast/combined_report.tsv', sep='\t', index=False)
"

# -------------------------------------------------------------------------
# 5. Run GTDB-Tk classification
# -------------------------------------------------------------------------
echo "Running GTDB-Tk classify_wf..."
gtdbtk classify_wf \
    --genome_dir "$DATA_DIR/clean_assembly" \
    --out_dir "$DATA_DIR/clean_gtdbtk_classify_wf" \
    -x fna \
    --cpus "$THREADS" \
    --keep_intermediates \
    --mash_db "$GTDBTK_MASH"

# -------------------------------------------------------------------------
# 6. Run CheckM
# -------------------------------------------------------------------------
echo "Running CheckM..."
CLEAN_CHECKM_OUT="$DATA_DIR/clean_checkm"

mkdir -p "$CLEAN_CHECKM_OUT/tax_wf_outs"
mkdir -p "$CLEAN_CHECKM_OUT/lin_wf_outs"
mkdir -p "$CLEAN_CHECKM_OUT/plots"

# Taxonomy workflow
checkm taxonomy_wf domain Bacteria \
    "$DATA_DIR/clean_assembly" \
    "$CLEAN_CHECKM_OUT/tax_wf_outs" \
    -x fna \
    --tab_table \
    -f "$CLEAN_CHECKM_OUT/tax_wf_outs/tax_wf_result.tsv" \
    -t "$THREADS"

# Lineage workflow
checkm lineage_wf \
    "$DATA_DIR/clean_assembly" \
    "$CLEAN_CHECKM_OUT/lin_wf_outs" \
    -x fna \
    --tab_table \
    -f "$CLEAN_CHECKM_OUT/lin_wf_outs/lin_wf_result.tsv" \
    -t "$THREADS" \
    --pplacer_threads 2

# Additional CheckM plots
checkm nx_plot -x fna "$DATA_DIR/clean_assembly" "$CLEAN_CHECKM_OUT/plots"
checkm len_hist -x fna "$DATA_DIR/clean_assembly" "$CLEAN_CHECKM_OUT/plots"
checkm marker_plot -x fna "$CLEAN_CHECKM_OUT/lin_wf_outs" \
    "$DATA_DIR/clean_assembly" "$CLEAN_CHECKM_OUT/plots"

echo "All workflows completed for cleaned assemblies."

```

## DRAM protein functional profile

DRAM works on assembled genomes.

### Prepare DRAM database

#### Fix an error in the source code.

In version 1.4.6, there is an error in `mag_annotator/database_processing.py` around line 310 in the function `process_vogdb`:

```python
def process_vogdb(vog_hmm_targz, output_dir='.', logger=LOGGER, version=DEFAULT_VOGDB_VERSION, threads=1, verbose=True):
    hmm_dir = path.join(output_dir, 'vogdb_hmms')
    mkdir(hmm_dir)
    vogdb_targz = tarfile.open(vog_hmm_targz)
    vogdb_targz.extractall(hmm_dir)
    vog_hmms = path.join(output_dir, f'vog_{version}_hmms.txt')
    merge_files(glob(path.join(hmm_dir, 'VOG*.hmm')), vog_hmms)
    run_process(['hmmpress', '-f', vog_hmms], logger, verbose=verbose) 
    LOGGER.info('VOGdb database processed')
    return {'vogdb': vog_hmms}
```

It needs to be modified like:

```python
def process_vogdb(vog_hmm_targz, output_dir='.', logger=LOGGER, version=DEFAULT_VOGDB_VERSION, threads=1, verbose=True):
    hmm_dir = path.join(output_dir, 'vogdb_hmms')
    mkdir(hmm_dir)
    vogdb_targz = tarfile.open(vog_hmm_targz)
    vogdb_targz.extractall(hmm_dir)
    vog_hmms = path.join(output_dir, f'vog_{version}_hmms.txt')
    merge_files(glob(path.join(hmm_dir, "hmm", 'VOG*.hmm')), vog_hmms)
    run_process(['hmmpress', '-f', vog_hmms], logger, verbose=verbose) 
    LOGGER.info('VOGdb database processed')
    return {'vogdb': vog_hmms}
```

Constructing sqlite database is an extremely inefficient step. It takes ~500GB RAM at peak and runs on single core. On a r5a.16xlarge AWS EC2 instance, it took 4~5 hours.

```bash
DRAM-setup.py prepare_databases --output_dir ~/data/dram --threads 16
```

The temporary directory `database_files` takes about 70GB, which can be safely deleted. `tmp` is very small. The remaining essential files are about 900GB.

#### Upload to AWS S3

Run `sync` command twice like this:

```bash
aws s3 sync <DIR_OF_DRAM_DB> s3://yiming-qu/data/dram \
    --exclude "*/*" \
    --exclude "*mmsdb.idx*" \
    --exclude "pfam.mmsmsa" \
    --exclude "*.mmsdb" \
    --delete
```

```bash
aws s3 sync <DIR_OF_DRAM_DB> s3://yiming-qu/data/dram --exclude "database_files/*" --exclude "tmp/*" --delete
```

In this way, most large but unnecessary files will not be uploaded.

### Start from existing DRAM database

```bash
aws s3 sync s3://yiming-qu/data/dram <DIR_OF_DRAM_DB> --exclude "database_files/*" --exclude "tmp/*"
```

### Run DRAM

We run DRAM on all fasta assemblies with

```bash
DRAM.py annotate -i "assembly/*.fna" -o ./dram_annot --threads 16 
```

assuming your fasta files are `assembly/*.fna`.

TODO: How long does the program take?

## HUMAnN functional annotation

HUMAnN works on unassembled reads.

### Download HUMAnN database

```bash
humann_databases --download uniref uniref90_diamond ~/data/humann/uniref90_diamond
humann_databases --download chocophlan full  ~/data/humann/chocophlan
humann_config --update database_folders chocophlan ~/data/humann/chocophlan
```

### Run HUMAnN

```bash
DATA_DIR=<your_data_dir>
# Number of threads for parallel operations
THREADS=16

mkdir -p "$DATA_DIR/fastp_concat" $DATA_DIR/humann

for READ1 in "$DATA_DIR/fastp/"*_1.fq.gz; do
    READ2="${READ1%_1.fq.gz}_2.fq.gz"
    BASENAME=$(basename "$READ1")
    SAMPLE="${BASENAME%_1.fq.gz}"

    # concat read1 and read2
    zcat $READ1 $READ2 | pigz > "$DATA_DIR/fastp_concat/$SAMPLE.fq.gz"
    humann --input "$DATA_DIR/fastp_concat/$SAMPLE.fq.gz" --output $DATA_DIR/humann --threads $THREADS --remove-temp-output
    humann_renorm_table --input $DATA_DIR/humann/$SAMPLE_pathabundance.tsv --output $DATA_DIR/humann/$SAMPLE_pathabundance_norm.tsv --units cpm
    humann_renorm_table --input $DATA_DIR/humann/$SAMPLE_genefamilies.tsv --output $DATA_DIR/humann/$SAMPLE_genefamilies_norm.tsv --units cpm
done
```

## Metabolomic modeling

Coming next.
