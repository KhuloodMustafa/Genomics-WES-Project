# Whole-Exome Sequencing (WES) Analysis Pipeline for Osteopetrosis Variant Detection

This repository contains a pipeline for processing whole-exome sequencing (WES) data from a family trio (proband, father, mother) to identify candidate variants associated with osteopetrosis, a rare genetic bone disorder. The pipeline starts with raw FASTQ files, performs alignment, variant calling, annotation, and filtering, and generates a report of potential causative variants under the autosomal recessive inheritance model.

## Overview

The pipeline consists of 13 steps, implemented as Bash scripts, to process WES data for a trio where the proband is affected by osteopetrosis, and the consanguineous parents are unaffected. The goal is to identify rare, high-impact variants in known osteopetrosis genes.

## Prerequisites

- **System**: Linux (tested on Ubuntu)
- **Conda**: For managing environments and dependencies
- **Docker**: For running GATK tools
- **Reference Files**:
  - Human genome: `hg19.fa` (e.g., `/media/m0hamed/files/files/bioinformatics/wes/ref/hg19.fa`)
  - VEP cache: GRCh37 (e.g., `/media/m0hamed/files/files/bioinformatics/wes/ref/vep_cache`)
  - Exome interval file: intervals.list file for the the chromosome interval
- **Tools** (installed via Conda or Docker):
  - FastQC
  - MultiQC
  - Trimmomatic
  - BWA
  - Samtools
  - Picard
  - GATK (via Docker: `broadinstitute/gatk:latest`)
  - VEP
  - bgzip/tabix
  - bcftools
- **Conda Environments**:
  - `wes_analysis`: For FastQC, MultiQC, Trimmomatic, BWA, Samtools, Picard, VEP, bcftools
  ```bash
  conda create -n wes_analysis -c bioconda fastqc multiqc trimmomatic bwa samtools picard vep bcftools
  ```

## Pipeline Steps

The pipeline processes FASTQ files through alignment, variant calling, annotation, and variant filtering. Each step is implemented as a Bash script, with quality control performed manually or scripted separately.

### Step 0: Index Reference Genome (`0_index.sh`)

- **Tools**: Samtools, Picard, BWA
- **Description**: Indexes the `hg19.fa` reference genome to generate `.fai`, `.dict`, and BWA index files (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`) required for alignment and variant calling.
- **Code**:
  ```bash
  #!/bin/bash
  # Downlowd the genome and extract it
  wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
  gunzip hg19.fa.gz

  # Set input FASTA file

  FASTA="hg19.fa"

  # Step 1: Create .fai index with samtools
  echo "Indexing FASTA with samtools..."
  conda activate samtools
  samtools faidx "$FASTA"

  # Step 2: Create sequence dictionary (.dict) with Picard
  echo "Creating sequence dictionary with Picard..."
  conda activate picard
  picard CreateSequenceDictionary R="$FASTA" O="${FASTA%.fa}.dict"

  # Step 3: Index with BWA to generate .amb, .ann, .bwt, .pac, .sa
  echo "Indexing FASTA with BWA..."
  conda activate bwa
  bwa index "$FASTA"

  echo "Reference genome indexing complete!"
  ```
- **Explanation**: Generates necessary index files for `hg19.fa` using `samtools faidx` (`.fai`), Picard (`CreateSequenceDictionary` for `.dict`), and `bwa index` (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`).
- **Output**: `hg19.fa.fai`, `hg19.dict`, `hg19.fa.amb`, `hg19.fa.ann`, `hg19.fa.bwt`, `hg19.fa.pac`, `hg19.fa.sa`
- **Command**:
  ```bash
  ./0_index.sh
  ```

### Step 1: Download Data (`1_data.sh`)

- **Description**: Downloads raw FASTQ files for the trio from Zenodo.
- **Code**:
  ```bash
  #!/bin/bash
  wget https://zenodo.org/record/3243160/files/father_R1.fq.gz
  wget https://zenodo.org/record/3243160/files/father_R2.fq.gz
  wget https://zenodo.org/record/3243160/files/mother_R1.fq.gz
  wget https://zenodo.org/record/3243160/files/mother_R2.fq.gz
  wget https://zenodo.org/record/3243160/files/proband_R1.fq.gz
  wget https://zenodo.org/record/3243160/files/proband_R2.fq.gz
  ```
- **Explanation**: Uses `wget` to fetch paired-end FASTQ files for each sample, storing them in the working directory.
- **Output**: `father_R1.fq.gz`, `father_R2.fq.gz`, `mother_R1.fq.gz`, `mother_R2.fq.gz`, `proband_R1.fq.gz`, `proband_R2.fq.gz`
- **Command**:
  ```bash
  ./1_data.sh
  ```

### Step 2: Quality Control

- **Tools**: FastQC, MultiQC
- **Description**: Assesses the quality of raw FASTQ files for the trio (`father_R1.fq.gz`, `father_R2.fq.gz`, `mother_R1.fq.gz`, `mother_R2.fq.gz`, `proband_R1.fq.gz`, `proband_R2.fq.gz`) to identify issues like low-quality bases or adapter contamination.
- **Command**:
  ```bash
  mkdir fastqc_results
  fastqc -o fastqc_results *.fq.gz
  multiqc fastqc_results -o multiqc_report
  ```
- **Output**: `fastqc_results/*.html`, `multiqc_report/multiqc_report.html`

### Step 3: Trim Reads (`2_trim.sh`)

- **Tool**: Trimmomatic
- **Description**: Removes low-quality bases and adapters from FASTQ files to improve alignment quality.
- **Code**:
  ```bash
  #!/bin/bash
  samples=(father mother proband)
  total=${#samples[@]}
  for sample in "${samples[@]}"; do
    echo "Processing $sample"
    trimmomatic PE -phred33 \
      ${sample}_R1.fq.gz ${sample}_R2.fq.gz \
      ${sample}_R1_trimmed.fq.gz ${sample}_R1_unpaired.fq.gz \
      ${sample}_R2_trimmed.fq.gz ${sample}_R2_unpaired.fq.gz \
      LEADING:32 TRAILING:32 MINLEN:36 2>&1 | tee ${sample}_trimming.log
    echo "Finished $sample"
  done
  echo "All samples trimmed!"
  ```
- **Explanation**: Loops through samples, running Trimmomatic in paired-end mode to trim reads with quality <32 and length <36, logging results.
- **Output**: `${sample}_R1_trimmed.fq.gz`, `${sample}_R2_trimmed.fq.gz`, `${sample}_R1_unpaired.fq.gz`, `${sample}_R2_unpaired.fq.gz`, `${sample}_trimming_LOG`
- **Command**:
  ```bash
  ./2_trim.sh
  ```

### Step 4: Align Reads (`3_align.sh`)

- **Tool**: BWA
- **Description**: Aligns trimmed reads to the hg19 reference genome, producing SAM files.
- **Code**:
  ```bash
  #!/bin/bash
  samples=(father mother proband)
  reference=/media/m0hamed/files/files/bioinformatics/wes/ref/hg19.fa
  threads=1
  for sample in "${samples[@]}"; do
    echo "Aligning $sample to hg19..."
    bwa mem -t $threads -R "@RG\tID:${sample}\tSM:${sample}" \
      $reference \
      ${sample}_R1_trimmed.fq.gz ${sample}_R2_trimmed.fq.gz > ${sample}.sam
    echo "Finished aligning $sample"
  done
  echo "Alignment complete for all samples!"
  ```
- **Explanation**: Uses BWA-MEM to align paired-end reads, adding read group metadata (`-R`) and outputting SAM files.
- **Output**: `${sample}.sam`
- **Command**:
  ```bash
  ./3_align.sh
  ```

### Step 5: Filter Reads (`4_filter.sh`)

- **Tool**: Samtools
- **Description**: Filters SAM files to retain paired, mapped reads with high mapping quality (`-q 20`), converting to BAM.
- **Code**:
  ```bash
  #!/bin/bash
  samples=(father mother proband)
  for sample in "${samples[@]}"; do
    echo "Filtering paired reads for $sample..."
    samtools view -b -F 4 -F 8 -q 20 \
      ${sample}.sam > ${sample}_filtered.bam
    echo "Finished filtering $sample"
  done
  echo "Filtering complete for all samples!"
  ```
- **Explanation**: Uses `samtools view` to filter out unmapped (`-F 4`), unpaired (`-F 8`), and low-quality (`-q 20`) reads, producing BAM files.
- **Output**: `${sample}_filtered.bam`
- **Command**:
  ```bash
  ./4_filter.sh
  ```

### Step 6: Sort BAM Files (`5_sort.sh`)

- **Tool**: Samtools
- **Description**: Sorts and indexes filtered BAM files for downstream processing.
- **Code**:
  ```bash
  for sample in father mother proband; do
    samtools sort ${sample}_filtered.bam -o ${sample}_filtered_sorted.bam
    samtools index ${sample}_filtered_sorted.bam
  done
  ```
- **Explanation**: Sorts BAM files by coordinate using `samtools sort` and creates index files with `samtools index`.
- **Output**: `${sample}_filtered_sorted.bam`, `${sample}_filtered_sorted.bam.bai`
- **Command**:
  ```bash
  ./5_sort.sh
  ```

### Step 7: Mark Duplicates (`6_dedup.sh`)

- **Tool**: Picard
- **Description**: Marks and removes duplicate reads from sorted BAM files to reduce PCR artifacts.
- **Code**:
  ```bash
  #!/bin/bash
  samples=(father mother proband)
  for sample in "${samples[@]}"; do
    echo "Processing $sample..."
    echo "  Marking duplicates for $sample..."
    picard MarkDuplicates \
      INPUT=${sample}_filtered_sorted.bam \
      OUTPUT=${sample}_dedup.bam \
      METRICS_FILE=${sample}_dedup_metrics.txt \
      REMOVE_DUPLICATES=true \
      CREATE_INDEX=true
    echo "Finished processing $sample"
  done
  echo "Duplicate marking complete for all samples!"
  ```
- **Explanation**: Uses Picard’s `MarkDuplicates` to remove duplicates, generating deduplicated BAMs and metrics files.
- **Output**: `${sample}_dedup.bam`, `${sample}_dedup_metrics.txt`
- **Command**:
  ```bash
  ./6_dedup.sh
  ```

### Step 8: Index BAM Files (`7_index.sh`)

- **Tool**: Samtools
- **Description**: Creates index files for deduplicated BAMs for efficient access.
- **Code**:
  ```bash
  #!/bin/bash
  SAMPLES=(father mother proband)
  for SAMPLE in "${SAMPLES[@]}"; do
    echo "Indexing ${SAMPLE}_dedup.bam..."
    samtools index ${SAMPLE}_dedup.bam
    echo "Indexing complete for $SAMPLE"
  done
  ```
- **Explanation**: Indexes deduplicated BAM files using `samtools index`.
- **Output**: `${sample}_dedup.bam.bai`
- **Command**:
  ```bash
  ./7_index.sh
  ```

### Step 9: Variant Calling (`8_calling.sh`)

- **Tool**: GATK (HaplotypeCaller, GenomicsDBImport, GenotypeGVCFs)
- **Description**: Performs per-sample variant calling, combines gVCFs, and conducts joint genotyping, using an exome interval file to restrict analysis to targeted regions.
- **Code**:
  ```bash
#!/bin/bash

# Absolute path to your working directory
WORKDIR=$(pwd)
REF=$WORKDIR/hg19.fa
SAMPLES=(father)

# HaplotypeCaller for each sample
for SAMPLE in "${SAMPLES[@]}"; do
  echo "Running HaplotypeCaller for $SAMPLE..."
  docker run -v $WORKDIR:/data -w /data \
    broadinstitute/gatk:latest gatk --java-options "-Xmx6g" HaplotypeCaller \
    -R /data/hg19.fa \
    -I /data/${SAMPLE}_dedup.bam \
    -O /data/${SAMPLE}.g.vcf.gz \
    -ERC GVCF
  echo "HaplotypeCaller complete for $SAMPLE"
done

# Joint Genotyping
echo "Combining gVCFs..."
docker run -v $WORKDIR:/data -w /data \
  broadinstitute/gatk:latest gatk --java-options "-Xmx6g" GenomicsDBImport \
  -R /data/hg19.fa \
  -V /data/father.g.vcf.gz \
  -V /data/mother.g.vcf.gz \
  -V /data/proband.g.vcf.gz \
  -L /data/intervals.list \
  --genomicsdb-workspace-path /data/trio_db

echo "Running joint genotyping..."
docker run -v $WORKDIR:/data -w /data \
  broadinstitute/gatk:latest gatk GenotypeGVCFs \
  -R /data/hg19.fa \
  -V gendb:///data/trio_db \
  -O /data/trio.vcf.gz

echo "Variant calling complete!"
  ```
- **Explanation**: Converts the exome BED file to a Picard-style `.interval_list` using `BedToIntervalList`, then uses GATK’s `HaplotypeCaller` to generate per-sample gVCFs, combines them with `GenomicsDBImport`, and performs joint genotyping with `GenotypeGVCFs`, restricting all steps to exome regions using `-L intervals.list`.
- **Output**: `intervals.list`, `${sample}.g.vcf.gz`, `${sample}.g.vcf.gz.tbi`, `trio_db`, `trio.vcf.gz`, `trio.vcf.gz.tbi`
- **Command**:
  ```bash
  ./8_calling.sh
  ```

### Step 10: Post-Process Variants (`9_post_proc.sh`)

- **Tool**: GATK (LeftAlignAndTrimVariants)
- **Description**: Splits multiallelic variants and normalizes indels in the joint VCF.
- **Code**:
  ```bash
  #!/bin/bash
  WORKDIR=$(pwd)
  echo "Splitting multiallelics and normalizing indels..."
  docker run --memory-swap -1 -v $WORKDIR:/data -w /data \
    broadinstitute/gatk:latest gatk --java-options "-Xmx6g" LeftAlignAndTrimVariants \
    -R /data/hg19.fa \
    -V /data/trio.vcf.gz \
    -O /data/trio_normalized_split.vcf.gz \
    --split-multi-allelics
  echo "Normalization and splitting complete!"
  ```
- **Explanation**: Normalizes indels and splits multiallelic variants using GATK, improving VCF consistency.
- **Output**: `trio_normalized_split.vcf.gz`
- **Command**:
  ```bash
  ./9_post_proc.sh
  ```

### Step 11: Annotate Variants (`10_annotate_VEP.sh`)

- **Tool**: VEP
- **Description**: Annotates variants with functional consequences, ClinVar significance, and population frequencies.
- **Code**:
  ```bash
  #!/bin/bash
  WORKDIR=$(pwd)
  CACHEDIR=/media/m0hamed/files/files/bioinformatics/wes/ref/vep_cache
  echo "Running VEP annotation..."
  vep --cache \
    --dir_cache $CACHEDIR \
    --assembly GRCh37 \
    --input_file $WORKDIR/trio_normalized_split.vcf.gz \
    --output_file $WORKDIR/trio_annotated_vep.vcf \
    --vcf \
    --fork 1 \
    --everything \
    --verbose \
    --force_overwrite \
    --no_stats \
    --offline \
    || { echo "VEP failed! Check logs above."; exit 1; }
  echo "Compressing output..."
  bgzip -c $WORKDIR/trio_annotated_vep.vcf > $WORKDIR/trio_annotated_vep.vcf.gz
  tabix -p vcf $WORKDIR/trio_annotated_vep.vcf.gz
  echo "Annotation and compression complete!"
  ```
- **Explanation**: Runs VEP to add annotations (e.g., IMPACT, SYMBOL, CLIN_SIG) to the VCF, then compresses and indexes the output.
- **Output**: `trio_annotated_vep.vcf.gz`, `trio_annotated_vep.vcf.gz.tbi`
- **Command**:
  ```bash
  ./10_annotate_VEP.sh
  ```

### Step 12: Find Candidate Variants (`11_find.sh`)

- **Tool**: bcftools
- **Description**: Filters variants for autosomal recessive (proband homozygous, parents heterozygous) and compound heterozygous patterns, focusing on high/moderate impact variants.
- **Code**:
  ```bash
  #!/bin/bash
  VCF="/media/m0hamed/files/files/bioinformatics/wes/ref/trio_annotated_vep.vcf.gz"
  PROBAND="proband"
  FATHER="father"
  MOTHER="mother"
  AR_VCF="autosomal_recessive_candidates.vcf"
  AR_TSV="autosomal_recessive_candidates.tsv"
  CH_TEMP="compound_het_temp.vcf"
  CH_TSV="compound_het_candidates.tsv"
  echo "Running autosomal recessive filtering..."
  bcftools view \
      -s "$FATHER,$MOTHER,$PROBAND" \
      -i "GT[2]='1/1' && GT[0]='0/1' && GT[1]='0/1' && INFO/CSQ !~ 'IMPACT=LOW'" \
      -o "$AR_VCF" \
      "$VCF" || {
      echo "Autosomal recessive filtering failed!"
      exit 1
  }
  bcftools query \
      -f '%CHROM\t%POS\t%REF\t%ALT\t%CSQ\t%ID\tINFO/max_aaf_all\n' \
      "$AR_VCF" > "$AR_TSV"
  echo "Autosomal recessive candidates saved to $AR_TSV"
  echo "Top candidates:"
  head -n 10 "$AR_TSV" | column -t
  echo "Candidate variant detection complete!"
  ```
- **Explanation**: Uses `bcftools` to filter variants based on inheritance patterns, outputting VCF and TSV files.
- **Output**: `autosomal_recessive_candidates.vcf`, `autosomal_recessive_candidates.tsv`
- **Command**:
  ```bash
  ./11_find.sh
  ```

## Output Files

### `autosomal_recessive_candidates.tsv`
* **Source**: Generated by the `12_find.sh` script, which queries the GEMINI database (`trio.db`).
* **Purpose**: Contains candidate variants filtered based on an autosomal recessive inheritance pattern (e.g., likely homozygous in the proband and heterozygous in both unaffected parents) and pre-filtered for rarity using population allele frequency databases (e.g., `max_aaf_all < 0.01`, depending on the exact query used in `12_find.sh`).
* **Format**: Tab-separated values (TSV) with 7 columns.
* **Column Description**:
    1. `CHROM`: Chromosome name.
    2. `POS`: Genomic position (1-based).
    3. `REF`: Reference allele.
    4. `ALT`: Alternate allele(s).
    5. `VEP_Annotation`: The full annotation string generated by Variant Effect Predictor (VEP). This is a complex field containing multiple pieces of information separated by pipes (`|`), such as `Allele|Consequence|IMPACT|SYMBOL|Gene|...`. There may be multiple annotations per variant (separated by commas `,`), typically one per affected transcript. **Note**: Key information like predicted impact, gene symbol, and detailed allele frequencies must be parsed from this string.
    6. `Quality/Placeholder`: Often contains a placeholder like `.` (potentially the QUAL field from the original VCF if not otherwise specified by the query).
    7. `Filter_Info`: Contains the name of the allele frequency filter field used by the GEMINI query (e.g., `INFO/max_aaf_all`). **Note**: This column does *not* contain the actual allele frequency value for the variant; frequency information should be extracted from the VEP Annotation string (column 5) if needed or checked in external databases.
* **Interpretation**: This file provides a list of rare variants fitting the specified inheritance model. However, it requires further analysis to prioritize candidates. Subsequent filtering based on VEP, checking against known disease genes (like `TCIRG1`, `CLCN7`, `CA2` etc.), assessing biological relevance, and consulting databases like ClinVar and gnomAD are necessary next steps to narrow down the list to the most promising candidates. Our analysis showed this file contained 289 variants.

## File Structure

```
wes_pipeline/
├── 0_index.sh
├── 1_data.sh
├── 2_trim.sh
├── 3_align.sh
├── 4_filter.sh
├── 5_sort.sh
├── 6_dedup.sh
├── 7_index.sh
├── 8_calling.sh
├── 9_post_proc.sh
├── 10_annotate_VEP.sh
├── 11_find.sh
├── fastqc_results/
├── multiqc_report/
├── hg19.fa
├── hg19.fa.fai
├── hg19.dict
├── hg19.fa.amb
├── hg19.fa.ann
├── hg19.fa.bwt
├── hg19.fa.pac
├── hg19.fa.sa
├── exome_intervals.bed
├── intervals.list
├── father_R1.fq.gz
├── father_R2.fq.gz
├── mother_R1.fq.gz
├── mother_R2.fq.gz
├── proband_R1.fq.gz
├── proband_R2.fq.gz
├── *.sam
├── *.bam
├── *.vcf.gz
├── trio.ped
├── trio.db
├── autosomal_recessive_candidates.tsv
```

## Steps
1. Ensure `hg19.fa`, exome interval file (`intervals.list`), and VEP cache are in `/media/m0hamed/files/files/bioinformatics/wes/ref/`.
2. Run scripts sequentially:
   ```bash
   bash 0_index.sh
   bash 1_data.sh
   # Perform quality control with FastQC
   bash 2_trim.sh
   # ... continue through 12_find.sh
   ```
3. Perform quality control separately before Step 3.

## Notes

- **Sample Names**: Hardcoded as `father`, `mother`, `proband`. Adjust in scripts if different.
- **Reference Genome**: Uses `hg19.fa`. Update paths in scripts for your setup.
- **Interval File**: Requires an exome capture kit-specific BED file (e.g., `exome_intervals.bed`). Convert to `.interval_list` in `8_calling.sh`.
- **VEP Cache**: Assumes GRCh37 cache at `/media/m0hamed/files/files/bioinformatics/wes/ref/vep_cache`.
