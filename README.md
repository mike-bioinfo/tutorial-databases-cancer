# Preparing Databases for Mapping and Variant Calling in Cancer Research

### (Whole Genome and Exome Workflows)

This guide describes how to prepare reference genomes, annotations, and
supporting resources for downstream analyses such as alignment, somatic
variant calling (e.g., Mutect2), germline calling, and structural
variant workflows.

## Software required

SEQKIT: https://bioinf.shenwei.me/seqkit/
PICARD: https://broadinstitute.github.io/picard/
SAMTOOLS: https://www.htslib.org/
BWA: https://github.com/lh3/bwa
GATK: https://github.com/broadinstitute/gatk/releases

## Download Human Reference Genome

``` bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh38.primary_assembly.genome.fa.gz
```

## Download Annotation Files (GTF & GFF3)

``` bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gff3.gz
```

## Install and Configure gsutil

Instructions: https://cloud.google.com/storage/docs/gsutil_install

### Fetch Panel of Normals (PON)

``` bash
gsutil -m cp gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz{,.tbi} .
```

### Fetch gnomAD Germline Resource (Mutect2)

``` bash
gsutil -m cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz{,.tbi} .
```

## Download Variant Databases

### dbSNP

``` bash
wget https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz
```

### ClinVar

``` bash
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
```

## Extract Major Chromosomes ID

``` bash
grep '>' GRCh38.primary_assembly.genome.fa | head -n25 > ids_chr.txt
```

## Extract Chromosomes Using SeqKit

``` bash
seqkit grep -n -f ids_chr.txt GRCh38.primary_assembly.genome.fa > GRCh38.primary_assembly.genome.final.fa
```

## Sort Genome by Chromosome Number

``` bash
seqkit sort -nN GRCh38.primary_assembly.genome.final.fa -o GRCh38.primary_sorted.genome.final.fa
```

## Clean FASTA Headers

``` bash
sed -i 's/ .*//' GRCh38.primary_sorted.genome.final.fa
```

## Index and Create Sequence Dictionary (GATK)

``` bash
samtools faidx GRCh38.primary_assembly.genome.final.fa

java -jar picard.jar CreateSequenceDictionary R=GRCh38.primary_assembly.genome.final.fa O=GRCh38.primary_assembly.genome.final.dict
```

## Prepare ID File for Annotation Filtering

``` bash
sed -i 's/ .*//' id_primary_ctg.txt
```

## Extract Annotation Subsets

``` bash
grep -w -f id_primary_ctg.txt gencode.v48.annotation.gtf > gencode.v48.annotation.modified.gtf
grep -w -f id_primary_ctg.txt gencode.v48.annotation.gff3 > gencode.v48.annotation.modified.gff3
```

## Restore Header to Modified GTF

``` bash
{ grep '^#' gencode.v48.annotation.gtf; cat gencode.v48.annotation.modified.gtf; } > tmp && mv tmp gencode.v48.annotation.modified.gtf
```

## Generate BED File from GTF

``` bash
awk '$3 == "gene"' gencode.v48.annotation.modified.gtf > gencode.v48.annotation_only_Genes.gtf
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' gencode.v48.annotation_only_Genes.gtf > gencode.v48.annotation_reformated.gtf
convert2bed -i gtf < gencode.v48.annotation_reformated.gtf > gencode.v48.annotation_reformated.bed
awk -F '[\t;]' '{print $1,$2,$3,$12}' gencode.v48.annotation_reformated.bed > gencode.v48.annotation_coords.bed
awk -F '[ ]' '{print $1,$2,$3,$6}' gencode.v48.annotation_coords.bed > gencode.v48.annotation_minimal.bed
sed 's/"//g' gencode.v48.annotation_minimal.bed | sed 's/ /\t/g' > GRCh38.annotation_final.bed
```

### One-line Alternative

``` bash
zcat gencode.v48.annotation.modified.gtf.gz | grep -P "\tgene\t" | cut -f1,4,5,7,9 | sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | awk -F $'\t' 'BEGIN { OFS=FS } { print $1,$2-1,$3,$6,".",$4,$10,$12,$14 }' | sort -k1,1 -k2,2n | cut -f1,2,3,7
```

## Index Reference Genome with BWA

``` bash
bwa index GRCh38.primary_assembly.genome.final.fa
```

## Changing the vcf annotation adding "chr" prefix to chromosome names in COSMIC vcf file (https://cancer.sanger.ac.uk/cosmic/download/cosmic)

``` bash
zcat Cosmic_CompleteTargetedScreensMutant_v102_GRCh38.vcf.gz | awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' > Cosmic_CompleteTargetedScreensMutant_v102_GRCh38_modified.vcf
zcat Cosmic_GenomeScreensMutant_v102_GRCh38.vcf.gz | awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' > Cosmic_GenomeScreensMutant_v102_GRCh38_modified.vcf
zcat Cosmic_NonCodingVariants_v102_GRCh38.vcf.gz | awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' > Cosmic_NonCodingVariants_v102_GRCh38_modified.vcf
```

## If it is necessary, you can also remove the "chr" from vcf file.
## Adding or removing prefix "chr" in vcf databases files are critical for a comprehensive annotation and interpretation of results.

``` bash
zcat Cosmic_CompleteTargetedScreensMutant_v102_GRCh38.vcf.gz | awk '{gsub(/^chr/,""); print}' > no_chr.vcf
```

## Notes

-   Please make sure in your vcf databases and annotations are following the same naming convention about chromosomes (chr or simply numbers, i.e chr1 or 1).
-   Ensure consistent genome versions across FASTA, GTF, GFF3, dbSNP,
    ClinVar, PON, and gnomAD.    
-   GATK requires `.fa`, `.fa.fai`, `.dict`, and indexed `.vcf.gz` +
    `.tbi`.
-   Always verify file checksums.
