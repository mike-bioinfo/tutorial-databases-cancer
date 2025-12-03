## This is a basic pipeline for pre-processing prior to variant calling process which can be improved depending on what we are trying to answer

```bash
for sample in "${SAMPLES[@]}"; do
    echo "Processing $sample"
    
    # 1. Alignment
    echo " - Aligning with BWA-MEM"
    $bwa mem -t $THREADS \
    -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" \
    "$REFERENCE" \
    "${WES_R1[$sample]}" "${WES_R2[$sample]}" \
    2> "$OUTDIR/bam/${sample}.bwa.log" | \
    samtools view -@ $THREADS -Sb - | \
    samtools sort -@ $THREADS -o "$OUTDIR/bam/${sample}.sorted.bam"
    samtools index "$OUTDIR/bam/${sample}.sorted.bam"
    samtools flagstat "$OUTDIR/bam/${sample}.sorted.bam" > "$OUTDIR/bam/${sample}.mapped.metrics"
    $pandepth -i "$OUTDIR/bam/${sample}.sorted.bam" -o "$OUTDIR/bam/${sample}.pandepth" -b "$BED"
    gunzip -c "$OUTDIR/bam/${sample}.pandepth.bed.stat.gz" | (head -n 1 && tail -n +2 | sort -k9,9nr) > "$OUTDIR/bam/${sample}.pandepth.sorted.xls"

    # 2. Mark duplicates
    echo " - Marking duplicates"
    $gatk MarkDuplicates \
        -I "$OUTDIR/bam/${sample}.sorted.bam" \
        -O "$OUTDIR/bam/${sample}.dedup.bam" \
        -M "$OUTDIR/bam/${sample}.metrics.txt"
    
    # 3. Base quality recalibration (WITH GNOMAD + DBSNP)
    echo " - Base quality recalibration"
    $gatk BaseRecalibrator \
        -R "$REFERENCE" \
        -I "$OUTDIR/bam/${sample}.dedup.bam" \
        --known-sites "$GNOMAD" \
        --known-sites "$DBSNP" \
        -O "$OUTDIR/bam/${sample}.recal.table"

    echo " - Applying BQSR"
    $gatk ApplyBQSR \
        -R "$REFERENCE" \
        -I "$OUTDIR/bam/${sample}.dedup.bam" \
        --bqsr-recal-file "$OUTDIR/bam/${sample}.recal.table" \
        -O "$OUTDIR/bam/${sample}.bqsr.bam"

    samtools index "$OUTDIR/bam/${sample}.bqsr.bam"
    samtools flagstat "$OUTDIR/bam/${sample}.bqsr.bam" > "$OUTDIR/bam/${sample}.bqsr.metrics"
    $pandepth -i "$OUTDIR/bam/${sample}.bqsr.bam" -o "$OUTDIR/bam/${sample}.pandepth.bqsr" -b "$BED"
    gunzip -c "$OUTDIR/bam/${sample}.pandepth.bqsr.bed.stat.gz" | (head -n 1 && tail -n +2 | sort -k9,9nr) > "$OUTDIR/bam/${sample}.pandepth.bqsr.xls"
    
    echo "Completed $sample"
    echo "---"
done
```
