# Variant Calling Workflow: Single-Sample GATK Mutect2 with Normalization
## Overview
This code snippet demonstrates a single-sample variant calling pipeline using GATK Mutect2, followed by filtering and normalization steps. 
It processes multiple tumor samples independently without using a panel of normals (PON), which can sometimes overfilter biologically relevant variants.

```bash
#!/bin/bash

# Single-sample calling with normalization (--panel-of-normals "$PON" removed 
# because can indeed lead to overfiltering and loss of biologically relevant variants)
echo " - Calling variants in each sample independently"
for sample in "${SAMPLES[@]}"; do
    echo " - Processing $sample"
    
    $gatk Mutect2 \
        -R "$REFERENCE" \
        -I "$OUTDIR/bam/${sample}.bqsr.bam" \
        -L "$BED" \
        --germline-resource "$GNOMAD" \
        -O "$OUTDIR/vcf/single_sample/${sample}.single.vcf.gz"
    
    $gatk FilterMutectCalls \
        -V "$OUTDIR/vcf/single_sample/${sample}.single.vcf.gz" \
        -R "$REFERENCE" \
        -L "$BED" \
        --contamination-table "$OUTDIR/metrics/contamination/${sample}.contamination.table" \
        -O "$OUTDIR/vcf/single_sample/${sample}.single.filtered.vcf.gz"
    
    # ADDED: Normalization
    bcftools norm -f "$REFERENCE" -m-both \
        "$OUTDIR/vcf/single_sample/${sample}.single.filtered.vcf.gz" \
        -O z -o "$OUTDIR/vcf/single_sample/${sample}.single.filtered.normalized.vcf.gz"
    
    tabix -p vcf "$OUTDIR/vcf/single_sample/${sample}.single.filtered.normalized.vcf.gz"

    echo "     Completed single-sample calling for $sample"
done
```

## Best Practices Notes

Quality Control:Ensure BAM files pass QC before variant calling <br>
Resource Files: Use appropriate version-matched reference files <br>
Memory: GATK Mutect2 may require significant memory for large genomes
