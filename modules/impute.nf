#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
=                                 h3achipimputation                                    =
========================================================================================
 h3achipimputation imputation functions and processes.
----------------------------------------------------------------------------------------
 @Authors

----------------------------------------------------------------------------------------
 @Homepage / @Documentation
  https://github.com/h3abionet/chipimputation
----------------------------------------------------------------------------------------

----------------------------------------------------------------------------------------
*/

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

process impute_minimac4 {
    tag "imp_${target_name}_${chrm}:${chunk_start}-${chunk_end}_${ref_name}_${tagName}"
    label "verylarge"
    input:
        tuple val(chrm), val(chunk_start), val(chunk_end), val(target_name), file(target_phased_vcf), val(ref_name), file(ref_vcf), file(ref_m3vcf), val(tagName)
    output:
        tuple val(chrm), val(chunk_start), val(chunk_end), val(target_name), val(ref_name), file("${base}_imputed.dose.vcf.gz"), file("${base}_imputed.info"), val(tagName)
    shell:
        base = "${file(target_phased_vcf.baseName).baseName}"
        """
        nblines=\$(zcat ${target_phased_vcf} | grep -v "^#" | wc -l)
        if (( \$nblines > 1 ))
        then
            minimac4 \
                --refHaps ${ref_m3vcf} \
                --haps ${target_phased_vcf} \
                --format GT,DS \
                --allTypedSites \
                --minRatio ${params.minRatio} \
                --chr ${chrm} --start ${chunk_start} --end ${chunk_end} --window ${params.buffer_size} \
                --prefix ${base}_imputed \
                --cpus ${task.cpus}
        else
             touch ${base}_imputed.dose.vcf && bgzip ${base}_imputed.dose.vcf
             touch ${base}_imputed.info
        fi
        """
}

process impute_minimac4_1 {
    tag "imp_${target_name1}_${chrm1}:${chunk_start1}-${chunk_end1}_${ref_name1}_${tagName1}"
    label "bigmem"
    input:
        tuple val(chrm1), val(chunk_start1), val(chunk_end1), val(target_name1), file(target_phased_vcf1), val(ref_name1), file(ref_vcf1), file(ref_m3vcf1), val(tagName1)
    output:
        tuple val(chrm1), val(chunk_start1), val(chunk_end1), val(target_name1), val(ref_name1), file("${base}_imputed.dose.vcf.gz"), file("${base}_imputed.info"), val(tagName1)
    shell:
        base = "${file(target_phased_vcf1.baseName).baseName}"
        """
        nblines=\$(zcat ${target_phased_vcf1} | grep -v "^#" | wc -l)
        if (( \$nblines > 1 ))
        then
            bcftools annotate  --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' ${target_phased_vcf1} | 
            vcftools --vcf - --keep-INFO-all --max-missing 0.05 --hwe 0.00001 --mac 1 --recode --stdout | bgzip > ${base}.id.vcf.gz
            minimac4 \
                --cpus ${task.cpus} \
                --refHaps ${ref_m3vcf1} \
                --haps ${base}.id.vcf.gz \
                --format GT,DS \
                --allTypedSites \
                --minRatio ${params.minRatio} \
                --chr ${chrm1} --start ${chunk_start1} --end ${chunk_end1} --window ${params.buffer_size} \
                --prefix ${base}_imputed
        else
             touch ${base}_imputed.dose.vcf && bgzip ${base}_imputed.dose.vcf
             touch ${base}_imputed.info
        fi
        """
}