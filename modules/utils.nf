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

process vcf_to_m3vcf {
    tag "m3vcf_${dataset}_${chrm}"
    publishDir "${params.outdir}/${dataset}/m3vcfs", overwrite: true, mode:'copy'
    label "bigmem10"

    input:
        tuple val(dataset), val(chrm), file(dataset_vcf), file(dataset_vcf_idx)

    output:
        tuple val(dataset), val(chrm), file(dataset_m3vcf)

    script:
        base = file(dataset_vcf.baseName).baseName
        dataset_m3vcf = "${base}.m3vcf.gz"
        params_ = ''
        if(params.genome_build == 'b38'){ params_ = " --mychromosome ${chrm}"}
        """
        minimac3 \
            --refHaps ${dataset_vcf} \
            --processReference \
            ${params_} \
            --cpus ${task.cpus} \
            --prefix ${base}
        """
}

process vcf_legend {
    tag "legend_${dataset}_${chrm}"
    publishDir "${params.outdir}/${chrm}/legends", mode:'copy', overwrite: true
    label "bigmem"

    input:
        tuple val(dataset), val(chrm), file(dataset_vcf), file(dataset_vcf_idx)

    output:
        tuple val(chrm), val(dataset), file("${base}.legend.gz")

    script:
        base = file(dataset_vcf.baseName).baseName
        """
        echo -e 'id position a0 a1 afr.aaf super_pop.aaf afr.maf super_pop.maf' > ${base}.legend
        bcftools annotate --set-id +'%CHROM\\_%POS\\_%REF\\_%ALT' ${dataset_vcf} -Ob -o ${base}_temp1.bcf
        bcftools +fill-tags ${base}_temp1.bcf -Ob -o ${base}_temp2.bcf -- -t AF,MAF
        bcftools query -f '%ID %POS %REF %ALT %INFO/AF %INFO/AF %INFO/MAF %INFO/MAF\\n' ${base}_temp2.bcf >> ${base}.legend
        bgzip ${base}.legend
        rm -f ${base}_temp*
        """
}

process vcf_to_bcf {
    tag "bcf_${dataset}_${chrm}"
    publishDir "${params.outdir}/${chrm}/bcfs", mode:'copy', overwrite: true
    label "bigmem"

    input:
        tuple val(dataset), val(chrm), file(dataset_vcf), file(dataset_vcf_idx)

    output:
        tuple val(chrm), val(dataset), file("${base}.bcf"), file("${base}.bcf.csi")

    script:
        base = file(dataset_vcf.baseName).baseName
        """
        bcftools view ${dataset_vcf} -Ob -o${base}.bcf
        tabix -f ${base}.bcf
        """
}