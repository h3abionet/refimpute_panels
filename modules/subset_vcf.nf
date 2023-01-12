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

process extract_site_from_vcf {
    tag "extract_site_${target_name}_${site_name}"
    label "bigmem"

    input:
        tuple val(target_name), val(chrm), val(start), val(end), file(target_vcf), val(site_name), file(site_file)
    output:
        tuple val(target_name), val(chrm), val(start), val(end), file(tag_target_vcf), val(site_name), file(site_file)
    script:
        base = "${target_name}_${site_name}"
        tag_target_vcf = "${base}.bcf"
        """
        tabix -f ${target_vcf}
        bcftools view --regions-file ${site_file} ${target_vcf}  -Ob -o ${tag_target_vcf}
        """
}

process subset_ref_panel_for_pca {
    tag "subset_pca_${ref_panel}"
    label "bigmem"

    input:
        tuple val(dataset), file(dataset_vcf), file(dataset_map), val(ref_panel), file(ref_panel_vcf)
    output:
        tuple val(dataset), file(dataset_vcf_subset), file(dataset_map), val(ref_panel), file(ref_panel_vcf_subset)
    script:
        base = file(ref_panel_vcf.baseName).baseName
        ref_panel_vcf_subset = "${base}_subset_target.bcf"
        ref_panel_map = "${base}_subset_target.map"
        dataset_vcf_subset = "${file(dataset_vcf.baseName).baseName}_subset_target.bcf"
        """
        tabix -f ${ref_panel_vcf}
        bcftools view --regions-file ${dataset_map} ${ref_panel_vcf} | bcftools sort -T . -Ob -o ${ref_panel_vcf_subset}
        tabix -f ${ref_panel_vcf_subset}
        bcftools query -f '%CHROM\\t%POS\\n' -i 'TYPE="snp"' ${ref_panel_vcf_subset} > ${ref_panel_map}
        tabix -f ${dataset_vcf}
        bcftools view --regions-file ${ref_panel_map} ${dataset_vcf} | bcftools sort -T . -Ob -o ${dataset_vcf_subset}
        tabix -f ${dataset_vcf_subset}
        """
}


process sort_vcf {
    tag "sort_${target_name}"
    label "bigmem"

    input:
        tuple val(target_name), val(chrm), val(start), val(end), file(target_vcf), val(site_name), file(site_file)
    output:
        tuple val(target_name), val(chrm), val(start), val(end), file(sorted_vcf), val(site_name), file(site_file)
    script:
        base = file(target_vcf.baseName).baseName
        sorted_vcf = "${base}_sorted.bcf"
        """
        bcftools sort ${target_vcf} -T . -Ob -o ${sorted_vcf}
        """
}

process fill_tags_VCF {
    tag "fill_tags_${dataset}_${chrm}"
    label "bigmem"

    input:
        tuple val(dataset), file(vcf), val(chrm)
    output:
        tuple val(dataset), file(out_vcf), val(chrm)
    script:
        base = file(vcf.baseName).baseName
        out_vcf = "${base}_AF.vcf.gz"
        """
        tabix -f ${vcf}
        bcftools +fill-tags ${vcf}  -Oz -o ${out_vcf} --threads ${task.cpus} -- -t AF,MAF
        """
}

process vcf_freq {
    tag "freq_${dataset}_${chrm}"
    label "bigmem"

    input:
        tuple val(dataset), file(vcf), val(chrm)
    output:
        tuple val(dataset), file(frq), val(chrm)
    script:
        base = file(vcf.baseName).baseName
        frq = "${base}.frq"
        """
        tabix -f ${vcf}
        echo -e 'CHROM\tPOS\tSNP\tREF\tALT\tMAF' > ${frq}
        bcftools annotate --set-id '%CHROM\\_%POS' ${vcf} | \
        bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/MAF\\n' >> ${frq}
        """
}

        // bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' ${vcf} | \

process vcf_map {
    tag "map_${target_name}"
    label "bigmem"

    input:
        tuple val(target_name), val(chrm), val(start), val(end), file(target_vcf), val(site_name), file(site_file)
    output:
        tuple val(target_name), val(chrm), val(start), val(end), file(target_vcf), file("${target_vcf}.csi"), val(site_name), file(map)
    script:
        base = file(target_vcf.baseName).baseName
        map = "${base}.map"
        """
        tabix -f ${target_vcf}
        bcftools annotate  --set-id +'%CHROM\\_%POS' ${target_vcf} | \
        bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\n' > ${map}
        """
}

process vcf_map_simple {
    tag "map_${target_name}"
    label "bigmem"

    input:
        tuple val(target_name), val(chrm), val(start), val(end), file(target_vcf), val(site_name), file(site_file)
    output:
        tuple val(target_name), val(chrm), val(start), val(end), file(target_vcf), file("${target_vcf}.csi"), val(site_name), file(map)
    script:
        base = file(target_vcf.baseName).baseName
        map = "${base}.map"
        """
        tabix ${target_vcf}
        bcftools query -f '%CHROM\\t%POS\\n' ${target_vcf} > ${map}
        """
}

process vcf_map_vcf {
    tag "map_${dataset}_${description}"
    label "bigmem"

    input:
        tuple val(dataset), val(description), file(dataset_vcf)
    output:
        tuple val(dataset), val(description), file(dataset_vcf), file(dataset_map)
    script:
        base = file(dataset_vcf.baseName).baseName
        dataset_map = "${base}.map"
        """
        tabix ${dataset_vcf}
        bcftools norm -c x -f ${params.reference_genome} ${dataset_vcf} |\
        bcftools +prune -l 0.5 -w 50kb |\
        bcftools query -f '%CHROM\\t%POS\\n' -i 'TYPE="snp"' > ${dataset_map}
        """
}

process info_freq {
    tag "info_freq_${dataset}_${ref_panel}"
    label "bigmem"

    input:
        tuple val(dataset), file(info), val(ref_panel)
    output:
        tuple val(dataset), val(ref_panel), file(frq), file(info)
    script:
        base = info.baseName
        frq = "${base}.frq"
        """
        awk '{print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6}' ${info} > ${frq}
        """
}

process get_vcf_sites {
    tag "get_sites_${target_name}_${chrm}:${chunk_start}-${chunk_end}_${tagName}"
    label "bigmem"

    input:
        tuple val(chrm), val(chunk_start), val(chunk_end), val(target_name), val(ref_name), file(target_vcf), val(tagName)
    output:
        tuple val(chrm), val(chunk_start), val(chunk_end), val(target_name), val(ref_name), file(sites), val(tagName)
    script:
        base = file(target_vcf.baseName).baseName
        sites = "${base}.sites"
        """
        tabix -f ${target_vcf}
        echo -e 'SNP' > ${sites}
        bcftools annotate  --set-id +'%CHROM\\_%POS' ${target_vcf} | \
        bcftools query  -f '%ID\\n' >> ${sites}
        """
}

process subset_vcf_sites {
    tag "subset_vcf_sites_${target_name}_${site_name}"
    label "bigmem"

    input:
        tuple val(target_name), val(chrm), val(start), val(end), file(target_vcf), val(site_name), file(site_file)
    output:
        tuple val(target_name), val(chrm), val(start), val(end), file(tag_target_vcf)
    script:
        base = "${target_name}_${site_name}"
        tag_target_vcf = "${base}.bcf"
        """
        tabix -f ${target_vcf}
        bcftools view --regions-file ${site_file} ${target_vcf} | \
        bcftools sort -T . -Ob -o ${tag_target_vcf}
        """
}

process get_vcf_sites_1 {
    tag "get_sites_${target_name}_${chrm}:${chunk_start}-${chunk_end}_${tagName}"
    label "bigmem"

    input:
        tuple val(chrm), val(chunk_start), val(chunk_end), val(target_name), val(ref_name), file(target_vcf), file(target_info), val(tagName)
    output:
        tuple val(chrm), val(chunk_start), val(chunk_end), val(target_name), val(ref_name), file(sites), file(target_info), val(tagName)
    script:
        base = file(target_vcf.baseName).baseName
        sites = "${base}.sites"
        """
        tabix -f ${target_vcf}
        echo -e 'SNP' > ${sites}
        bcftools annotate  --set-id +'%CHROM\\_%POS' ${target_vcf} | \
        bcftools query  -f '%ID\\n' >> ${sites}
        """
}

process generate_chunks {
    tag "generate_chunks_${target_name}"
    label "small"

    input:
        tuple val(target_name), file(mapFile), val(chunk_size)
    output:
        tuple val(target_name), file(chunkFile)
    script:
        if(params.chunk){ chunk = params.chunk } else{ chunk='' }
        chromosomes = ''
        chunkFile = "chunks.txt"
        template "generate_chunks.py"
}

process generate_chunks_vcf {
    tag "generate_chunks_${target_name}"
    label "small"

    input:
        tuple val(target_name), file(vcf), file(vcf_csi), file(mapFile), val(chunk_size)
    output:
        tuple val(target_name), file(vcf), file(vcf_csi), file(chunkFile)
    script:
        if(params.chunk){ chunk = params.chunk } else{ chunk='' }
        chromosomes = ''
        chunkFile = "chunks.txt"
        template "generate_chunks.py"
}


process split_target_to_chunk {
    tag "split_${target_name}_${chrm}:${chunk_start}-${chunk_end}_${tagName}"
    label "bigmem"
    maxForks 30

    input:
        tuple val(target_name), val(chrm), val(chunk_start), val(chunk_end), val(tagName), file(tag_target_vcf), file(tag_target_vcf_csi)
    output:
        tuple val(target_name), val(chrm), val(chunk_start), val(chunk_end), val(tagName), file(vcf_chunk_out)
    script:
        base = file(tag_target_vcf.baseName).baseName
        vcf_chunk_out = "${base}.${chrm}_${chunk_start}-${chunk_end}_${tagName}.bcf"
        start = chunk_start - params.buffer_size
        if(chunk_start.toInteger() - params.buffer_size.toInteger() <= 0){ end = 1 }
        end = chunk_end.toInteger() + params.buffer_size.toInteger()
        """
        bcftools view --regions ${chrm}:${start}-${end} ${tag_target_vcf} -Ob -o ${vcf_chunk_out}
        """
}


process extract_pop {
    tag "extract_pop_vcf_${target_name}_${chrm}_${ref_name}"
    // publishDir "${params.outDir}/imputed/vcfs/${ref_name}/${prefix}/${target_name}", overwrite: true, mode:'link', pattern: '*vcf.gz*'
    label "bigmem"

    input:
        tuple val(chrm), val(target_name), val(ref_name), val(tagName), file(vcf), file(vcf_tbi), file(sample_to_extract), val(prefix)
    output:
        tuple val(chrm), val(target_name), val(ref_name), val(tagName), file(vcf_out), file("${vcf_out}.tbi")
    script:
        base = file(vcf.baseName).baseName
        vcf_out = "${base}_${prefix}.vcf.gz"
        """
        bcftools \
            view \
             \
            --force-samples \
            --samples-file ${sample_to_extract} \
            ${vcf} \
            -Oz -o ${vcf_out} 
        tabix ${vcf_out} 
        """
}


process remove_pop {
    tag "remove_pop_vcf_${dataset}_${chrm}"
    label "bigmem5"

    input:
        tuple val(dataset), val(chrm), val(start), val(end), file(vcf), file(sample_to_remove), val(prefix)
    output:
        tuple val(dataset), val(chrm), val(start), val(end), file(vcf_out)
    script:
        base = file(vcf.baseName).baseName
        vcf_out = "${base}_${prefix}.bcf"
        """
        awk '{print \$1}' ${sample_to_remove} > ${sample_to_remove}.samples
        bcftools \
            view \
            --force-samples \
            --samples-file ^${sample_to_remove}.samples \
            --min-ac 1 \
            ${vcf} |\
        bcftools +fill-tags --threads ${task.cpus} -Ob -o ${vcf_out} -- -t AF,MAF
        tabix ${vcf_out} 
        """
}

process keep_samples {
    tag "remove_pop_vcf_${dataset}_${chrm}"
    label "bigmem5"

    input:
        tuple val(dataset), val(chrm), val(start), val(end), file(vcf), file(sample_to_remove), val(prefix)
    output:
        tuple val(dataset), val(chrm), val(start), val(end), file(vcf_out)
    script:
        base = file(vcf.baseName).baseName
        vcf_out = "${base}_${prefix}.bcf"
        """
        awk '{print \$1}' ${sample_to_remove} > ${sample_to_remove}.samples
        bcftools \
            view \
            --force-samples \
            --samples-file ${sample_to_remove}.samples \
            --min-ac 1 \
            ${vcf} |\
        bcftools +fill-tags --threads ${task.cpus} -Ob -o ${vcf_out} -- -t AF,MAF
        tabix ${vcf_out} 
        """
}

process get_sample1 {
    tag "get_sample1_${dataset}_${ref_panel}"
    label "bigmem"

    input:
        tuple val(dataset), file(dataset_vcf), val(ref_panel), file(ref_panel_sample_file)
    output:
        tuple val(dataset), file(dataset_sample), val(ref_panel), file(ref_panel_sample_file)
    script:
        dataset_sample = "${dataset}.samples"
        """
        bcftools query -l ${dataset_vcf} > ${dataset_sample}
        """
}

"""
Step: Add group to evec file from smartpca
"""
process update_evec {
    tag "update_evec_${group_name}"
    label "bigmem"
    input:
        tuple val(group_name), file(group_evec), file(group_eval), file(group_grmjunk), file(group_sample)
    output:
        tuple val(group_name), file(group_evec_update), file(group_eval), file(group_grmjunk), file(group_sample)
    script:
        group_evec_update = "${file(group_evec).baseName}_update.evec"
        evec_file = group_evec
        evec_out = group_evec_update
        annot_file = group_sample
        template "update_evec.py"
}

"""
Step: Plot PCA analysis for group
"""
process plot_pca_group {
    tag "plot_pca_group_${group_name}"
    label "rplot"
    publishDir "${params.outdir}/pca", overwrite: true, mode:'copy', pattern: '*tiff'
    
    input:
        tuple val(group_name), file(group_evec), file(group_eval), file(group_grmjunk), file(group_sample)
    output:
        tuple val(group_name), file(group_evec), file(group_sample), file("${output_pdf}*tiff")
    script:
        output_pdf = "${group_name}"
        input_evec = group_evec
        template "plot_pca.R"
}

process merge_dataset_vcfs {
    tag "merge_dataset_vcfs_${dataset}"
    label "bigmem"
    
    input:
        tuple val(dataset), val(vcfs)
    
    output:
        tuple val(dataset), file(vcf_out_plink)
    
    script:
        vcf_out_plink = "${dataset}_plink.vcf.gz"
        """
        bcftools merge ${vcfs} -m none --force-samples |\
        bcftools view -m2 -M2 -v snps |\
        bcftools sort -T . |\
        bcftools annotate -x INFO,^FORMAT --set-id '%CHROM\\_%POS' -Oz -o ${vcf_out_plink}
        tabix -f ${vcf_out_plink}
        """
}

process vcf_to_plink {
    tag "vcf_to_plink_${dataset}"
    label "bigmem"
    input:
        tuple val(dataset), file(dataset_vcf)
    output:
        tuple val(dataset), file("${base}_pruned.ped"), file("${base}_pruned.map")
    script:
        base = file(dataset_vcf.baseName).baseName
        """
        vcftools --gzvcf ${dataset_vcf} --plink --out ${base}_pruned
        """
}

process reduce_sample_id_size {
    tag "reduce_sample_id_size_${dataset}"
    label "bigmem"
    input:
        tuple val(dataset), file(dataset_ped), file(dataset_map)
    output:
        tuple val(dataset), file(new_dataset_ped), file(dataset_map)
    script:
        new_dataset_ped = "${dataset_ped.baseName}_new.ped"
        in_ped = dataset_ped
        out_ped = new_dataset_ped
        template "reduce_sample_id_size.py"       
}

process combine_sample_files {
    tag "combine_sample_files_${dataset}_${ref_panel}"
    label "bigmem"
    input:
        tuple val(dataset), file(dataset_sample), val(ref_panel), file(ref_panel_sample_file)
    output:
        tuple val(group_name), file(group_sample_file)
    script:
        sample_name1 = dataset
        sample_file1 = dataset_sample
        sample_name2 = ref_panel
        sample_file2 = ref_panel_sample_file
        group_name = "${dataset}-${ref_panel}"
        group_sample_file = "${group_name}.samples"
        sample_file_out = group_sample_file
        template "combine_sample_files.py"       
}



process smartpca_dataset {
    tag "smartpca_dataset_${dataset}"
    label "smartpca"
    input:
        tuple val(dataset), file(dataset_ped), file(dataset_map)
    output:
        tuple val(dataset), file(dataset_evec), file(dataset_eval), file(dataset_grmjunk)
    script:
        dataset_evec = "${dataset}.evec"
        dataset_eval = "${dataset}.eval"
        dataset_grmjunk = "${dataset}.evec_grmjunk"
        """
        ## Create parameter file for smartpca
        echo -e \
        "genotypename:   ${dataset_ped}
        snpname:         ${dataset_map}
        indivname:       ${dataset_ped}
        evecoutname:     ${dataset_evec}
        evaloutname:     ${dataset_eval}
        altnormstyle:    NO
        numoutevec:      5
        numoutlieriter:  5
        familynames:     NO
        grmoutname:      ${dataset_grmjunk}" > ${dataset}.EIGENSTRAT.par
        ## Run smartpca
        smartpca \
                -p ${dataset}.EIGENSTRAT.par \
                > ${dataset}.EIGENSTRAT.log
        """
}