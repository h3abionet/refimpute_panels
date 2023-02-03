#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { split_vcf_chromosome; split_vcf_chunk; filter_min_ac; filter_f_missing; split_multi_allelic; qc_dupl; get_chromosome; get_chromosome_vcf; check_chromosome; check_files; check_chromosome_vcf; check_mismatch; no_mismatch ; target_qc as target_qc; target_qc as target_qc1; qc_site_missingness as qc_site_missingness1; qc_site_missingness as qc_site_missingness2; sites_only ; combine_vcfs ; combine_infos; combine_csvs as combine_freqs; combine_vcfs_chrm; filter_simple_snps_only; combine_vcfs_sites_only; combine_vcf_sites as combine_vcf_sites_1; combine_vcf_sites as combine_vcf_sites_2; combine_vcf_sites as combine_vcf_sites_3; combine_vcf_sites as combine_vcf_sites_4; combine_vcf_sites as combine_vcf_sites_5;} from './modules/qc'
include { generate_chunks_vcf; split_target_to_chunk; vcf_map; vcf_map_simple } from './modules/subset_vcf'
include { phasing_vcf_no_ref_chunk } from './modules/phasing'
include { vcf_to_m3vcf; vcf_to_bcf; vcf_legend } from './modules/utils'

// def check_files(file_list) {
//     file_list.each { myfile ->
//         if (!file(myfile).exists() && !file(myfile).isFile()) exit 1, "|-- ERROR: File ${myfile} not found. Please check your config file."
//     }
// }

process phasing_vcf {
    tag "phase_${dataset}_${chrm}"
    publishDir "${params.outdir}/${dataset}/vcfs_phased", overwrite: true, mode:'copy'
    label "bigmem"

    input:
        tuple chrm, dataset, file(dataset_vcf), file(dataset_sample), refpanel, file(refpanel_vcf), file(eagle_genetic_map)

    output:
        tuple chrm, dataset, file("${file_out}.vcf.gz")

    script:
        base = file(dataset_vcf.baseName).baseName
        file_out = "${base}_phased"
        """
        tabix ${dataset_vcf}
        tabix ${refpanel_vcf}
        eagle \
            --numThreads=${task.cpus} \
            --vcfTarget=${dataset_vcf} \
            --geneticMapFile=${eagle_genetic_map} \
            --vcfRef=${refpanel_vcf} \
            --vcfOutFormat=z \
            --chrom=${chrm} \
            --bpStart=${chunk_start} \
            --bpEnd=${chunk_end} \
            --outPrefix=${file_out} 2>&1 | tee ${file_out}.log
        """
}


process phasing_vcf_no_ref {
    tag "phase_${dataset}_${chrm}"
    publishDir "${params.outdir}/${dataset}/vcfs_phased", overwrite: true, mode:'copy'
    label "bigmem10"

    input:
        tuple dataset, chrm, file(dataset_vcf), file(eagle_genetic_map)

    output:
        tuple chrm, dataset, file("${file_out}.vcf.gz")

    script:
        base = file(dataset_vcf.baseName).baseName
        file_out = "${base}_phased"
        """
        eagle \
            --numThreads=${task.cpus} \
            --vcf=${dataset_vcf} \
            --geneticMapFile=${eagle_genetic_map} \
            --vcfOutFormat=z \
            --chrom=${chrm} \
            --outPrefix=${file_out} 2>&1 | tee ${file_out}.log
        tabix "${file_out}.vcf.gz"
        """
}

process tabix_phasing_vcf_no_ref {
    tag "tabix_${dataset}_${chrm}"
    publishDir "${params.outdir}/${dataset}/vcfs_phased_no_ref", overwrite: true, mode:'copy'
    label "bigmem"

    input:
        tuple chrm, dataset, file(dataset_vcf)

    output:
        tuple chrm, dataset, file(dataset_vcf), file("${dataset_vcf}.tbi")

    script:
        """
        tabix ${dataset_vcf}
        """
}

process convertToImp5 {
    tag "convertToImp5_${dataset}_${chrm}"
    publishDir "${params.outdir}/${dataset}/imp5", overwrite: true, mode:'copy'
    label "impute5"

    input:
        tuple val(dataset), val(chrm), file(dataset_vcf), file(dataset_vcf_idx)

    output:
        tuple val(dataset), val(chrm), file(imp5), file("${imp5}.idx")

    script:
        imp5 = "${dataset_vcf.getSimpleName()}.imp5"
        // TODO rename imp5Converter_1.1.5_static to imp5Converter in container
        """
        imp5Converter_1.1.5_static --h ${dataset_vcf} --thread ${task.cpus} --r ${chrm} --o ${imp5}
        """
}

workflow preprocess{
    take:

    main:
        check_files([params.eagle_genetic_map])

        datasets = []
        params.datasets.each { dataset, dataset_vcf, dataset_sample ->
            datas = []
            dataset_vcfs = file(dataset_vcf)
            if (dataset_vcfs instanceof List){
                dataset_vcfs.each{ vcf ->
                    datas = [ dataset ]
                    if( file(vcf).getExtension() == "gz" ){
                        vcf_index = "${vcf}.tbi"
                    }
                    else if( file(vcf).getExtension() == "bcf" ){
                        vcf_index = "${vcf}.csi"
                    }
                    else{
                        vcf_index = "${vcf}.tbi"
                    }
                    check_files([ vcf, vcf_index ])
                    datas << file(vcf)
                    datas << file(vcf_index)
                    datasets << datas
                }
            }
            else{
                datas = [ dataset ]
                if( file(dataset_vcf).getExtension() == "gz" ){
                    dataset_vcf_index = "${dataset_vcf}.tbi"
                }
                else if( file(dataset_vcf).getExtension() == "bcf" ){
                    dataset_vcf_index = "${dataset_vcf}.csi"
                }
                else{
                    dataset_vcf_index = "${dataset_vcf}.tbi"
                }
                check_files([ dataset_vcf, dataset_vcf_index ])
                datas << file(dataset_vcf)
                datas << file(dataset_vcf_index)
                datasets << datas
            }
        }
        
        datasets_ch = Channel.from(datasets)

        filter_simple_snps_only(datasets_ch)

        get_chromosome(filter_simple_snps_only.out)

        generate_chunks_vcf(get_chromosome.out.map{ dataset, vcf, vcf_idx, map_file -> [ dataset, file(vcf), file(vcf_idx), file(map_file), params.chunk_size ] })
        chunks_datas = generate_chunks_vcf.out.flatMap{ dataset, vcf, vcf_idx, chunk_file ->
            datas = []
            chunks = file(chunk_file).text.split()
            chunks.each{ chunk_data ->
                data = chunk_data.split(',')
                chrm = data[0]
                chunk_start = data[1]
                chunk_end = data[2]
                datas << [dataset, chrm, chunk_start, chunk_end, dataset, file(vcf), file(vcf_idx)]
            }
            return datas
        }
        //TODO: put this in a workflow so that combine_vcf_sites is not called multiple times
        split_target_to_chunk(chunks_datas)
        combine_vcf_sites_1(split_target_to_chunk.out.map{ dataset, chrm, start, end, tagname, vcf -> [ dataset, "${chrm}_split_chunks", file(vcf) ] })

        // Checkk REF mismacthes 
        check_mismatch(split_target_to_chunk.out.map{ dataset, chrm, start, end, tagname, vcf -> [ dataset, chrm, start, end, file(vcf), file(params.reference_genome) ] })
        check_mismatch.out.map{ dataset, vcf, chrm, start, end, warn, summary -> no_mismatch(dataset, warn, summary) }

        // QC
        filter_f_missing(split_target_to_chunk.out.map{ dataset, chrm, start, end, tagname, vcf -> [ dataset, chrm, start, end, file(vcf), '' ] })
        combine_vcf_sites_2(filter_f_missing.out.map{ dataset, chrm, start, end, vcf -> [ dataset, "${chrm}_filter_missingness", file(vcf) ] })

        // qc_dupl(split_target_to_chunk.out.map{ dataset, chrm, start, end, tagname, vcf -> [ dataset, chrm, start, end, file(vcf)] })

        split_multi_allelic(filter_f_missing.out)
        combine_vcf_sites_3(split_multi_allelic.out.map{ dataset, chrm, start, end, vcf -> [ dataset, "${chrm}_split_multi_allelic", file(vcf) ] })

        filter_min_ac(split_multi_allelic.out.map{ dataset, chrm, start, end, vcf -> [ dataset, chrm, start, end, file(vcf), " --min-ac ${params.min_ac} --max-alleles ${params.max_alleles} --min-alleles ${params.min_alleles} -v snps "  ] })
        combine_vcf_sites_4(filter_min_ac.out.map{ dataset, chrm, start, end, vcf -> [ dataset, "${chrm}_filter_min_ac", file(vcf) ] })
        
        // filter_min_ac(split_multi_allelic.out.map{ dataset, chrm, start, end, vcf -> [ dataset, chrm, start, end, file(vcf), " --min-ac ${params.min_ac} "  ] })
        vcf_map_simple(filter_min_ac.out.map{ dataset, chrm, start, end, vcf -> [ dataset, chrm, start, end, file(vcf), '', '' ] })
        // vcf_map_simple(split_multi_allelic.out.map{ dataset, chrm, start, end, vcf -> [ dataset, chrm, start, end, file(vcf), '', '' ] })

    emit:
        qc_data = vcf_map_simple.out
}

workflow phasing {
    take: data

    main:
        // // Phasing without a reference
        phasing_ch = data.flatMap{ dataset, chrm, start, end, vcf_file, vcf_idx, chip_name, map_file -> 
            phasing_ch_data = [ dataset, chrm, start, end, file(vcf_file), file(params.eagle_genetic_map) ]
            if(!file(map_file).isEmpty()){
                return [ phasing_ch_data ]
            }
        }
        phasing_vcf_no_ref_chunk(phasing_ch)
        combine_vcfs(phasing_vcf_no_ref_chunk.out.groupTuple(by:[0,1]).map{ dataset, chrm, starts, ends, vcfs -> [ dataset, chrm, vcfs ] })
        combine_vcf_sites_5(phasing_vcf_no_ref_chunk.out.map{ dataset, chrm, start, end, vcf -> [ dataset, "${chrm}_phasing_vcf_no_ref_chunk", file(vcf) ] })

    emit:
        phased_data = combine_vcfs.out
}

workflow impute5convert{
    take: data

    main:
        convertToImp5(data).view()
    emit:
        data
}

workflow{
    
    preprocess()
    phasing(preprocess.out.qc_data)

    ///// For Minimac4 
    vcf_to_m3vcf(phasing.out.phased_data)
    vcf_to_bcf(phasing.out.phased_data)
    vcf_legend(phasing.out.phased_data)

    ///// For IMPUTE5
    // impute5convert(phasing.out.phased_data)
}