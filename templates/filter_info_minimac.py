#!/usr/bin/env python2.7

import argparse, sys

parser = argparse.ArgumentParser()
parser.add_argument("--infoFiles", default="${infos}", help="")
parser.add_argument("--outWell_imputed", default="${well_out}", help="")
parser.add_argument("--outSNP_acc", default="${acc_out}", help="")
parser.add_argument("--infoCutoff", default="${impute_info_cutoff}", help="")

args = parser.parse_args()


def filter_info(infoFiles, infoCutoff, outWell_imputed, outSNP_acc):
    """
    Return:
        well_imputed: certainty >= 1
        SNP_concordance: concord_type0 != -1
    """
    well_imputed = {}
    SNP_concordance = {}
    count = 0
    infoFiles = infoFiles.split(',')
    header = []
    outWell_imputed_out = open(outWell_imputed + ".tsv", 'w')
    outWell_imputed_snp_out = open(outWell_imputed + "_snp.tsv", 'w')
    outWell_bad_out = open(outWell_imputed + "_bad_info.tsv", 'w')
    outWell_qc_out = open(outWell_imputed + "_bad_qc.tsv", 'w')
    outSNP_accuracy_out = open(outSNP_acc + ".tsv", 'w')
    for infoFile in infoFiles:
        infoFile = infoFile.strip().split('==')
        dataset = infoFile[0]
        info = infoFile[1]
        chip = infoFile[2]
        sites = [ snp.strip() for snp in open(infoFile[3]).readlines() ]
        well_imputed[dataset] = []
        SNP_concordance[dataset] = []
        # print info
        for line in open(info):
            data = line.strip().split()
            if "SNP" in line and "Rsq" in line:
                if len(header) == 0:
                    header = data
                    snp_idx = header.index("SNP")
                    info_idx = header.index("Rsq")
                    conc_idx = header.index("EmpRsq")
                    outWell_imputed_out.writelines('\\t'.join([dataset] + data) + '\\n')
                    outWell_imputed_snp_out.writelines(data[1] + '\\n')
                    outSNP_accuracy_out.writelines('\\t'.join([dataset] + data) + '\\n')
            else:
                # print info_idx, data
                r2 = data[info_idx]
                snp = data[snp_idx]
                conc = data[conc_idx]
                try:
                    if float(r2) >= float(infoCutoff):
                        if snp in sites:
                            outWell_imputed_out.writelines('\\t'.join([dataset] + data) + '\\n')
                            outWell_imputed_snp_out.writelines(
                                '\\t'.join([dataset] + data) + '\\n')
                        else:
                            outWell_qc_out.writelines(
                                '\\t'.join([dataset] + data) + '\\n')
                    if float(r2) <= 0.3:
                        # print r2, infoCutoff
                        outWell_bad_out.writelines(snp + '\\n')
                except:
                    print r2, infoCutoff
                    sys.exit(1)
                if conc != '-':
                    outSNP_accuracy_out.writelines('\\t'.join([dataset] + data) + '\\n')
                count += 1
    outWell_imputed_out.close()
    outWell_imputed_snp_out.close()
    outSNP_accuracy_out.close()
    outWell_bad_out.close()
    outWell_qc_out.close()

if args.infoFiles and args.infoCutoff:
    filter_info(args.infoFiles, args.infoCutoff, args.outWell_imputed, args.outSNP_acc)

