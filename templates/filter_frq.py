#!/usr/bin/env python3
'''

'''
import argparse,sys,time

parser = argparse.ArgumentParser()
parser.add_argument("--dataset_frq",
                    default="${dataset_frq}", help="dataset freq 1")
parser.add_argument("--dataset", default="${dataset}", help="dataset name 1")
parser.add_argument("--ref_panel",
                    default="${ref_panel}", help="dataset name 2")
parser.add_argument(
    "--ref_panel_frq", default="${ref_panel_frq}", help="dataset freq 2")
parser.add_argument("--frq_output", default="${frq_output}", help="")

args = parser.parse_args()


def filter_frq(dataset, dataset_frq, ref_panel, ref_panel_frq, frq_output):
    """
    :param bedFile:
    :param frq_file:
    :param frq_output:
    :return:
    """
    dataset_freqs = {}
    out = open(frq_output, 'w')
    nline = 1
    for line in open(dataset_frq):
        data = line.strip().split('\t')
        if nline == 1:
            snp_idx = data.index("SNP")
            maf_idx = data.index("MAF")
        snp = data[snp_idx]
        maf = data[maf_idx]
        if snp not in dataset_freqs:
            # if ',' not in maf:  # and len(maf.split(',')) == 1: # Only take bi-allelic
            dataset_freqs[snp] = maf
            #     try:
            #         if float(maf) <= 0.5:
            #             dataset_freqs[snp] = maf
            #         else:
            #             dataset_freqs[snp] = str(float(maf)-0.5)
            #     except:
            #         print(line)
        nline += 1
    nline = 1
    ref_panel_freqs = {}
    for line in open(ref_panel_frq):
        data = line.replace("#", '').strip().split('\t')
        if nline == 1:
            snp_idx = data.index("SNP")
            maf_idx = data.index("MAF")
        snp = data[snp_idx]
        maf = data[maf_idx]
        # print(snp)
        # time.sleep(4)
        if nline == 1:
            out.writelines(snp+"\\t"+dataset+"_AF"+"\\t"+ref_panel+"_AF\\n")
        else:
            try:
                if snp in dataset_freqs:
                    if snp not in ref_panel_freqs:  # To avoid duplicates
                        # and len(maf.split(',')) == 1: # Only take bi-allelic
                        # if ',' not in maf:
                        ref_panel_freqs[snp] = maf
                        #     try:
                        #         if float(maf) <= 0.5:
                        #             ref_panel_freqs[snp] = maf
                        #         else:
                        #             ref_panel_freqs[snp] = str(float(maf)-0.5)
                        #     except:
                        #         print(line)
            except:
                print(line)
        nline += 1
    for snp in ref_panel_freqs:
        out.writelines(
            snp+"\\t"+dataset_freqs[snp]+"\\t"+ref_panel_freqs[snp]+"\\n")
    out.close()


filter_frq(args.dataset, args.dataset_frq, args.ref_panel,
           args.ref_panel_frq, args.frq_output)
