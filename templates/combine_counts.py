#!/usr/bin/env python2.7

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--csv_counts", default="${csv_counts}", help="")
parser.add_argument("--label", default="${dataset}", help="")
parser.add_argument("--csv_chrms", default="${csv_chrms}", help="")
parser.add_argument("--csv_out", default="${csv_out}", help="")

args = parser.parse_args()

def combine_csv_counts(csv_counts, label, csv_chrms, csv_out):
    """
    :param pop_freqs:
    :return:
    """
    files = csv_counts.split(',')
    # labels = csv_labels.split(',')
    chrms = csv_chrms.split(',')
    classes = ['Group', 'CHROM']
    datas = {}
    for i in range(len(sorted(chrms))):
        data = []
        # pop = labels[i]
        chrm = chrms[i]
        datas[i] = [label, chrm]
        for line in open(files[i]):
            line = line.split(';')
            classe = line[0].strip()
            count = line[1]
            datas[i].append(count.strip())
            if classe not in classes:
                classes.append(classe)

    out = open(csv_out, 'w')
    out.writelines('\\t'.join(classes)+'\\n')
    for i in datas:
        out.writelines('\\t'.join(datas[i])+'\\n')
    out.close()

if __name__ == '__main__':
    combine_csv_counts(args.csv_counts, args.label,
                       args.csv_chrms, args.csv_out)

