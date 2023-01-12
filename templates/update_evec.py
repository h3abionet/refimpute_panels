#!/usr/bin/env python3

import argparse,sys
import time

parser = argparse.ArgumentParser()
parser.add_argument("--evec_file", default="${evec_file}", help="")
parser.add_argument("--evec_out", default="${evec_out}", help="")
parser.add_argument("--annot_file", default="${annot_file}", help="")
parser.add_argument("--shapes", default="", help="")
args = parser.parse_args()

def update_group_smartpca_evec(evec_file, evec_out, annot_file, shapes=''):
    """
    Receive an ind EIGENSTRAT file and update phenotype to group using annotation file provided
    """
    groups = {}
    # pop_shapes = {}
    # for pops in shapes.split(","):
    #     pops = pops.split(":")
    #     pop_shapes[pops[0].strip()] = pops[1].strip()
    for line in open(annot_file):
        line = line.strip().split()
        groups[line[0]] = [ line[1], line[2], line[3] ]#, pop_shapes[line[1]]]
    out = open(evec_out, 'w')
    for line in open(evec_file):
        try:
            line = line.strip().split()
            ind_id = line[0]
            line[-1] = '\\t'.join(groups[ind_id])
            out.writelines('\\t'+'\\t'.join(line)+'\\n')
        except:
            pass
    out.close()

if __name__ == "__main__":
    update_group_smartpca_evec(args.evec_file, args.evec_out, args.annot_file, args.shapes)
