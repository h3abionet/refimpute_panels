#!/usr/bin/env python3

import argparse,time

parser = argparse.ArgumentParser()
parser.add_argument(
    "--in_ped", default="${in_ped}", help="")
parser.add_argument("--out_ped", default="${out_ped}", help="comma separated")

args = parser.parse_args()


def reduce_sample(in_ped, out_ped):
    """
    :param :
    :return:
    """
    out = open(out_ped, 'w')
    for it in open(in_ped):
        it = it.split('\t')
        if len(it[0]) > 20:
            it[0] = it[0][:10]+"_"+it[0][-5:]
            it[1] = it[1][:10]+"_"+it[1][-5:]
        it[5] = '1'
        out.writelines('\t'.join(it))
    out.close()
  
if __name__ == '__main__':
    reduce_sample(args.in_ped, args.out_ped)
