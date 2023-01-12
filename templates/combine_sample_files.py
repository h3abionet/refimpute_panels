#!/usr/bin/env python3

import argparse,time

parser = argparse.ArgumentParser()
parser.add_argument(
    "--sample_file1", default="${sample_file1}", help="Sample file")
parser.add_argument(
    "--sample_file2", default="${sample_file2}", help="Sample file")
parser.add_argument(
    "--sample_name1", default="${sample_name1}", help="Sample name")
parser.add_argument(
    "--sample_name2", default="${sample_name2}", help="Sample name")
parser.add_argument("--sample_file_out",
                    default="${sample_file_out}", help="Sample file")

args = parser.parse_args()


def combine_sample_files(sample_name1, sample_file1, sample_name2, sample_file2, sample_file_out):
    """
    :param :
    :return:
    """
    out = open(sample_file_out, 'w')
    sample_data1 = [it.strip().split('\\t')
                    for it in open(sample_file1).readlines()]
    sample_data2 = [it.strip().split('\\t')
                    for it in open(sample_file2).readlines()]

    datas = {}

    for data in sample_data1:
        sample_id = data[0]
        if len(sample_id) > 20:
            sample_id = data[0][:10]+"_"+data[0][-5:]
        pop = sample_name1
        project = sample_name1
        dataset = sample_name1
        datas[sample_id] = [sample_id, pop, project, dataset]

    for data in sample_data2:
        sample_id = data[0]
        if len(sample_id) > 20:
            sample_id = data[0][:10]+"_"+data[0][-5:]
        pop = data[1]
        project = data[3]
        dataset = sample_name2
        datas[sample_id] = [sample_id, pop, project, dataset]

    for data in datas:
        out.writelines('\\t'.join(datas[data])+'\\n')
    
    out.close()
  
if __name__ == '__main__':
    combine_sample_files(args.sample_name1, args.sample_file1,
                  args.sample_name2, args.sample_file2, args.sample_file_out)
