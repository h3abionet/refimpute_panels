#!/usr/bin/env python3

import argparse,time

parser = argparse.ArgumentParser()
parser.add_argument(
    "--meta_file", default="/cbio/projects/001/scratch/gerrit/all.v6.meta", help="")
parser.add_argument("--pop_include", default="${pop_include}", help="comma separated")

args = parser.parse_args()

def get_sample(meta_file, pop_include):
    """
    :param :
    :return:
    """
    not_to_include = ['GBR', 'FIN', 'CHS', 'PUR', 'CDX', 'CLM', 'IBS', 'PEL', 'PJL', 'KHV', 'BEB', 'STU', 'ITU', 'CEU', 'CHB', 'JPT', 'MXL', 'TSI', 'GIH', 'Han', 'NorthernHan', 'Dai', 'Lahu', 'Naxi', 'She', 'Yi', 'Miao', 'Tujia', 'Tu', 'Xibo', 'Hezhen', 'Mongolian', 'Daur', 'Oroqen', 'Cambodian', 'Japanese', 'Yakut', 'Uygur', 'Hazara','Burusho', 'Pathan', 'Sindhi', 'Kalash', 'Brahui', 'Balochi', 'Makrani', 'Adygei', 'Russian', 'French', 'Basque', 'Orcadian', 'BergamoItalian', 'Tuscan', 'Sardinian', 'Druze', 'Palestinian', 'Bedouin', 'Mozabite', 'Bougainville', 'PapuanSepik', 'PapuanHighlands', 'Pima', 'Maya', 'Colombian', 'Surui', 'Karitiana', 'unknown']   
    kg = ['GBR', 'FIN', 'CHS', 'PUR', 'CDX', 'CLM', 'IBS', 'PEL', 'PJL', 'KHV', 'BEB', 'STU', 'ITU', 'CEU', 'CHB', 'JPT', 'MXL', 'TSI', 'GIH', 'ACB', 'GWD', 'ESN', 'MSL', 'YRI', 'LWK', 'ASW']   
    pops = {}
    pop_size = {}

    path_ = "/cbio/dbs/refpanels/H3AR6x/meta/all.v6.meta"
    out = open(f'{path_}.african-only.tsv', 'w')
    out_no_aa = open(f'{path_}.non-african-only.tsv', 'w')
    out_ext = open(f'{path_}.extended.tsv', 'w')
    out_ext_pops = open(f'{path_}.extended.pops.tsv', 'w')
    out_kg = open(f'{path_}.kg.tsv', 'w')
    out_not_kg = open(f'{path_}.not_kg.tsv', 'w')
    
    pop_include = pop_include.split(',')
    for line in open(meta_file):
        line1 = line.strip().split('\t')
        pop = line.strip().split('\t')[1]
        project = line.strip().split('\t')[3]
        if pop not in not_to_include:
            out.writelines(line)
            out_ext.writelines('\t'.join(line1+["African-Ancestry"])+'\n')
            if pop not in pops:
                pops[pop] = [pop, project, "African-Ancestry"]
                pop_size[pop] = 0
            pop_size[pop] += 1
        else:
            out_no_aa.writelines(line)
            out_ext.writelines('\t'.join(line1+["Non-African-Ancestry"])+'\n')
            if pop not in pops:
                pops[pop] = [pop, project, "Non-African-Ancestry"]
                pop_size[pop] = 0
            pop_size[pop] += 1
        
        # For 1KG
        if pop in kg:
            out_kg.writelines(line)
        else:
            out_not_kg.writelines(line)

    for pop in pops:
        out_ext_pops.writelines('\t'.join(pops[pop]+[str(pop_size[pop])])+'\n')

    out.close()
    out_no_aa.close()
    out_ext.close()
    out_ext_pops.close()
    out_kg.close()
    out_not_kg.close()
    # time.sleep(4)
  
if __name__ == '__main__':
    get_sample(args.meta_file, args.pop_include)
