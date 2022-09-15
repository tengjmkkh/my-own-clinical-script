
import os
import argparse
import csv
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input_file", nargs='?', action="store" , type=str, help="cnv file")
parser.add_argument("-bed", "--bed_file", nargs='?', action="store" , type=str, help="bed file")
parser.add_argument("-pqcnv", "--pqcnv", nargs='?', action="store" , type=str, help="pqcnv file")
parser.add_argument("-output", "--output_dir", nargs='?', action="store" , type=str, help="output_dir")
parser.add_argument("-s", "--sample", nargs='?', action="store" , type=str, help="sample name")
parser.add_argument("-n", "--num", nargs='?', action="store" , type=str, help="data id")
args = parser.parse_args()

input_file=args.input_file
bed_file=args.bed_file
pqcnv=args.pqcnv
output_dir=args.output_dir
sample=args.sample
num=args.num


gene_info = {}
pq_info = {}
output_result = []
title = ["data_id","gene","result","tumor","cfdna","type","chr","nm","region"]
data_id = num
output_result.append(title)
with open(input_file,'r') as f:
    data = csv.reader(f)
    for i in data:
        if i[1] == 'cnv':
            continue
        gene_info[i[0]] = {'cnv':i[1],'chrs':i[2],'nm':i[3],'region':i[4]}
with open(pqcnv,'r') as f:
    data = f.readlines()
    sex = 'female'
    for i in data:
        [pq,cnv] = i.split('\t')
        if pq[0] not in ['p','q']:
            pq = pq[0] + pq[1:]
            if pq[0] == 'Y':
                sex = 'male'
        else:
            pq = pq[1:] + pq[0]
        pq_info[pq] = cnv[:-1]
print(sex)
with open(bed_file,'r') as f:
    data = f.readlines()
    for i in data:
        bed_row = i.split(' ')
        if bed_row[1] in gene_info:
            gene_info[bed_row[1]]['pq_site'] = bed_row[2]
for i in gene_info:
    for j in pq_info:
        if j in gene_info[i]['pq_site'][:len(j)]:
            amplification,loss = 3,1.3
            if 'X' in j or 'Y' in j:
                if sex == 'male':
                    amplification,loss = 1.5,0.6
            if float(gene_info[i]['cnv']) / float(pq_info[j]) >= 1.5 and float(gene_info[i]['cnv']) >= amplification:
                if sample =='ffpe':
                    tmp = [data_id,i,'扩增',gene_info[i]['cnv'],'','cnv',gene_info[i]['chrs'],gene_info[i]['nm'],gene_info[i]['region']]
                else:
                    tmp = [data_id,i,'扩增','',gene_info[i]['cnv'],'cnv',gene_info[i]['chrs'],gene_info[i]['nm'],gene_info[i]['region']]
                output_result.append(tmp)          
            elif float(gene_info[i]['cnv']) / float(pq_info[j]) <= 0.75 and float(gene_info[i]['cnv']) <= loss:
                if sample =='ffpe':
                    tmp = data_id,i,'缺失',gene_info[i]['cnv'],'','cnv',gene_info[i]['chrs'],gene_info[i]['nm'],gene_info[i]['region']
                else:
                    tmp = [data_id,i,'缺失','',gene_info[i]['cnv'],'cnv',gene_info[i]['chrs'],gene_info[i]['nm'],gene_info[i]['region']]
                output_result.append(tmp)
with open(input_file[:-4] + '_gain.csv','w') as w:
    writer = csv.writer(w)
    writer.writerows(output_result)
#os.system('curl haplab.haplox.net/api/report/csv?type=cnv -F "import_file=@{S1}" '.format(S1=input_file[:-4] + '_gain.csv'))
