
#coding:utf-8
import csv
import argparse
import subprocess
from unittest import result
import os


def get_merge(mutation_list,mutation_index):
    Done_row=[]
    new_row=[]
    for mutation in  mutation_list:
        if mutation not in Done_row:
            Ref_row=mutation
            new_result,Done_result=every_merge(new_row,Done_row,Ref_row,mutation_index,mutation_list)
            Done_row+=Done_result
            new_row.append(new_result)              
    return new_row

def every_merge(new_row,Done_row,Ref_row,mutation_index,mutation_list):
    new_row=[]
    for mutation in mutation_list:
        if mutation not in Done_row and mutation!=Ref_row:
            if check_merge(mutation,Ref_row,mutation_index)==False:
                continue
            else:
                Done_row+=[mutation,Ref_row]
                VAF=str(min(float(mutation[mutation_index['VAF_percent']]),float(Ref_row[mutation_index['VAF_percent']])))
                start,end=str(Ref_row[mutation_index['start']]),str(mutation[mutation_index['end']])
                if int(mutation[mutation_index['start']])-int(Ref_row[mutation_index['end']])==1:
                    ref,alt=str(Ref_row[mutation_index['ref']])+str(mutation[mutation_index['ref']]),str(Ref_row[mutation_index['alt']])+str(mutation[mutation_index['alt']])
                else:
                    gap_position=mutation[mutation_index['chr']]+':'+str(Ref_row[mutation_index['end']])+'-'+str(mutation[mutation_index['start']])
                    (status1, output) = subprocess.getstatusoutput('samtools faidx /thinker/net/ctDNA/hg19/hg19.fa {S1}'.format(S1=gap_position))
                    sequence=output.split('\n')[-1]
                    print(sequence)
                    if str(mutation[mutation_index['ref']])!='-' and str(Ref_row[mutation_index['ref']])!='-':
                        ref,alt=str(Ref_row[mutation_index['ref']])+sequence[1:][:-1]+str(mutation[mutation_index['ref']]),str(Ref_row[mutation_index['alt']])+sequence[1:][:-1]+str(mutation[mutation_index['alt']])
                    if str(mutation[mutation_index['ref']])=='-' and str(Ref_row[mutation_index['ref']])=='-':
                        ref,alt=sequence,sequence[0]+str(Ref_row[mutation_index['alt']])+sequence[1:]+str(mutation[mutation_index['alt']])
                    if str(mutation[mutation_index['ref']])=='-' and str(Ref_row[mutation_index['ref']])!='-':
                        ref,alt=str(Ref_row[mutation_index['ref']])+sequence[1:],str(Ref_row[mutation_index['alt']])+sequence[1:]+str(mutation[mutation_index['alt']])
                    if str(mutation[mutation_index['ref']])!='-' and str(Ref_row[mutation_index['ref']])=='-':
                        ref,alt=sequence[:-1]+str(mutation[mutation_index['ref']]),sequence[0]+str(Ref_row[mutation_index['alt']])+sequence[1:][:-1]+str(mutation[mutation_index['alt']])
                ref=ref.strip(' -')
                alt=alt.strip(' -')
                if len(ref)>len(alt):
                    if abs(len(ref)-len(alt)) % 3 == 0:
                        ExonicFunc='nonframeshift deletion'
                    else:
                        ExonicFunc='frameshift deletion'
                else:
                    if abs(len(ref)-len(alt)) % 3 == 0:
                        ExonicFunc='nonframeshift insertion'
                    else:
                        ExonicFunc='frameshift insertion'
                if ref=='':
                    ref='-'
                if alt=='':
                    alt='-'
                new_row=[mutation[0],'indel']+mutation[2:6]+['','',VAF,'-','-','-','-',ExonicFunc,start,end,ref,alt]
                Ref_row=new_row
    if new_row !=[]:
        return new_row,Done_row
    else:
        return Ref_row,Ref_row                  


def check_merge(mutation,Ref_row,mutation_index):
    if  mutation[mutation_index['chr']]!= Ref_row[mutation_index['chr']]:
        return False
    typelist=[mutation[mutation_index['type']],Ref_row[mutation_index['type']]]
    if 'indel' in typelist and int(mutation[mutation_index['start']])-int(Ref_row[mutation_index['end']]) > 6:
        return False
    elif 'indel' not in typelist and int(mutation[mutation_index['start']])-int(Ref_row[mutation_index['end']]) > 3:
        return False
    if int(mutation[mutation_index['start']])-int(Ref_row[mutation_index['end']])<=0:
        return False
    if abs(float(mutation[mutation_index['VAF_percent']])-float(Ref_row[mutation_index['VAF_percent']]))>max(float(mutation[mutation_index['VAF_percent']]),float(Ref_row[mutation_index['VAF_percent']]))*0.1:
        return False
    return True

    #=================解析MrBam  title===================    
def analyse_title(row,target):
    analyse_param = {}
    for target_title in target:
        for title in range(len(row.split(','))):
            if row.split(',')[title] == target_title:
                analyse_param[target_title] = title
    return analyse_param
def get_mutation_index(file):
    with open(file,'r',encoding='gbk') as f:
        data = f.readlines()
        target = ['Data_id','type','gene','chr','NM','exon','base','AA','VAF_percent','cosmic','cosmic_id','cosmic_occu','dbsnp','ExonicFunc','start','end','ref','alt']
        title_analyse_bool = True
        for line in data:
            if title_analyse_bool == True:
                analyse_param = analyse_title(line[:-1],target)
            title_analyse_bool = False
    return analyse_param

def main():
    #获取参数
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputfile", nargs='?', action="store" , type=str, help="file name")
    parser.add_argument("-o", "--outputdir", nargs='?', action="store" , type=str, help="outputdir")
    args = parser.parse_args()
    file = args.inputfile
    outputdir = args.outputdir
    mutation_list=[]
    mutation_index=get_mutation_index(file)
    title=[['Data_id','type','gene','chr','NM','exon','base','AA','VAF_percent','cosmic','cosmic_id','cosmic_occu','dbsnp','ExonicFunc','start','end','ref','alt']]
    with open(file,'r',encoding='utf-8')as f:
        next(f)
        for line in f.readlines():
            line=line[:-1]
            col=line.split(',')
            mutation_list.append(col)
    if  mutation_list!=[]:
        mutation_list.sort(key=lambda mutation_list:int(mutation_list[mutation_index['start']]),reverse=False) #start升序排列
        mutation_list=get_merge(mutation_list,mutation_index)
    Result=title+mutation_list
    print(Result)
    with open(outputdir +'/mutation_result.csv','w',encoding='gbk') as w:
        writer = csv.writer(w)
        writer.writerows(Result)

                
if __name__ == "__main__":
    main()



