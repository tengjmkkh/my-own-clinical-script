#Get Information From MrBam File to nobias-GB1803(-baseline) File
#Author Tengjm
#version:1.3-2022.09.13 ===filt rows
import os
import subprocess
import csv
import re
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--snv_MrBam_file", action="store" , type=str, help="snv_MrBam.txt")
    parser.add_argument("-i", "--indel_MrBam_file", action="store" , type=str, help="indel_MrBam.txt")   
    parser.add_argument("-o", "--outputdir", action="store" , type=str, help="output dir")
    parser.add_argument("-t", "--sampletype", action="store" , type=str, help="pair/single")
    args = parser.parse_args()
    #=================获取参数===================
    snv_MrBam_file = args.snv_MrBam_file
    indel_MrBam_file = args.indel_MrBam_file
    outputdir=args.outputdir
    sampletype=args.sampletype
    sample=snv_MrBam_file.split('/')[-1].split('.')[0]
    dataid=sample.split('_')[-1]
    if sampletype=='pair':
        title=[["Data_id","gene","chr","NM","exon","base","AA","VAF_percent","cosmic","cosmicid","chromosome","position","genotype","total","sr","sv","sr_for","sr_rev","sv_for","sv_rev","Multiple_Overlap_Alt","Multiple_Nonoverlap_Alt","Multiple_Single_Alt","One_Overlap_Alt","One_Nonoverlap_Alt","One_Single_Alt","unique_alt","unique_alt_overlap","Gunique_alt","Gunique_alt_overlap","cosDetails","cosCount","chrs","start","end","ref","alt","FuncRefGene","GeneDetail","ExonicFunc","g1000aug2015_all","snp138","sift","gtotal","gsr","gsv","gpercent","gMultiple_Overlap_Alt","gMultiple_Nonoverlap_Alt","gMultiple_Single_Alt","gOne_Overlap_Alt","gOne_Nonoverlap_Alt","gOne_Single_Alt"]]
    elif sampletype=='single':
        title=[["Data_id","gene","chr","NM","exon","base","AA","VAF_percent","cosmic","cosmicid","chromosome","position","genotype","total","sr","sv","sr_for","sr_rev","sv_for","sv_rev","Multiple_Overlap_Alt","Multiple_Nonoverlap_Alt","Multiple_Single_Alt","One_Overlap_Alt","One_Nonoverlap_Alt","One_Single_Alt","unique_alt","unique_alt_overlap","cosDetails","cosCount","chrs","start","end","ref","alt","FuncRefGene","GeneDetail","ExonicFunc","g1000aug2015_all","snp138","sift"]]
    for MrBamfile in [snv_MrBam_file,indel_MrBam_file]:
        if 'snv' in MrBamfile:
            type='snv'
        elif 'indel' in MrBamfile:
            type='indel'
        MrBam_info=get_MrBam_info(MrBamfile,dataid,sampletype)
        #===================按丰度排序====================================
        MrBam_info.sort(key=lambda MrBam_info:float(MrBam_info[7]),reverse=True)
        #=================================================================
        result=title+MrBam_info
        get_result(outputdir,result,sample,type)
        
    #=================解析MrBam  title===================    
def analyse_title(row,target):
    analyse_param = {}
    for target_title in target:
        for title in range(len(re.split('[\t  ]',row))):
            if re.split('[\t  ]',row)[title] == target_title:
                analyse_param[target_title] = title
    return analyse_param
def get_MrBam_index(MrBamfile):
    with open(MrBamfile,'r',encoding='utf-8') as f:
        data = f.readlines()
        target = ['Chr','Start','End','Ref','Alt','Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene','cytoBand','genomicSuperDups','esp6500siv2_all','1000g2015aug_all','1000g2015aug_afr','1000g2015aug_eas','1000g2015aug_eur','snp138','SIFT_score','SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_pred','Polyphen2_HVAR_score','Polyphen2_HVAR_pred','LRT_score','LRT_pred','MutationTaster_score','MutationTaster_pred','MutationAssessor_score','MutationAssessor_pred','FATHMM_score','FATHMM_pred','RadialSVM_score','RadialSVM_pred','LR_score','LR_pred','VEST3_score','CADD_raw','CADD_phred','GERP++_RS','phyloP46way_placental','phyloP100way_vertebrate','SiPhy_29way_logOdds','cosmic87','CLNALLELEID','CLNDN','CLNDISDB',' CLNREVSTAT','CLNSIG','Otherinfo']
        title_analyse_bool = True
        for line in data:
            if title_analyse_bool == True:
                analyse_param = analyse_title(line,target)
                title_analyse_bool = False
                return analyse_param
    #======================获取MrBam 信息==========================
def get_MrBam_info(MrBamfile,dataid,sampletype):
    MrBam_index=get_MrBam_index(MrBamfile)
    MrBam_info=[]
    with open(MrBamfile,'r',encoding='utf-8') as f:
        next(f)
        data = f.readlines()
        for line in data:
            col=line.split('\t')
            cosCount=0
            chr,chromosome,position,chrs,start,end,ref,alt=col[MrBam_index['Chr']],col[MrBam_index['Chr']],col[MrBam_index['Start']],col[MrBam_index['Chr']],col[MrBam_index['Start']],col[MrBam_index['End']],col[MrBam_index['Ref']],col[MrBam_index['Alt']]
            FuncRefGene,GeneDetail,ExonicFunc,g1000aug2015_all,snp138,sift=col[MrBam_index['Func.refGene']],col[MrBam_index['GeneDetail.refGene']],col[MrBam_index['ExonicFunc.refGene']],col[MrBam_index['1000g2015aug_all']],col[MrBam_index['snp138']],col[MrBam_index['SIFT_score']]
            if sampletype=='pair':
                tumor_info,normal_info=col[-3],col[-4]
                [gtotal,gsr,gsv,gpercent]=normal_info.split(':')[2:5]+[str(float(normal_info.split(':')[5][:-1])+0.001)]
                [gMultiple_Overlap_Alt,gMultiple_Nonoverlap_Alt,gMultiple_Single_Alt,gOne_Overlap_Alt,gOne_Nonoverlap_Alt,gOne_Single_Alt]=normal_info.split(':')[-1].split(',')[6:12]
                Gunique_alt,Gunique_alt_overlap=str(int(gOne_Overlap_Alt)+int(gOne_Nonoverlap_Alt)+int(gMultiple_Overlap_Alt)+int(gMultiple_Nonoverlap_Alt)),str(int(gOne_Overlap_Alt)+int(gMultiple_Overlap_Alt))
                [genotype,total,sr,sv,percent,sr_for,sr_rev,sv_for_2,sv_rev_2]=[tumor_info.split(':')[0]]+tumor_info.split(':')[2:5]+[tumor_info.split(':')[5][:-1]]+tumor_info.split(':')[6].split(',') 
            elif sampletype=='single':
                tumor_info=col[-1]
                [genotype,total,sr,sv,percent,sr_for,sr_rev,sv_for_2,sv_rev_2]=[tumor_info.split(':')[0]]+[tumor_info.split(':')[2]]+tumor_info.split(':')[4:6]+[tumor_info.split(':')[6][:-1]]+tumor_info.split(':')[10:14]
            [Multiple_Overlap_Alt,Multiple_Nonoverlap_Alt,Multiple_Single_Alt,One_Overlap_Alt,One_Nonoverlap_Alt,One_Single_Alt]=tumor_info.split(':')[-1].split(',')[6:11]+[tumor_info.split(':')[-1].split(',')[11]]
            [unique_alt,unique_alt_overlap]=[str(int(One_Overlap_Alt)+int(One_Nonoverlap_Alt)+int(Multiple_Overlap_Alt)+int(Multiple_Nonoverlap_Alt)),str(int(One_Overlap_Alt)+int(Multiple_Overlap_Alt))]

            VAF_percent=percent
            #====cosmic 信息===============
            if col[MrBam_index['cosmic87']]=='.':
                cosmicid,cosDetails,cosmic='-','-','-'
            else:
                cosmicid,cosDetails=col[MrBam_index['cosmic87']].split(';')[0].split('=')[-1],col[MrBam_index['cosmic87']].split('=')[-1]
                Details=cosDetails.split(',')
                IDs=cosmicid.split(',')
                cosmicid='-'.join(IDs)
                cosDetails='-'.join(Details)
                for count in Details:
                    cosCount+=int(count.split('(')[0])
                L=[]

                for id in IDs:
                    L.append(int(id[4:]))
                    mincosmic=min(L)
                    if str(mincosmic) in id:
                        cosmic=id                
                cosCount=str(cosCount)
            #=============gene,nm,exon,base,aa============
            TRANSCRIPT_FILE='/haplox/tools/clinical/bed/Gene_transcript_new.list'
            TRANSCRIPT_DICT={}
            with open(TRANSCRIPT_FILE,'r')as t:
                for line in t.readlines():
                    transcript_info=line.split('\t')
                    TRANSCRIPT_DICT[transcript_info[0]]=transcript_info[1]
            mutation_details=col[MrBam_index['AAChange.refGene']]
            regGene=col[MrBam_index['Gene.refGene']]
            call=before_filt_mutation(sampletype,FuncRefGene,g1000aug2015_all,regGene,alt,ExonicFunc,percent,GeneDetail)
            if call ==False:
                continue
            gene,nm,exon,base,aa='','','','',''
            if mutation_details !='.' and mutation_details!='UNKNOWN':
                if ',' not in mutation_details and ';' not in mutation_details :
                    [gene,nm,exon,base,aa]=mutation_details.split(':')
                else:
                    mutations=re.split('[, ;]',mutation_details)
                    gene,nm,exon,base,aa=[],[],[],[],[]
                    for mutation in mutations:
                        if 'wholegene' not in mutation:
                            [gene1,nm1,exon1,base1,aa1]=mutation.split(':')
                            if nm1==TRANSCRIPT_DICT[gene1]:
                                gene.append(gene1)
                                nm.append(nm1)
                                exon.append(exon1)
                                base.append(base1)
                                aa.append(aa1)
                    gene,nm,exon,base,aa=','.join(gene),','.join(nm),','.join(exon),','.join(base),','.join(aa)  
            
            elif mutation_details =='.'  and GeneDetail !='.':   
                mutations=re.split('[, ;]',GeneDetail)
                for mutation in mutations:
                    if 'exon' in mutation:
                        [nm1,exon1,base1]=mutation.split(':')
                    else:
                        [nm1,base1]=mutation.split(':')
                    aa1='p.?'
                    gene1=regGene.split(',')
                    for gene2 in gene1:
                        gene2=gene2.split(';')[0]
                        (tmp_error,tmp_data) = subprocess.getstatusoutput('grep "{S1}" {S2}'.format(S1=gene2,S2=TRANSCRIPT_FILE))
                        if nm1==tmp_data.split('\t')[1]:
                            gene,nm,exon,base,aa=gene2,nm1,exon1,base1,aa1
                       
            elif mutation_details =='.'  and GeneDetail =='.':
                if 'MET' in re.split('[, ;]',regGene) and 'intronic' in FuncRefGene.split(';'):
                    gene,nm,exon,base,aa='MET','NM_000245','intron','c.unknow','p.?'
                    if 116411709 <= int(start) <= 116411902:
                        exon='intron13'
                if 'TERT' in re.split('[, ;]',regGene) and 'upstream' in FuncRefGene.split(';'):
                    gene,nm,exon,aa='TERT','NM_198253','upstream','p.?'
                    base='c.'+str((1295104-int(start)))+ref+'>'+alt
                if 'ncRNA_exonic' in re.split('[, ;]',regGene) and 'EGFR-AS1' in re.split('[, ;]',gene) and alt=='TCCAGGAAGCCT':
                    gene,nm,exon,base,aa='EGFR','NM_005228','exon20','c.2290_2291insTCCAGGAAGCCT','p.A763_Y764insFQEA'
            #==========genedetail=========
            if GeneDetail !='.':
                for everydetail in GeneDetail.split(','):
                    if nm in everydetail:
                        GeneDetail=everydetail
          
            #==========整合结果=============
            if [gene,nm,exon,base,aa] != ['','','','','']:
                    if sampletype=='pair' :
                        if pait_after_filt(sv_for_2,sv_rev_2,cosCount,gpercent,percent,sv)==True:
                            MrBam_info.append([dataid,gene,chr,nm,exon,base,aa,VAF_percent,cosmic,cosmicid,chromosome,position,genotype,total,sr,sv,sr_for,sr_rev,sv_for_2,sv_rev_2,Multiple_Overlap_Alt,Multiple_Nonoverlap_Alt,Multiple_Single_Alt,One_Overlap_Alt,One_Nonoverlap_Alt,One_Single_Alt,unique_alt,unique_alt_overlap,Gunique_alt,Gunique_alt_overlap,cosDetails,cosCount,chrs,start,end,ref,alt,FuncRefGene,GeneDetail,ExonicFunc,g1000aug2015_all,snp138,sift,gtotal,gsr,gsv,gpercent,gMultiple_Overlap_Alt,gMultiple_Nonoverlap_Alt,gMultiple_Single_Alt,gOne_Overlap_Alt,gOne_Nonoverlap_Alt,gOne_Single_Alt])
                    elif sampletype=='single':
                        MrBam_info.append([dataid,gene,chr,nm,exon,base,aa,VAF_percent,cosmic,cosmicid,chromosome,position,genotype,total,sr,sv,sr_for,sr_rev,sv_for_2,sv_rev_2,Multiple_Overlap_Alt,Multiple_Nonoverlap_Alt,Multiple_Single_Alt,One_Overlap_Alt,One_Nonoverlap_Alt,One_Single_Alt,unique_alt,unique_alt_overlap,cosDetails,cosCount,chrs,start,end,ref,alt,FuncRefGene,GeneDetail,ExonicFunc,g1000aug2015_all,snp138,sift])
    return MrBam_info

     #=====过滤=======#
def before_filt_mutation(sampletype,FuncRefGene,g1000aug2015_all,gene,alt,ExonicFunc,percent,GeneDetail):
    if  g1000aug2015_all!='.' and float(g1000aug2015_all)>0.001:
        return False 
    if ExonicFunc=='synonymous SNV' or ExonicFunc=='unknow'  or percent =='0':
        return False
    if sampletype=='single' and gene== "MEF2BNB-MEF2B":
        return False
    if sampletype=='pair':
        if gene== "MDC1" or gene== "ZNF717":
            return False
    if ',' in gene or ';' in gene:
        for one in re.split('[, ;]',gene):
            if one in ['UGT1A1','UGT1A10','UT1A3','UGT1A4','UGT1A5','UGT1A6','UGT1A7','UGT1A8','UGT1A9']:
                return False   
    if 'exonic' in FuncRefGene.split(';'):
        return True
    if 'splicing' in FuncRefGene.split(';')  and 'exon' in GeneDetail:
            return True
    if 'UTR5' in FuncRefGene.split(';')  or 'upstream' in FuncRefGene.split(';'):
        if 'TERT' in re.split('[, ;]',gene):  
            return True
    if 'intronic' in FuncRefGene.split(';') and 'MET' in re.split('[, ;]',gene):  
        return True
    if 'ncRNA_exonic' in FuncRefGene.split(';') and 'EGFR-AS1' in re.split('[, ;]',gene) and alt=='TCCAGGAAGCCT':
        return True
    return False
def pait_after_filt(sv_for_2,sv_rev_2,cosCount,gpercent,percent,sv):
    if float(gpercent)>2 or int(sv)<2 or float(percent)<0.05:
        return False
    if int(sv_for_2)>=int(sv_rev_2) and int(sv_rev_2)/(int(sv_for_2)+1e-05)<0.01:
        return False
    if int(sv_for_2)<int(sv_rev_2) and int(sv_for_2)/(int(sv_rev_2)+1e-05)<0.01:
        return False
    if int(cosCount)>=10 and float(percent)/float(gpercent) >2:
        return True
    if int(cosCount)<10 and float(percent)/float(gpercent) >=5:
        return True
    return False
 

    #======输出结果=====#
def get_result(outputdir,result,sample,type):
    with open(outputdir +'/'+sample+'.'+type+'-test-nobias-GB18030.csv','w',encoding='GB18030',newline='') as w:
        writer = csv.writer(w)
        writer.writerows(result)
    os.system('sed -i  "s/\r//" {S1}'.format(S1=outputdir +'/'+sample+'.'+type+'-test-nobias-GB18030.csv'))
    os.system('/haplox/thinker/net/ctDNA/baselineanno4 -anno 1,3,4,5,6,9,10,11 -redis Haplox2022*#@10.6.12.2:6379 -i {S1} -o {S2}'.format(S1=outputdir +'/'+sample+'.'+type+'-test-nobias-GB18030.csv',S2=outputdir +'/'+sample+'.'+type+'-test-nobias-GB18030-baseline.csv'))
    
if __name__ =='__main__':
    main()


'''
def get_result(outputdir,result,sample,type):
    with open(outputdir +'/'+sample+'.'+type+'-nobias-GB18030.csv','w',encoding='GB18030',newline='') as w:
        writer = csv.writer(w)
        writer.writerows(result)
    os.system('sed -i  "s/\r//" {S1}'.format(S1=outputdir +'/'+sample+'.'+type+'-nobias-GB18030.csv'))
    os.system('/haplox/thinker/net/ctDNA/baselineanno4 -anno 1,3,4,5,6,9,10,11 -redis Haplox2022*#@10.6.12.2:6379 -i {S1} -o {S2}'.format(S1=outputdir +'/'+sample+'.'+type+'-nobias-GB18030.csv',S2=outputdir +'/'+sample+'.'+type+'-nobias-GB18030-baseline.csv'))
    
if __name__ =='__main__':
    main()

'''
