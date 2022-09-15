#coding=utf-8
import sys
import csv
import datetime
import os
import subprocess
import time
import requests
import json

  
def dorequest(taskName):
    Authorization ='Bearer eyJhbGciOiJSUzI1NiIsInR5cCI6IkpXVCJ9.eyJlbWFpbCI6Inh1a3hAaGFwbG94LmNvbSIsImV4cCI6MTY2MzE1ODg2NCwiaWF0IjoxNjYzMTE1NjY0LCJpaWQiOjIyLCJpc3MiOiJodHRwczovL2FwaS5oYXB5dW4uY29tL2FjY291bnQvbG9naW4iLCJuYW1lIjoi5b6Q5oG65qyjIiwibmJmIjoxNjYzMTE1NjY0LCJyb2xlIjoiYWRtaW4iLCJyb2xlcyI6Im1hbGwsZW1wbG95ZWUsYmlvaW5mbyIsInN1YiI6NTA3NDAsInVzZXJuYW1lIjoieHVreCJ9.PVdHn2BKKaoivNJMbANoiZ_XWTCof1-eM7vvoVGZThoHkn18jnyNDk8yo39VFkCdqISdR1qP9p1g_d-x3xlAujnZoRqdm5bZJSEtCbzQYd00mxehiucZtcQyL0FE8pSVKVgBciS-pwCqPluYLuo2bMfDUAVSjQaKlM97UGs96kc'
    url = "https://api.hapyun.com/sz-api/tasks?page=1&page_size=10&taskName=" + taskName
    response = requests.request("GET", url, headers={'Authorization': Authorization}, data={})
    if response.json()['error'] == 'crypto/rsa: verification error':
        return -1
    if response.json()['data'][0]['status'] == 'Finished':
        return 1
    elif response.json()['data'][0]['status'] == 'Failed' and response.json()['data'][0]['totalTools'] - response.json()['data'][0]['finishedTools'] == 1:
        return 1
    elif response.json()['data'][0]['status'] == 'Running':
        return 0
    else:
        print(response.json())
        return -1


def download_sh(batch,json_file):
    with open(json_file, 'r') as f:
        stat = json.load(f)
    sample_num = 0 
    while sample_num != len(stat):
        print(sample_num,len(stat))
        for sample in stat:
            if stat[sample] == 'working':
                print(sample,dorequest(sample))
                if dorequest(sample) == 1:
                    if 'cfdna' not in sample:
                        print('nohup sh /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/download.sh &'.format(S1=batch,S2=sample))
                        os.system('nohup sh /haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}/download.sh &'.format(S1=batch,S2=sample))
                    else:
                        print('nohup sh /haplox/rawout/{S1}/{S2}/download.sh &'.format(S1=batch,S2=sample))
                        os.system('nohup sh /haplox/rawout/{S1}/{S2}/download.sh &'.format(S1=batch,S2=sample))
                    print(sample)
                    stat[sample] = 'finish'
                    sample_num += 1
                elif dorequest(sample) == -1:
                    stat[sample] = 'error'
                    sample_num += 1
        time.sleep(300)

    return stat

def write_json(stat,json_file):
    print(json_file)
    with open(json_file, 'w') as f:
        json.dump(stat, f)

def make_sh(batch,upload_csv):
    stat = {}
    with open(upload_csv,'r') as f:
        data = csv.reader(f)
        month = str(datetime.datetime.now().month)
        year = str(datetime.datetime.now().year)
        if len(month) == 1:
            now_date = year + '0' + month
        else:
            now_date = year + month

        for i in data:
            if i[2] == 'nR1':
                continue
            output_sh = ''
            sample = i[11]
            stat[sample] = 'working'
            order = i[10]
            data_id = i[9]
            for gdna in i[2].split('/')[-1].split('_'):
                if gdna.isdigit():
                    gdna_id = gdna
            gdna_sample = i[2].split('/')[-1].split(gdna_id)[0] + gdna_id
            if 'cfdna' in sample:
                basedir = '/haplox/rawout/{S1}/{S2}'.format(S1=batch,S2=sample)
            else:
                basedir = '/haplox/rawout/{S1}/ffpedna_vs_gdna/{S2}'.format(S1=batch,S2=sample)

            os.system('mkdir -p ' + basedir)
            output_sh += 'mkdir -p {S1}/fastp/QC/ \n'.format(S1=basedir)
            output_sh += 'mkdir -p {S1}/MutScan/ \n'.format(S1=basedir)
            output_sh += 'mkdir -p {S1}/cnv/ \n'.format(S1=basedir)
            output_sh += 'mkdir -p {S1}/1p19q/ \n'.format(S1=basedir)
            #===fastpQc====
            output_sh += '#===fastpQc====\n\n'
            output_sh += 'ossutil cp -ruf oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S1}/tmp/fastp-16core64g_node_2/reports/ --include S* {S2}/fastp/QC/ \n\n'.format(S1=sample,S2=basedir,S3=now_date)
            output_sh += 'ossutil cp -ruf oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S1}/tmp/fastp-16core64g_node_1/reports/ --include S* {S2}/fastp/QC/ \n\n'.format(S1=sample,S2=basedir,S3=now_date)
            #===MrBam====
            output_sh += '#===MrBam====\n\n'
            output_sh += 'ossutil cp -ruf oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S1}/tmp/MrBam_pair_node_16/mrbam/{S4}_{S5}.snv_MrBam.txt {S2}\n'.format(S1=sample,S2=basedir,S3=now_date,S4=order,S5=data_id)
            output_sh += 'mv {S1}/{S4}_{S5}.snv_MrBam.txt {S1}/{S4}_{S5}.snv_MrBam.txt\n'.format(S1=basedir,S4=order,S5=data_id)
            output_sh += 'ossutil cp -ruf oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S1}/tmp/MrBam_pair_node_17/mrbam/{S4}_{S5}.indel_MrBam.txt {S2}\n'.format(S1=sample,S2=basedir,S3=now_date,S4=order,S5=data_id)
            #===Mutscan====
            output_sh += '#===Mutscan====\n\n'
            output_sh += 'ossutil cp -ruf oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S1}/tmp/MutScan_node_24/ {S2}/MutScan/\n'.format(S1=sample,S2=basedir,S3=now_date)
            #===cnv====
            output_sh += '#===cnv====\n\n'
            output_sh += 'ossutil cp -ruf oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S1}/tmp/cnv-pair_node_25/ {S2}/cnv/\n'.format(S1=sample,S2=basedir,S3=now_date)
            output_sh += 'mv {S2}/cnv/0001.sample_rg_cnv_result.txt    {S2}/cnv/{S1}_rg_cnv_result.txt\n'.format(S1=sample,S2=basedir)
            output_sh += 'sh /haplox/users/liaowt/Script/Cnv/getCnv.sh {S2} {S1}\n'.format(S1=sample,S2=basedir)
            output_sh += 'sh /haplox/users/liaowt/Script/Cnv/get_1p19q/get_1p19q_cnr_clinic_pair.sh {S2} {S1}\n'.format(S1=sample,S2=basedir)
            output_sh += 'sh /haplox/users/liaowt/Script/circos/circos_cmd.sh {S2} {S1}\n'.format(S1=sample,S2=basedir)
            output_sh += 'Rscript /haplox/users/liaowt/Script/Cnv/get_genes_cns_broad/get_genes_cns_broad.R {S2}/cnv/0001.sample_rg.cns {S2}/cnv/{S1}_get_genes_cns_broad.txt {S3}\n'.format(S1=sample,S2=basedir,S3=data_id)
            output_sh += 'Rscript /haplox/users/huang/myR/pqCnv.R -i {S2}/cnv/0001.sample_rg.cns\n'.format(S1=sample,S2=basedir)
            output_sh += 'mv {S2}/cnv/0001.sample_rg.txt {S2}/cnv/{S1}_pqcnv.txt\n'.format(S1=sample,S2=basedir)
            if 'HP-IO-55' in sample:
                output_sh += 'Rscript /haplox/users/liaowt/Script/Cnv/last_HP-IO-55_cnv.R {S1}/\n'.format(S1=basedir)
            #---------------fusion-------------------------------------------------------
            output_sh += '#===fusion====\n\n'
            output_sh += 'mkdir -p {S2}/fusionscan/\n'.format(S1=sample,S2=basedir)
            output_sh += 'ossutil cp -ruf oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S1}/tmp/genefuse_huang_test_node_22/ {S2}/fusionscan/\n'.format(S1=sample,S2=basedir,S3=now_date)

            output_sh += 'mv {S2}/fusionscan/0001.genefuse_huang_test_out_html_output_22.html {S2}/fusionscan/{S1}_fusion.html\n'.format(S1=sample,S2=basedir)
            output_sh += 'mv {S2}/fusionscan/0001.genefuse_huang_test_out_json_output_22.json {S2}/fusionscan/{S1}_fusion.json\n'.format(S1=sample,S2=basedir)

            output_sh += 'mkdir -p {S2}/fusion/\n'.format(S1=sample,S2=basedir)
            output_sh += 'ossutil cp -ruf oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S1}/tmp/factera_fusion_node_35/fusion/0001.sentieon-bwa-gencore-pileup-no-umi_cfdna_rgfactera.fusions.txt_output_4.factera.fusions.txt {S2}/fusion/\n'.format(S1=sample,S2=basedir,S3=now_date)
            output_sh += 'ossutil cp oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S1}/tmp/factera_fusion_node_35/fusion/0001.sentieon-bwa-gencore-pileup-no-umi_cfdna_rgfactera.discordantpair.depth.txt_output_4.factera.discordantpair.depth.txt {S2}/fusion/\n'.format(S1=sample,S2=basedir,S3=now_date)
            output_sh += 'ossutil cp oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S1}/tmp/factera_fusion_node_35/fusion/0001.sentieon-bwa-gencore-pileup-no-umi_cfdna_rgfactera.discordantpair.details.txt_output_4.factera.discordantpair.details.txt {S2}/fusion/\n'.format(S1=sample,S2=basedir,S3=now_date)

            output_sh += 'mv {S2}/fusion/0001.sentieon-bwa-gencore-pileup-no-umi_cfdna_rgfactera.fusions.txt_output_4.factera.fusions.txt {S2}/fusion/{S1}.fusions.txt\n'.format(S1=sample,S2=basedir)
            output_sh += 'mv {S2}/fusion/0001.sentieon-bwa-gencore-pileup-no-umi_cfdna_rgfactera.discordantpair.depth.txt_output_4.factera.discordantpair.depth.txt  {S2}/fusion/{S1}.discordantpair.depth.txt\n'.format(S1=sample,S2=basedir)
            output_sh += 'mv {S2}/fusion/0001.sentieon-bwa-gencore-pileup-no-umi_cfdna_rgfactera.discordantpair.details.txt_output_4.factera.discordantpair.details.txt {S2}/fusion/{S1}.discordantpair.details.txt\n'.format(S1=sample,S2=basedir)
            
            #-------------------------germline----------------------------------
            output_sh += '#===germline====\n\n'
            output_sh += 'mkdir -p {S2}/germline/result/\n'.format(S1=sample,S2=basedir)
            output_sh += 'ossutil cp -ruf oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S1}/tmp/germline-sentieon-nofastp_node_36/result/ --exclude *bam {S2}/germline/result/\n'.format(S1=sample,S2=basedir,S3=now_date)

            output_sh += 'mv {S2}/germline/result/sample.information.txt {S2}/germline/result/{S1}.information.txt\n'.format(S1=gdna_sample,S2=basedir)
            output_sh += 'mv {S2}/germline/result/sample.Target_451.txt  {S2}/germline/result/{S1}.Target_451.txt\n'.format(S1=gdna_sample,S2=basedir)

            output_sh += 'perl /haplox/users/liaowt/Script/germline/new_germline/PathCall.pl -i {S2}/germline/result/sample.filter.hg19_multianno.txt -g /haplox/users/liaowt/Script/germline/new_germline/Heart_Cancer_Incidental_v2.txt -o {S2}/germline/result/{S1}.germline.txt\n'.format(S1=gdna_sample,S2=basedir)


            #-----------------------Visualmsi-----------------------------------------------
            output_sh += '#===msi====\n\n'
            output_sh += 'mkdir -p {S2}/Visual_msi/\n'.format(S1=sample,S2=basedir)
            output_sh += 'ossutil cp -ruf oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S1}/tmp/Visual_msi_node_41/ {S2}/Visual_msi/\n'.format(S1=sample,S2=basedir,S3=now_date)
            output_sh += 'mv {S2}/Visual_msi/0001.Visual_msi_html_output_41.html {S2}/Visual_msi/{S1}.html\n'.format(S1=sample,S2=basedir)
            output_sh += 'mv {S2}/Visual_msi/0001.Visual_msi_json_output_41.json {S2}/Visual_msi/{S1}.json\n'.format(S1=sample,S2=basedir)
            output_sh += 'perl /haplox/users/liaowt/Script/Visual_msi/visual_msi_curl.pl {S2}/Visual_msi/{S1}.json >{S2}/Visual_msi/{S1}_visual_msi_curl.txt\n'.format(S1=sample,S2=basedir)
            output_sh += 'perl /haplox/users/liaowt/Script/Visual_msi/visual_msi_stat.pl {S2}/Visual_msi/{S1}.json >{S2}/Visual_msi/{S1}_visual_msi_stat.txt\n'.format(S1=sample,S2=basedir)
             #-----------------------Virus-----------------------------------------------
            output_sh += 'mkdir -p {S1}/virus/\n'.format(S2=sample,S1=basedir,S3=now_date)
            output_sh += 'ossutil cp -ru oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S2}/tmp/virus_node_37/0001.virus_resulttxt_cfdna_output_37.txt {S1}/virus/\n'.format(S2=sample,S1=basedir,S3=now_date)
            output_sh += 'ossutil cp -ru oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S2}/tmp/virus_node_37/virusBaseline_curl.txt {S1}/virus/\n'.format(S2=sample,S1=basedir,S3=now_date)
            output_sh += 'mv {S1}/virus/0001.virus_resulttxt_cfdna_output_37.txt  {S1}/virus/{S2}_virus_result.txt\n'.format(S2=sample,S1=basedir,S3=now_date)
            output_sh += 'mv {S1}/virus/virusBaseline_curl.txt {S1}/virus/{S2}_virus_result_curl.txt\n'.format(S2=sample,S1=basedir,S3=now_date)

            output_sh += 'ossutil cp -ru oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S2}/tmp/virus_node_38/0001.virus_resulttxt_gdna_output_38.txt {S1}/virus/\n'.format(S2=sample,S1=basedir,S3=now_date)
            output_sh += 'ossutil cp -ru oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S2}/tmp/virus_node_38/virusBaseline_curl.txt {S1}/virus/\n'.format(S2=sample,S1=basedir,S3=now_date)
            output_sh += 'mv {S1}/virus/0001.virus_resulttxt_gdna_output_38.txt  {S1}/virus/{S2}_virus_result.txt\n'.format(S2=sample,S1=basedir,S3=now_date)
            output_sh += 'mv {S1}/virus/virusBaseline_curl.txt {S1}/virus/{S2}_virus_result_curl.txt\n'.format(S2=sample,S1=basedir,S3=now_date)

            output_sh += 'perl /haplox/users/liaowt/Script/virus/virus_tumor_normal.pl {S1}/virus/{S2}_virus_result.txt {S1}/virus/{S2}_virus_result.txt > {S1}/virus/virus_tumor_normal_result.txt\n'.format(S2=sample,S1=basedir,S3=now_date)

            #------------------------BamQc----------------------------------------------
            output_sh += '#===BamQc====\n\n'
            output_sh += 'mkdir -p {S2}/BamQC/{S1}\n'.format(S1=sample,S2=basedir)
            output_sh += 'ossutil cp -ruf oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S1}/tmp/rgbam_captureByBase_depth_maprate_node_43/ {S2}/BamQC/{S1}\n'.format(S1=sample,S2=basedir,S3=now_date)
            output_sh += 'mv {S2}/BamQC/{S1}/depth/0001.sentieon-bwa-gencore-pileup-no-umi_cfdna_rgbam_output_4.bam.depth {S2}/BamQC/{S1}/depth/{S1}_bam.depth\n'.format(S1=sample,S2=basedir)
            output_sh += 'mv {S2}/BamQC/{S1}/depth/0001.sentieon-bwa-gencore-pileup-no-umi_cfdna_rgbam_output_4.bam_ave_depth.stat {S2}/BamQC/{S1}/depth/{S1}_bam_ave_depth.stat\n'.format(S1=sample,S2=basedir)
            output_sh += 'mv {S2}/BamQC/{S1}/capture/0001.sentieon-bwa-gencore-pileup-no-umi_cfdna_rgbam_output_4.bam_capture_stat.txt {S2}/BamQC/{S1}/capture/{S1}_capture_stat.txt\n'.format(S1=sample,S2=basedir)
            output_sh += 'mv {S2}/BamQC/{S1}/bam_pileup_QC/0001.sentieon-bwa-gencore-pileup-no-umi_cfdna_rgbam_output_4.bam_qc.csv {S2}/BamQC/{S1}/bam_pileup_QC/{S1}_bam_qc.csv\n'.format(S1=sample,S2=basedir)
            output_sh += 'mv {S2}/BamQC/{S1}/bam_pileup_QC/0001.sentieon-bwa-gencore-pileup-no-umi_cfdna_rgbam_output_4.bam_temp_len.png {S2}/BamQC/{S1}/bam_pileup_QC/{S1}_bam_temp_len.png\n'.format(S1=sample,S2=basedir)

            output_sh += 'mkdir -p {S2}/BamQC/{S1}\n'.format(S1=gdna_sample,S2=basedir)
            output_sh += 'ossutil cp -ruf oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S1}/tmp/rgbam_captureByBase_depth_maprate_node_42/ {S2}/BamQC/{S4}\n'.format(S1=sample,S2=basedir,S3=now_date,S4=gdna_sample)
            output_sh += 'mv {S2}/BamQC/{S1}/depth/0001.sentieon-bwa-gencore-pileup-no-umi_gdna_rgbam_output_3.bam.depth {S2}/BamQC/{S1}/depth/{S1}_bam.depth\n'.format(S1=gdna_sample,S2=basedir)
            output_sh += 'mv {S2}/BamQC/{S1}/depth/0001.sentieon-bwa-gencore-pileup-no-umi_gdna_rgbam_output_3.bam_ave_depth.stat {S2}/BamQC/{S1}/depth/{S1}_bam_ave_depth.stat\n'.format(S1=gdna_sample,S2=basedir)
            output_sh += 'mv {S2}/BamQC/{S1}/capture/0001.sentieon-bwa-gencore-pileup-no-umi_gdna_rgbam_output_3.bam_capture_stat.txt {S2}/BamQC/{S1}/capture/{S1}_capture_stat.txt\n'.format(S1=gdna_sample,S2=basedir)
            output_sh += 'mv {S2}/BamQC/{S1}/bam_pileup_QC/0001.sentieon-bwa-gencore-pileup-no-umi_gdna_rgbam_output_3.bam_qc.csv {S2}/BamQC/{S1}/bam_pileup_QC/{S1}_bam_qc.csv\n'.format(S1=gdna_sample,S2=basedir)
            output_sh += 'mv {S2}/BamQC/{S1}/bam_pileup_QC/0001.sentieon-bwa-gencore-pileup-no-umi_gdna_rgbam_output_3.bam_temp_len.png {S2}/BamQC/{S1}/bam_pileup_QC/{S1}_bam_temp_len.png\n'.format(S1=gdna_sample,S2=basedir)
            output_sh += "sed -i 's/{S3}/{S4}/' {S2}/BamQC/{S1}/bam_pileup_QC/{S1}_bam_qc.csv\n".format(S1=gdna_sample,S2=basedir,S3=data_id,S4=gdna_id)
            #----------------------------------------------------------------------
            output_sh += '#===bam====\n\n'
            output_sh += '#ossutil cp oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S1}/tmp/sentieon-bwa-gencore-pileup-no-umi_node_3/0001.sentieon-bwa-gencore-pileup-no-umi_gdna_rgbam_output_3.bam {S2}\n'.format(S1=sample,S2=basedir,S3=now_date)
            output_sh += '#mv {S2}/0001.sentieon-bwa-gencore-pileup-no-umi_gdna_rgbam_output_3.bam {S2}/{S1}_rg.bam\n'.format(S1=gdna_sample,S2=basedir)
            
            output_sh += '#ossutil cp oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S1}/tmp/sentieon-bwa-gencore-pileup-no-umi_node_3/0001.sentieon-bwa-gencore-pileup-no-umi_gdna_rgbam_output_3.bam.bai {S2}\n'.format(S1=sample,S2=basedir,S3=now_date)
            output_sh += '#mv {S2}/0001.sentieon-bwa-gencore-pileup-no-umi_gdna_rgbam_output_3.bam.bai  {S2}/{S1}_rg.bai\n'.format(S1=gdna_sample,S2=basedir)
            output_sh += 'ossutil cp oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S1}/tmp/sentieon-bwa-gencore-pileup-no-umi_node_4/0001.sentieon-bwa-gencore-pileup-no-umi_cfdna_rgbam_output_4.bam {S2}\n'.format(S1=sample,S2=basedir,S3=now_date)
            output_sh += 'mv {S2}/0001.sentieon-bwa-gencore-pileup-no-umi_cfdna_rgbam_output_4.bam {S2}/{S1}_rg.bam\n'.format(S1=sample,S2=basedir)
            output_sh += 'ossutil cp oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S1}/tmp/sentieon-bwa-gencore-pileup-no-umi_node_4/0001.sentieon-bwa-gencore-pileup-no-umi_cfdna_rgbam_output_4.bam.bai {S2}/{S1}_rg.bam.bai\n'.format(S1=sample,S2=basedir,S3=now_date)
            output_sh += '#mv {S2}/0001.sentieon-bwa-gencore-pileup-no-umi_cfdna_rgbam_output_4.bam.bai {S2}/{S1}_rg.bam.bai \n'.format(S1=sample,S2=basedir)
            output_sh += 'ossutil cp oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S1}/tmp/sentieon-bwa-gencore-pileup-no-umi_node_4/gencore.json {S2} \n'.format(S1=sample,S2=basedir,S3=now_date)
            #----------------------------------------------------------------------
            output_sh += '#===annovar====\n\n'
            output_sh += 'ossutil cp -ruf oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S1}/tmp/varscan-annovar-pair_node_18/ {S2}\n'.format(S1=sample,S2=basedir,S3=now_date)

            output_sh += 'mv {S2}/sample.snp.vcf {S2}/{S1}_varscan2.snp.vcf\n'.format(S1=sample,S2=basedir)
            output_sh += 'mv {S2}/sample.indel.vcf {S2}/{S1}_varscan2.indel.vcf\n'.format(S1=sample,S2=basedir)
            output_sh += 'mv {S2}/sample_snv_annovar.avinput {S2}/{S1}_snv_annovar.avinput\n'.format(S1=sample,S2=basedir)
            output_sh += 'mv {S2}/sample_indel_annovar.avinput {S2}/{S1}_indel_annovar.avinput\n'.format(S1=sample,S2=basedir)
            output_sh += 'mv {S2}/sample_snv_annovar.hg19_multianno.vcf {S2}/{S1}_snv_annovar.hg19_multianno.vcf\n'.format(S1=sample,S2=basedir)
            output_sh += 'mv {S2}/sample_indel_annovar.hg19_multianno.vcf {S2}/{S1}_indel_annovar.hg19_multianno.vcf\n'.format(S1=sample,S2=basedir)
            output_sh += 'mv {S2}/0001.varscan-annovar-pair_snv_txt_output_18.txt {S2}/{S1}_snv_annovar.hg19_multianno.txt\n'.format(S1=sample,S2=basedir)
            output_sh += 'mv {S2}/0001.varscan-annovar-pair_indel_txt_output_18.txt {S2}/{S1}_indel_annovar.hg19_multianno.txt\n'.format(S1=sample,S2=basedir)
            #----------------------------------------------------------------------
            output_sh += '#===germline_cnv====\n\n'
            output_sh += 'mkdir -p {S2}/germline_cnv\n'.format(S1=sample,S2=basedir)
            output_sh += 'ossutil cp -ruf oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S1}/tmp/germline_cnv_node_44/ {S2}/germline_cnv/\n'.format(S1=sample,S2=basedir,S3=now_date)
            output_sh += 'mv {S2}/germline_cnv/0001.sample_rg.cnr {S2}/germline_cnv/{S1}_rg.cnr\n'.format(S1=gdna_sample,S2=basedir)
            output_sh += 'mv {S2}/germline_cnv/0001.sample_rg.cns {S2}/germline_cnv/{S1}_rg.cns\n'.format(S1=gdna_sample,S2=basedir)
            output_sh += 'mv {S2}/germline_cnv/0001.sample_rg_cnv_result.txt {S2}/germline_cnv/{S1}_rg_cnv_result.txt\n'.format(S1=gdna_sample,S2=basedir)
            output_sh += '#grep -E "APC|MUTYH|STK11|BMPR1A|SMAD4|MLH1|PMS2|PMS2|MSH6|EPCAM|TP53|SDHA|SDHB|SDHC|SDHD|BRCA1|BRCA2|CDH1|TSC2|TSC1|TFE3|MSH2" {S1}/germline_cnv/{S2}_rg.cnr > {S1}/germline_cnv/{S2}_important_genes_rg_cnr.txt\n'.format(S2=sample,S1=basedir,S3=now_date)
            output_sh += '#sh /haplox/users/liaowt/Script/Cnv/get_germline_Cnv.sh {S2} {S1}\n'.format(S1=gdna_sample,S2=basedir)

            #----------------------------------------------------------------------
            output_sh += '#===median_depth====\n\n'
            output_sh += 'ossutil cp -ruf oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S1}/tmp/median_depth_node_46/0001.median_depth_median_depth_stat_output_46_gdna.txt.txt {S2}/BamQC/{S1}/depth/\n'.format(S1=sample,S2=basedir,S3=now_date)
            output_sh += 'mv {S2}/BamQC/{S1}/depth/0001.median_depth_median_depth_stat_output_46_gdna.txt.txt {S2}/BamQC/{S1}/depth/{S1}_bam_median_depth.txt\n'.format(S1=sample,S2=basedir)

            output_sh += 'ossutil cp -ruf oss://sz-hapres/haplox/hapyun/{S3}/pair_IDT_clinical_v2_{S1}/tmp/median_depth_node_45/0001.median_depth_median_depth_stat_output_45_cfdna.txt.txt  {S2}/BamQC/{S4}/depth/\n'.format(S1=sample,S2=basedir,S3=now_date,S4=gdna_sample)
            output_sh += 'mv {S2}/BamQC/{S1}/depth/0001.median_depth_median_depth_stat_output_45_cfdna.txt.txt {S2}/BamQC/{S1}/depth/{S1}_bam_median_depth.txt\n'.format(S1=gdna_sample,S2=basedir)

            #----------------------------------------------------------------------
            output_sh += 'Rscript /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/nobias_getMrBam_txt_V3.R {S2} {S3}_{S4}\n'.format(S1=sample,S2=basedir,S3=order,S4=data_id)
            output_sh += 'Rscript /haplox/users/huang/mypy/data-analysis/ctdna_exome_pipeline/idSNP.R {S2} {S1}\n'.format(S1=sample,S2=basedir)
            output_sh += 'sh /haplox/users/liaowt/Script/Cnv/get_1p19q/get_1p19q_idSNP_clinic_pair.sh {S2} {S1}\n'.format(S1=sample,S2=basedir)
            output_sh += 'perl /thinker/net/ctDNA/check.pl {S2}/{S1}_snv_annovar.hg19_multianno.txt > {S2}/{S1}_check.pl.txt\n'.format(S1=sample,S2=basedir)
         
            output_sh += '#mkdir -p {S1}/tsv/\n'.format(S2=order,S1=basedir,S3=data_id)
            output_sh += '/haplox/thinker/net/tools/extract_vcf_nofilter {S1}/tsv/ {S1}/{S2}_{S3}.snv_MrBam.txt\n'.format(S2=order,S1=basedir,S3=data_id)
            output_sh += '/haplox/thinker/net/tools/extract_vcf_nofilter {S1}/tsv/ {S1}/{S2}_{S3}.indel_MrBam.txt\n'.format(S2=order,S1=basedir,S3=data_id)

            output_sh += '/haplox/thinker/net/ctDNA/samplemutationimport -i {S1}/tsv{S1}/{S2}_{S3}.snv_MrBam.tsv\n'.format(S2=order,S1=basedir,S3=data_id)
            output_sh += '/haplox/thinker/net/ctDNA/samplemutationimport -i {S1}/tsv{S1}/{S2}_{S3}.indel_MrBam.tsv\n'.format(S2=order,S1=basedir,S3=data_id)
            #----------------------------------------------------------------------
            if 'ffpedna_vs_gdna' in basedir: 
                #----------------------FineMSI-------------------------------------------
                output_sh += 'python3 /haplox/users/chenyr/bin/MSI/src/msi.py msi {S2}/{S1}_rg.bam -o {S2}/{S1}.mono6.v3 -d /haplox/users/chenyr/bin/MSI/sup/605.msi.sites.v3.txt --baseline /haplox/users/chenyr/bin/MSI/sup/msi.mono6_reference_gdna_1_100.baseline \n'.format(S1=sample,S2=basedir)
                output_sh += 'Rscript /haplox/users/huang/myR/convertFineMSI.R --incsv={S2}/{S1}.mono6.v3.total.msi.csv \n'.format(S1=sample,S2=basedir)
                output_sh += 'python3 /haplox/users/tengjm/script/cnv_updata.py -input {S2}/cnv/{S1}_rg_cnv.genes188in680.csv -bed /haplox/tools/clinical/bed/chrs_gene_cytoband_hg19.txt -pqcnv {S2}/cnv/{S1}_pqcnv.txt -output {S2}/cnv/ -s ffpe -n {S3} \n'.format(S1=sample,S2=basedir,S3=data_id)
            else:
                #----------------------FineBMSI-------------------------------------------
                output_sh += 'python3 /haplox/users/yangbo/bmsi/bmsi_test.py {S2}/{S1}_rg.bam  \n'.format(S1=sample,S2=basedir)
                output_sh += 'python3 /haplox/users/tengjm/script/cnv_updata.py -input {S2}/cnv/{S1}_rg_cnv.genes188in680.csv -bed /haplox/tools/clinical/bed/chrs_gene_cytoband_hg19.txt -pqcnv {S2}/cnv/{S1}_pqcnv.txt -output {S2}/cnv/ -s cfdna -n {S3} \n'.format(S1=sample,S2=basedir,S3=data_id)
           

            #----------------------------------------------------------------------
            output_sh += 'sh /haplox/users/liaowt/Script/germline/new_Chem/NewPan_chem.sh {S2} {S1}\n'.format(S1=gdna_sample,S2=basedir)
            output_sh += 'sh /haplox/users/liaowt/Script/germline/new_germline/New_germline.sh {S2} {S1}\n'.format(S1=gdna_sample,S2=basedir)
            output_sh += 'Rscript /haplox/users/liaowt/Script/germline/pair_86gene_chem_germline.R {S2} {S1}\n'.format(S1=gdna_sample,S2=basedir)
            output_sh += '#sh /haplox/tools/clinical/analysis/ave_depth_rate.sh {S1} \n'.format(S1=basedir)
            #----------------------------------------------------------------------
            output_sh += '#===upload====\n\n'
            output_sh += 'curl haplab.haplox.net/api/report/median-depth?data_id={S3} -F import_file=@{S2}/BamQC/{S1}/depth/{S1}_bam_median_depth.txt\n'.format(S1=sample,S2=basedir,S3=data_id)
            output_sh += 'curl haplab.haplox.net/api/report/median-depth?data_id={S3} -F import_file=@{S2}/BamQC/{S1}/depth/{S1}_bam_median_depth.txt\n'.format(S1=gdna_sample,S2=basedir,S3=gdna_id)
            output_sh += 'curl haplab.haplox.net/api/report/depth-new -F "import_file=@{S2}/BamQC/{S1}/bam_pileup_QC/{S1}_bam_qc.csv"\n'.format(S1=sample,S2=basedir)
            output_sh += 'curl haplab.haplox.net/api/report/depth-new -F "import_file=@{S2}/BamQC/{S1}/bam_pileup_QC/{S1}_bam_qc.csv"\n'.format(S1=gdna_sample,S2=basedir)
            output_sh += 'curl haplab.haplox.net/api/report/fragment/{S3} -F "picture=@{S2}/BamQC/{S1}/bam_pileup_QC/{S1}_bam_temp_len.png"\n'.format(S1=sample,S2=basedir,S3=data_id)
            output_sh += 'curl haplab.haplox.net/api/report/fragment/{S3} -F "picture=@{S2}/BamQC/{S1}/bam_pileup_QC/{S1}_bam_temp_len.png"\n'.format(S1=gdna_sample,S2=basedir,S3=gdna_id)
            output_sh += 'curl 192.168.1.10/api/report/hrr?data_id={S3} -F "import_file=@{S2}/germline/result/{S1}_trans.germline_DDR.cancer.txt"\n'.format(S1=gdna_sample,S2=basedir,S3=data_id)
            output_sh += 'curl haplab.haplox.net/api/report/ngs-cnv -F "import_file=@{S2}/cnv/{S1}_get_genes_cns_broad.txt"\n'.format(S1=sample,S2=basedir)
            output_sh += 'curl haplab.haplox.net/api/report/ngs-msi?data_id={S1} -F "import_file=@{S2}/Visual_msi/{S3}_visual_msi_curl.txt"\n'.format(S3=sample,S2=basedir,S1=data_id)
            output_sh += 'curl 192.168.1.10/api/report/msi-stat?data_id={S1} -F "import_file=@{S2}/Visual_msi/{S3}_visual_msi_stat.txt"\n'.format(S3=sample,S2=basedir,S1=data_id)
            output_sh +='curl haplab.haplox.net/api/report/msi-html/{S1} -F "html=@{S2}/Visual_msi/{S3}.html"\n'.format(S3=sample,S2=basedir,S1=data_id)
            output_sh += 'curl http://haplab.haplox.net/api/report/ngs-finemsi?data_id={S3}  -F "import_file=@{S2}/{S1}.mono6.v3.total.msi.txt" \n'.format(S1=sample,S2=basedir,S3=data_id)
            output_sh += 'curl haplab.haplox.net/api/report/glioma?data_id={S3} -F "import_file=@{S2}/cnv/{S1}_rg_cnv.genes12.csv"\n'.format(S1=sample,S2=basedir,S3=data_id)
            output_sh += 'curl haplab.haplox.net/api/report/glioma-files/circos/{S3} -F "file=@{S2}/cnv/{S1}_circos.png" \n'.format(S1=sample,S2=basedir,S3=data_id)
            output_sh += 'curl haplab.haplox.net/api/report/glioma-files/heatmap/{S3} -F "file=@{S2}/cnv/{S1}_heatmap.pdf" \n'.format(S1=sample,S2=basedir,S3=data_id)
            output_sh += 'curl haplab.haplox.net/api/report/ngs-pqcnv?data_id={S3}  -F "import_file=@{S2}/cnv/{S1}_pqcnv.txt" \n'.format(S1=sample,S2=basedir,S3=data_id)
            output_sh += 'python3 /haplox/users/yangbo/script/pollution_qc.py {S2}/{S1}_snv-idSNP-filter-genotype_1_1.csv {S3} \n'.format(S1=sample,S2=basedir,S3=data_id) #上传污染质控
            #----------------------------------------------------------------------
            output_sh += 'perl /haplox/users/liaowt/Script/Cnv/cnv_gender/cnv_gender.pl {S2}/cnv/{S1}_rg_cnv.chrX_genes.csv > {S2}/cnv/{S1}_cnv_gender_judge.txt\n'.format(S1=sample,S2=basedir)
             #----------------------------------------------------------------------
            output_sh += 'Rscript /haplox/users/liaowt/Script/Cnv/cnv_gender/sex_warning.R {S2}/ 220221_A00250_0097_AH357TDSX3\n'.format(S2=basedir)
            output_sh += 'grep -Ew "ALK|FGFR1|NAB2|NTRK3|RARA|BRAF|FGFR2|NCOA4|RET|CCDC6|FGFR3|NPM1|PAX3|ROS1|CD74|FUS|NR4A3|PAX7|SLC34A2|DDIT3|KIF5B|NRG1|PAX8|SS18|EML4|KMT2A|NTRK1|PDGFB|TFE3|EWSR1|MYC|NTRK2|PRKACA|TMPRSS2|CAMTA1|ERG|KRAS|NRAS|TP53" {S1}/fusion/{S2}.fusions.txt >{S1}/fusion/{S2}.factera_important_605.fusions.txt \n'.format(S2=sample,S1=basedir)
            output_sh +='python /haplox/users/yangbo/futionbase.py -f {S1}/fusionscan/{S2}_fusion.json -b {S1}/{S2}_rg.bam \n'.format(S2=sample,S1=basedir)
        
            #----------------------------------------------------------------------

            #------------------------------------------------------------------------
            output_sh += 'chmod -R 777 {S2}/\n'.format(S2=basedir)
            output_sh += "echo 'Hapyun pipeline FINISH !'"

            with open(basedir + '/download.sh','w') as w:
                w.writelines(output_sh)

    return stat

def main():
    batch = sys.argv[1]
    upload_csv = sys.argv[2]
    json_file = sys.argv[3]
    if os.path.isfile(json_file):
        stat = download_sh(batch,json_file)
    else:
        stat = make_sh(batch,upload_csv)
    write_json(stat,json_file)
    

if __name__ == "__main__":
    main()




