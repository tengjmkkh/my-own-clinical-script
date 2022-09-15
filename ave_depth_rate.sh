#!/bin/bash
route=$1
name=`echo ${route##*/}`
data_id=`echo ${name##*_}`
csv_file=$route/BamQC/$name/bam_pileup_QC/*bam_qc.csv
bam_file=$route/$name'_rg.bam'
var=`awk -F "," 'NR==2{print $5}' $csv_file`
target_ave_depth=${var%(X)}
bed=/haplox/tools/clinical/bed/HapOnco680_V2_Boke_primary.bed
samtools depth -b $bed $bam_file >$route/depth.txt
m=`awk 'END{print NR}' $route/depth.txt`
ave=`echo "scale=2;$target_ave_depth/5" |bc`
list=`awk -F '\t' '{print $3}' $route/depth.txt`
count1=0
for i in $list
    do
    if [ `echo "$i > $ave"|bc` -eq 1 ]
        then
        count1=$(($count1+1))
        continue
    fi
done
rate=`awk 'BEGIN{printf "%.2f%\n", ('$count1'/'$m')*100}'`
curl -d "data_id=$data_id" -d "uniformity=$rate" "haplab.haplox.net/api/sequence-qc/uniformity" 

