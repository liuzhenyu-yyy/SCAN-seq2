#!/bin/sh

CL=$1

Pipeline_dir=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/Pipeline/SCAN-Seq2

gtfToGenePred=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/UCSC/gtfToGenePred
genePredToBed=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/UCSC/genePredToBed
bedToGenePred=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/UCSC/bedToGenePred
genePredToGtf=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/UCSC/genePredToGtf

python2=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/miniconda3/envs/python2/bin/python2
tama_merge=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/tama/tama_merge.py
SQANTI3_dir=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/SQANTI3-4.1



########################################################
# S01. convert gtf to bed12
########################################################
mkdir bed
while read cell;
do \
    $gtfToGenePred ../../../batch/9CL_1000/04.assemble/${cell}/${cell}.stringtie_classification.filtered_lite.gtf ${cell}.GenePred; \
    $genePredToBed ${cell}.GenePred ${cell}.bed12; \
    cat ${cell}.bed12 | awk 'BEGIN{FS=OFS="\t"} {split($4,a,"."); $4=a[1]"."a[2]";"$4; $1="chr"$1; print $0}' > 9CL.${cell}.stringtie.filtered.bed; \
    rm ${cell}.GenePred ${cell}.bed12; \
done < ../../0.sample.list/sample.list.9CL.${CL}

while read cell;
do \
    $gtfToGenePred ../../../batch/9Mix_1000/04.assemble/${cell}/${cell}.stringtie_classification.filtered_lite.gtf ${cell}.GenePred; \
    $genePredToBed ${cell}.GenePred ${cell}.bed12; \
    cat ${cell}.bed12 | awk 'BEGIN{FS=OFS="\t"} {split($4,a,"."); $4=a[1]"."a[2]";"$4; $1="chr"$1; print $0}' > 9Mix.${cell}.stringtie.filtered.bed; \
    rm ${cell}.GenePred ${cell}.bed12; \
done < ../../0.sample.list/sample.list.9Mix.${CL}

if [ -f "../../0.sample.list/sample.list.4CL.${CL}" ]; then
    while read cell;
    do \
        $gtfToGenePred ../../../batch/M_4CL/04.assemble/${cell}/${cell}.stringtie_classification.filtered_lite.gtf ${cell}.GenePred; \
        $genePredToBed ${cell}.GenePred ${cell}.bed12; \
        cat ${cell}.bed12 | awk 'BEGIN{FS=OFS="\t"} {split($4,a,"."); $4=a[1]"."a[2]";"$4; $1="chr"$1; print $0}' > 4CL.${cell}.stringtie.filtered.bed; \
        rm ${cell}.GenePred ${cell}.bed12; \
    done < ../../0.sample.list/sample.list.4CL.${CL}
fi

########################################################
# S02. generate annotation file
########################################################
if [ -f "temp" ]; then
    rm temp
fi

ls *.bed | awk '{gsub(".stringtie.filtered.bed","",$1); print $1}' > sample.list

while read cell; 
do \
echo -e "capped\t1,1,1" >> temp; \
done < sample.list

ls *.bed > temp2

paste temp2 temp sample.list > file.list

rm temp2 temp

########################################################
# S03. run tama merge 
########################################################
$python2  ${tama_merge} -f file.list -p ${CL} -a 100 -m 20 -z 50 -d merge_dup

mv *stringtie.filtered.bed bed/

$bedToGenePred ${CL}.bed ${CL}.GenePred 
$genePredToGtf "file" ${CL}.GenePred ${CL}.gtf

rm ${CL}.GenePred