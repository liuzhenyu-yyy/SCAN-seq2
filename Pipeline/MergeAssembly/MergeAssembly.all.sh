#!/bin/sh

Pipeline_dir=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/Pipeline/SCAN-Seq2
database_dir=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/database

gtfToGenePred=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/UCSC/gtfToGenePred
genePredToBed=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/UCSC/genePredToBed
bedToGenePred=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/UCSC/bedToGenePred
genePredToGtf=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/UCSC/genePredToGtf

python2=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/miniconda3/envs/python2/bin/python2
tama_merge=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/tama/tama_merge.py
SQANTI3_dir=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/SQANTI3-4.1

########################################################
# S01. generate annotation file
########################################################
if [ -f "temp" ]; then
    rm temp
fi

ls bed/*.bed | awk '{gsub("bed/|.c10.bed",""); print $1}' > sample.list

while read cell; 
do \
echo -e "capped\t1,1,1" >> temp; \
done < sample.list

ls /gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/project/HTP_SCANSeq/assemble/2.all/bed/*.bed  > temp2
paste temp2 temp sample.list > file.list
rm temp2 temp

cat sample.list | grep -E 'AtT20|MEF|3T3' > sample.list.mm10
cat sample.list | grep -E -v 'AtT20|MEF|3T3' > sample.list.hg38
cat file.list | grep -E 'AtT20|MEF|3T3' > file.list.mm10
cat file.list | grep -E -v 'AtT20|MEF|3T3' > file.list.hg38

rm sample.list*

########################################################
# S02. run tama merge 
########################################################
$python2  ${tama_merge} -f file.list.mm10 -p merge_mouse_3CL -a 10 -m 10 -z 10 -d merge_dup
$python2  ${tama_merge} -f file.list.hg38 -p merge_human_6CL -a 10 -m 10 -z 10 -d merge_dup

$bedToGenePred merge_mouse_3CL.bed merge_mouse_3CL.raw.GenePred 
cat merge_mouse_3CL.raw.GenePred  | awk 'BEGIN{FS=OFS="\t"} {gsub("^.+?;","",$1); gsub("chr","",$2); print $0}' > merge_mouse_3CL.GenePred
$genePredToGtf "file" merge_mouse_3CL.GenePred merge_mouse_3CL.raw.gtf
cat merge_mouse_3CL.raw.gtf | awk 'BEGIN{FS=OFS="\t"} {split($9,a,";"); gene=a[1]; gsub("\\..+?$","\"",gene); gsub(a[1],gene,$9); print $0}' > merge_mouse_3CL.gtf
rm -f merge_mouse_3CL.raw.GenePred merge_mouse_3CL.GenePred merge_mouse_3CL.raw.gtf

$bedToGenePred merge_human_6CL.bed merge_human_6CL.raw.GenePred 
cat merge_human_6CL.raw.GenePred  | awk 'BEGIN{FS=OFS="\t"} {gsub("^.+?;","",$1); gsub("chr","",$2); print $0}' > merge_human_6CL.GenePred
$genePredToGtf "file" merge_human_6CL.GenePred merge_human_6CL.raw.gtf
cat merge_human_6CL.raw.gtf | awk 'BEGIN{FS=OFS="\t"} {split($9,a,";"); gene=a[1]; gsub("\\..+?$","\"",gene); gsub(a[1],gene,$9); print $0}' > merge_human_6CL.gtf
rm -f merge_human_6CL.raw.GenePred merge_human_6CL.GenePred merge_human_6CL.raw.gtf

mkdir mouse 
mkdir human
mv merge_mouse_3CL* mouse
mv merge_human_6CL* human

########################################################
# S03. run SQANTI QC
########################################################
source /gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/miniconda3/etc/profile.d/conda.sh
conda activate SQANTI3.env

export PYTHONPATH=$PYTHONPATH:/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/cDNA_Cupcake/sequence/
export PYTHONPATH=$PYTHONPATH:/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/cDNA_Cupcake/ 

# mm10
cd mouse
/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/miniconda3/envs/SQANTI3.env/bin/python ${SQANTI3_dir}/sqanti3_qc.py \
    --skipORF \
    -t 4 \
    --report pdf \
    -o merge_mouse_3CL.SQANTI3 \
    merge_mouse_3CL.gtf \
    ${database_dir}/mm10/Mus_musculus.GRCm38.101.gtf \
    ${database_dir}/mm10/minimap2/Mus_musculus.GRCm38.dna.primary_assembly.fa

/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/miniconda3/bin/python ${SQANTI3_dir}/sqanti3_RulesFilter.py \
    --report pdf \
    merge_mouse_3CL.SQANTI3_classification.txt \
    merge_mouse_3CL.SQANTI3_corrected.fasta \
    merge_mouse_3CL.SQANTI3_corrected.gtf
cd ..

# hg38
cd human
/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/miniconda3/envs/SQANTI3.env/bin/python ${SQANTI3_dir}/sqanti3_qc.py \
    --skipORF \
    -t 4 \
    --report pdf \
    -o merge_human_6CL.SQANTI3 \
    merge_human_6CL.gtf \
    ${database_dir}/hg38/Homo_sapiens.GRCh38.101.gtf \
    ${database_dir}/hg38/minimap2/Homo_sapiens.GRCh38.dna.primary_assembly.fa

/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/miniconda3/bin/python ${SQANTI3_dir}/sqanti3_RulesFilter.py \
    --report pdf \
    merge_human_6CL.SQANTI3_classification.txt \
    merge_human_6CL.SQANTI3_corrected.fasta \
    merge_human_6CL.SQANTI3_corrected.gtf
cd ..
