#!/bin/sh

#######################################################
#
#  Mapping to Genome & reference guided assemble
#  Author: Liu Zhenyu
#  Latest modification: Feb 18. 2022
#
#######################################################

cell=$1
ref=$2

# softwares
minimap2=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/minimap2/minimap2
samtools=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/samtools/bin/samtools
stringtie=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/stringtie/stringtie

# directories
Pipeline_dir=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/Pipeline/SCAN-Seq2
database_dir=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/database 
SQANTI3_dir=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/SQANTI3-4.1


########################################################
# S01. Mapping to Genome
########################################################

$minimap2 -t 3\
	 -ax splice\
	 -uf -k14\
	 --secondary=no\
	 ${database_dir}/${ref}/minimap2/*.dna.primary_assembly.mmi \
         ../../02.process/${cell}/${cell}.dedup.fastq.gz	|\
         $samtools view -ShuF 2308 -q 30 -| \
         $samtools sort - -o ${cell}_sorted.bam

$samtools index ${cell}_sorted.bam

########################################################
# S02. Trascriptome Assemble
########################################################
$stringtie ${cell}_sorted.bam \
  -L\
  --rf -v\
  -G ${database_dir}/${ref}/*.101.gtf\
  -l ${cell} \
  -p 4\
  -o ${cell}.stringtie.gff

 cat ${cell}.stringtie.gff | awk '$7!="." {print $0}' > ${cell}.stringtie.filter.gff

########################################################
# S03. SQANTI3  Trascriptome QC
########################################################
source /gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/miniconda3/etc/profile.d/conda.sh
conda activate SQANTI3.env

export PYTHONPATH=$PYTHONPATH:/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/cDNA_Cupcake/sequence/
export PYTHONPATH=$PYTHONPATH:/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/cDNA_Cupcake/ 

/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/miniconda3/envs/SQANTI3.env/bin/python ${SQANTI3_dir}/sqanti3_qc.py \
    --skipORF \
    -t 4 \
    --report pdf \
    -o ${cell}.stringtie \
    ${cell}.stringtie.filter.gff \
    ${database_dir}/${ref}/*.101.gtf \
    ${database_dir}/${ref}/minimap2/*.dna.primary_assembly.fa

/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/miniconda3/bin/python ${SQANTI3_dir}/sqanti3_RulesFilter.py \
    --report pdf \
    ${cell}.stringtie_classification.txt \
    ${cell}.stringtie_corrected.fasta \
    ${cell}.stringtie_corrected.gtf

rm  -rf ${cell}.stringtie_classification.txt ${cell}.stringtie_corrected.fasta ${cell}.stringtie.filter.gff ${cell}.stringtie_corrected.gtf \
    ${cell}.stringtie_corrected.genePred ${cell}.stringtie_corrected.gtf.cds.gff ${cell}.stringtie_junctions.txt \
    ${cell}.stringtie_classification.filtered_lite.fasta ${cell}.stringtie_classification.filtered_lite_junctions.txt \
    ${cell}.stringtie.params.txt GMST RTS ${cell}.stringtie.gff refAnnotation*

########################################################
# S04. Summarize statistics
########################################################
Total_Trans=$(wc -l ${cell}.stringtie_classification.filtered_lite_classification.txt | awk '{print $1-1}')
FSM=$(awk '$6=="full-splice_match" {print $1}' ${cell}.stringtie_classification.filtered_lite_classification.txt | wc -l)
ISM=$(awk '$6=="incomplete-splice_match" {print $1}' ${cell}.stringtie_classification.filtered_lite_classification.txt | wc -l)
NIC=$(awk '$6=="novel_in_catalog" {print $1}' ${cell}.stringtie_classification.filtered_lite_classification.txt | wc -l)
NNC=$(awk '$6=="novel_not_in_catalog" {print $1}' ${cell}.stringtie_classification.filtered_lite_classification.txt | wc -l)
Genic_Genomic=$(awk '$6=="genic" {print $1}' ${cell}.stringtie_classification.filtered_lite_classification.txt | wc -l)
Antisense=$(awk '$6=="antisense" {print $1}' ${cell}.stringtie_classification.filtered_lite_classification.txt | wc -l)
Fusion=$(awk '$6=="fusion" {print $1}' ${cell}.stringtie_classification.filtered_lite_classification.txt | wc -l)
Intergenic=$(awk '$6=="intergenic" {print $1}' ${cell}.stringtie_classification.filtered_lite_classification.txt | wc -l)
Genic_intron=$(awk '$6=="genic_intron" {print $1}' ${cell}.stringtie_classification.filtered_lite_classification.txt | wc -l)

NIC_CJ=$(awk '$6=="novel_in_catalog" && $15=="combination_of_known_junctions" {print $1}' ${cell}.stringtie_classification.filtered_lite_classification.txt | wc -l)
NIC_CS=$(awk '$6=="novel_in_catalog" && $15=="combination_of_known_splicesites" {print $1}' ${cell}.stringtie_classification.filtered_lite_classification.txt | wc -l)
NIC_IR=$(awk '$6=="novel_in_catalog" && $15=="intron_retention" {print $1}' ${cell}.stringtie_classification.filtered_lite_classification.txt | wc -l)


echo -e ${cell}'\t'${Total_Trans}'\t'${FSM}'\t'${ISM}'\t'${NIC}'\t'${NNC}'\t'${Genic_Genomic}'\t'${Antisense}'\t'${Fusion}'\t'${Intergenic}'\t'${Genic_intron}'\t'${NIC_CJ}'\t'${NIC_CS}'\t'${NIC_IR} > ${cell}_Isoform.category.txt
