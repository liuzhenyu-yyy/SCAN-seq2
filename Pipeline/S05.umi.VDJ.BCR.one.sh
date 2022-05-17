#!/bin/sh

cell=$1

export IGDATA=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/ncbi-igblast-1.17.1

# softwares
minimap2=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/minimap2/minimap2
samtools=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/samtools/bin/samtools
usearch=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/usearch/usearch11.0.667_i86linux32
seqkit=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/bin/seqkit
racon=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/miniconda3/bin/racon
medaka_consensus=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/miniconda3/envs/medaka/bin/medaka_consensus
igblastn=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/ncbi-igblast-1.17.1/bin/igblastn

# directories
Pipeline_dir=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/Pipeline/SCAN-Seq2
database_dir=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/database/hg38/
igblast_data=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/ncbi-igblast-1.17.1/database/human_germline


########################################################
# S01. Subset reads from locus
########################################################
$samtools view ../../bam/${cell}_sorted.bam \
    -L ${database_dir}/Igblast/IGL.bed -b > ${cell}.IGL.bam

$samtools view ../../bam/${cell}_sorted.bam \
    -L ${database_dir}/Igblast/IGH.bed -b > ${cell}.IGH.bam

$samtools fastq ${cell}.IGL.bam  > ${cell}.IGL.all.fastq
$samtools fastq ${cell}.IGH.bam  > ${cell}.IGH.all.fastq

IGL_reads=$(wc -l ${cell}.IGL.all.fastq | awk '{print $1/4}')
IGH_reads=$(wc -l ${cell}.IGH.all.fastq | awk '{print $1/4}')

########################################################
# S02. cluster reads, find representitve reads
########################################################
# IGL
$usearch \
    -cluster_fast ${cell}.IGL.all.fastq \
    -id 0.75 \
    -centroids centroids.fa \
    -sizeout\
    -clusters ./c_
echo > c_1000
cluster_max=$(grep -c '>' c_* | sed 's/:/\t/g' | sort -k 2 -n -r | head -n 1 | awk '{print $1}')

$seqkit head -n 1 ${cluster_max} | awk '{if ($1~/^>/) $1=">IGL_draft"; print $1}' > ${cell}.IGL.draft.fa
cat ${cluster_max} | grep '>' | awk '{gsub(">","",$1); print $0}'> reads.id
IGL_cluster=$(wc -l reads.id | awk '{print $1}')
$seqkit grep -f reads.id ${cell}.IGL.all.fastq > ${cell}.IGL.reads.fastq
rm c_*  ${cell}.IGL.all.fastq reads.id centroids.fa -f

$usearch \
    -cluster_fast ${cell}.IGH.all.fastq \
    -id 0.75 \
    -centroids centroids.fa \
    -sizeout\
    -clusters ./c_
echo > c_1000
cluster_max=$(grep -c '>' c_* | sed 's/:/\t/g' | sort -k 2 -n -r | head -n 1 | awk '{print $1}')

$seqkit head -n 1 ${cluster_max} | awk '{if ($1~/^>/) $1=">IGH_draft"; print $1}' > ${cell}.IGH.draft.fa
cat ${cluster_max} | grep '>' | awk '{gsub(">","",$1); print $0}'> reads.id
IGH_cluster=$(wc -l reads.id | awk '{print $1}')
$seqkit grep -f reads.id ${cell}.IGH.all.fastq > ${cell}.IGH.reads.fastq
rm c_*  ${cell}.IGH.all.fastq reads.id centroids.fa -f

echo -e "${cell}\t${IGL_reads}\t${IGH_reads}\t${IGL_cluster}\t${IGH_cluster}" > ${cell}.vdj.stat

########################################################
# S03. four round of Racon consensus
########################################################
# IGL
mkdir racon_IGL
mv ${cell}.IGL.draft.fa racon_IGL/IGL.draft.fa
mv ${cell}.IGL.reads.fastq racon_IGL/IGL.reads.fastq
cd racon_IGL

DRAFT=IGL.draft.fa
for ROUND in {01..04}; do
  READS2TIGS=reads2contigs_${ROUND}.paf
  NEWDRAFT=racon_${ROUND}.fa

  $minimap2 \
    -x map-ont \
    -t 1 \
    ${DRAFT} \
    IGL.reads.fastq > ${READS2TIGS}

  $racon \
    -w 200 -m 8 \
    -x -6 -g -8  \
    -t 1 \
    -q 7 \
    IGL.reads.fastq \
    ${READS2TIGS} \
    ${DRAFT} > ${NEWDRAFT}

  DRAFT=${NEWDRAFT}
done 2> consensus.log

cd ..
cat racon_IGL/$DRAFT |awk '{if ($1~/^>/) $1=">IGL_racon_contig_1"; print $1}' > ${cell}.IGL.racon.fa

# IGH
mkdir racon_IGH
mv ${cell}.IGH.draft.fa racon_IGH/IGH.draft.fa
mv ${cell}.IGH.reads.fastq racon_IGH/IGH.reads.fastq
cd racon_IGH

DRAFT=IGH.draft.fa
for ROUND in {01..04}; do
  READS2TIGS=reads2contigs_${ROUND}.paf
  NEWDRAFT=racon_${ROUND}.fa

  $minimap2 \
    -x map-ont \
    -t 1 \
    ${DRAFT} \
    IGH.reads.fastq > ${READS2TIGS}

  $racon \
    -w 200 -m 8 \
    -x -6 -g -8  \
    -t 1 \
    -q 7 \
    IGH.reads.fastq \
    ${READS2TIGS} \
    ${DRAFT} > ${NEWDRAFT}

  DRAFT=${NEWDRAFT}
done 2> consensus.log

cd ..
cat racon_IGH/$DRAFT |awk '{if ($1~/^>/) $1=">IGH_racon_contig_1"; print $1}' > ${cell}.IGH.racon.fa

########################################################
# S05. Run Meddaka polish
########################################################
source /gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/miniconda3/etc/profile.d/conda.sh
conda activate medaka

$medaka_consensus \
    -i racon_IGH/IGH.reads.fastq \
    -d ${cell}.IGH.racon.fa \
    -o medaka_IGH -t 1 -m r941_min_high_g360

$medaka_consensus \
    -i racon_IGL/IGL.reads.fastq \
    -d ${cell}.IGL.racon.fa \
    -o medaka_IGL -t 1 -m r941_min_high_g360

conda deactivate

rm *fa.fai *fa.map-ont.mmi
cp medaka_IGL/consensus.fasta ./${cell}.IGL.medaka.fa
cp medaka_IGH/consensus.fasta ./${cell}.IGH.medaka.fa

########################################################
# S04. Run Igblast
########################################################
$igblastn \
  -germline_db_V ${igblast_data}/IGLV \
  -germline_db_J ${igblast_data}/IGLJ \
  -germline_db_D ${igblast_data}/IGHD \
  -organism human \
  -query ${cell}.IGL.medaka.fa \
  -auxiliary_data ${igblast_data}/../../optional_file/human_gl.aux \
  -show_translation \
  -outfmt 19  > ${cell}.IgBlast.IGL.tsv

$igblastn \
  -germline_db_V ${igblast_data}/IGHV \
  -germline_db_J ${igblast_data}/IGHJ \
  -germline_db_D ${igblast_data}/IGHD \
  -organism human \
  -query ${cell}.IGH.medaka.fa \
  -auxiliary_data ${igblast_data}/../../optional_file/human_gl.aux \
  -show_translation \
  -outfmt 19  > ${cell}.IgBlast.IGH.tsv
