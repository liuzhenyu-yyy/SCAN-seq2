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
    -L ${database_dir}/Igblast/TRA.bed -b > ${cell}.TRA.bam

$samtools view ../../bam/${cell}_sorted.bam \
    -L ${database_dir}/Igblast/TRB.bed -b > ${cell}.TRB.bam

$samtools fastq ${cell}.TRA.bam  > ${cell}.TRA.all.fastq
$samtools fastq ${cell}.TRB.bam  > ${cell}.TRB.all.fastq

TRA_reads=$(wc -l ${cell}.TRA.all.fastq | awk '{print $1/4}')
TRB_reads=$(wc -l ${cell}.TRB.all.fastq | awk '{print $1/4}')

########################################################
# S02. cluster reads, find representitve reads
########################################################
# TRA
$usearch \
    -cluster_fast ${cell}.TRA.all.fastq \
    -id 0.75 \
    -centroids centroids.fa \
    -sizeout\
    -clusters ./c_
echo > c_1000
cluster_max=$(grep -c '>' c_* | sed 's/:/\t/g' | sort -k 2 -n -r | head -n 1 | awk '{print $1}')

$seqkit head -n 1 ${cluster_max} | awk '{if ($1~/^>/) $1=">TRA_draft"; print $1}' > ${cell}.TRA.draft.fa
cat ${cluster_max} | grep '>' | awk '{gsub(">","",$1); print $0}'> reads.id
TRA_cluster=$(wc -l reads.id | awk '{print $1}')
$seqkit grep -f reads.id ${cell}.TRA.all.fastq > ${cell}.TRA.reads.fastq
rm c_*  ${cell}.TRA.all.fastq reads.id centroids.fa -f

$usearch \
    -cluster_fast ${cell}.TRB.all.fastq \
    -id 0.75 \
    -centroids centroids.fa \
    -sizeout\
    -clusters ./c_
echo > c_1000
cluster_max=$(grep -c '>' c_* | sed 's/:/\t/g' | sort -k 2 -n -r | head -n 1 | awk '{print $1}')

$seqkit head -n 1 ${cluster_max} | awk '{if ($1~/^>/) $1=">TRB_draft"; print $1}' > ${cell}.TRB.draft.fa
cat ${cluster_max} | grep '>' | awk '{gsub(">","",$1); print $0}'> reads.id
TRB_cluster=$(wc -l reads.id | awk '{print $1}')
$seqkit grep -f reads.id ${cell}.TRB.all.fastq > ${cell}.TRB.reads.fastq
rm c_*  ${cell}.TRB.all.fastq reads.id centroids.fa -f

echo -e "${cell}\t${TRA_reads}\t${TRB_reads}\t${TRA_cluster}\t${TRB_cluster}" > ${cell}.vdj.stat

########################################################
# S03. four round of Racon consensus
########################################################
# TRA
mkdir racon_TRA
mv ${cell}.TRA.draft.fa racon_TRA/TRA.draft.fa
mv ${cell}.TRA.reads.fastq racon_TRA/TRA.reads.fastq
cd racon_TRA

DRAFT=TRA.draft.fa
for ROUND in {01..04}; do
  READS2TIGS=reads2contigs_${ROUND}.paf
  NEWDRAFT=racon_${ROUND}.fa

  $minimap2 \
    -x map-ont \
    -t 1 \
    ${DRAFT} \
    TRA.reads.fastq > ${READS2TIGS}

  $racon \
    -w 200 -m 8 \
    -x -6 -g -8  \
    -t 1 \
    -q 7 \
    TRA.reads.fastq \
    ${READS2TIGS} \
    ${DRAFT} > ${NEWDRAFT}

  DRAFT=${NEWDRAFT}
done 2> consensus.log

cd ..
cat racon_TRA/$DRAFT |awk '{if ($1~/^>/) $1=">TRA_racon_contig_1"; print $1}' > ${cell}.TRA.racon.fa

# TRB
mkdir racon_TRB
mv ${cell}.TRB.draft.fa racon_TRB/TRB.draft.fa
mv ${cell}.TRB.reads.fastq racon_TRB/TRB.reads.fastq
cd racon_TRB

DRAFT=TRB.draft.fa
for ROUND in {01..04}; do
  READS2TIGS=reads2contigs_${ROUND}.paf
  NEWDRAFT=racon_${ROUND}.fa

  $minimap2 \
    -x map-ont \
    -t 1 \
    ${DRAFT} \
    TRB.reads.fastq > ${READS2TIGS}

  $racon \
    -w 200 -m 8 \
    -x -6 -g -8  \
    -t 1 \
    -q 7 \
    TRB.reads.fastq \
    ${READS2TIGS} \
    ${DRAFT} > ${NEWDRAFT}

  DRAFT=${NEWDRAFT}
done 2> consensus.log

cd ..
cat racon_TRB/$DRAFT |awk '{if ($1~/^>/) $1=">TRB_racon_contig_1"; print $1}' > ${cell}.TRB.racon.fa

########################################################
# S05. Run Meddaka polish
########################################################
source /gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/miniconda3/etc/profile.d/conda.sh
conda activate medaka

$medaka_consensus \
    -i racon_TRB/TRB.reads.fastq \
    -d ${cell}.TRB.racon.fa \
    -o medaka_TRB -t 1 -m r941_min_high_g360

$medaka_consensus \
    -i racon_TRA/TRA.reads.fastq \
    -d ${cell}.TRA.racon.fa \
    -o medaka_TRA -t 1 -m r941_min_high_g360

conda deactivate

rm *fa.fai *fa.map-ont.mmi
cp medaka_TRA/consensus.fasta ./${cell}.TRA.medaka.fa
cp medaka_TRB/consensus.fasta ./${cell}.TRB.medaka.fa

########################################################
# S04. Run Igblast
########################################################
$igblastn \
  -germline_db_V ${igblast_data}/TRAV \
  -germline_db_J ${igblast_data}/TRAJ \
  -germline_db_D ${igblast_data}/TRBD \
  -organism human \
  -ig_seqtype TCR \
  -query ${cell}.TRA.medaka.fa \
  -auxiliary_data ${igblast_data}/../../optional_file/human_gl.aux \
  -show_translation \
  -outfmt 19  > ${cell}.IgBlast.TRA.tsv

$igblastn \
  -germline_db_V ${igblast_data}/TRBV \
  -germline_db_J ${igblast_data}/TRBJ \
  -germline_db_D ${igblast_data}/TRBD \
  -organism human \
  -ig_seqtype TCR \
  -query ${cell}.TRB.medaka.fa \
  -auxiliary_data ${igblast_data}/../../optional_file/human_gl.aux \
  -show_translation \
  -outfmt 19  > ${cell}.IgBlast.TRB.tsv
