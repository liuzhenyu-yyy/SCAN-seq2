#!/bin/sh

#######################################################
#
#  Merge umi count matrix, generate statistics
#  Author: Liu Zhenyu
#  Latest modification: Aug 24, 2021
#
#######################################################

Pipeline_dir=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/Pipeline/SCAN-Seq2

mkdir 03.summary

# sammarize statistics
echo -e "cell\tRaw_Reads\tQ7_Reads\tQ7_Percent\tQ10_Reads\tQ10_Percent\tQ15_Reads\tQ15_Percent\tQC_Reads\tMapped_Reads\tUMI_count\tMapped_percent" > 03.summary/summary.txt
cat 02.process/B*/Salmon_Quant/*_ReadCount.txt >> 03.summary/summary.txt

my_sbatch.l.sh S03.SCAN-Seq2_Summarize 4 "${Pipeline_dir}/merge_counts_tsvs.R sample.list"

rm -f 01.demultiplex/*fastq

