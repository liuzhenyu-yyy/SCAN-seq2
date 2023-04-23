#!/bin/sh

#######################################################
#
#  Process Each Single Cell for SCAN-Seq2 
#  Author: Liu Zhenyu
#  Latest modification: Mar 11. 2022
#
#######################################################

cell=$1
ref=$2

# softwares
salmon=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/salmon-latest_linux_x86_64/bin/salmon
minimap2=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/miniconda3/bin/minimap2
NanoFilt=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/miniconda3/bin/NanoFilt
NanoStat=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/miniconda3/bin/NanoStat
seqkit=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/bin/seqkit
samtools=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/samtools/bin/samtools
cdna_classifier=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/miniconda3/bin/cdna_classifier.py
cutadapt=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/miniconda3/bin/cutadapt

# directories
Pipeline_dir=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/Pipeline/SCAN-Seq2
database_dir=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/database
primer_dir=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/database/HTP_SCANSeq/primers
soft_dir=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/bin
umi_tools=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/miniconda3/bin/umi_tools


########################################################
# S01. QC of Reads
########################################################
$NanoFilt -q 7 -l 500 ../../01.demultiplex/${cell}.fastq > ${cell}_clean.fastq
$NanoStat --fastq ${cell}_clean.fastq -t 4 -o STAT/ -n ${cell}_stat.txt

########################################################
# S02. Trim anchor, identify and orient full-length cDNA 
########################################################
primer_3plus=${cell#*_}

# identify full length; orient(+/-) and make reverse complement of (-); trim primer 
$cdna_classifier -m edlib \
	-b ${primer_dir}/${primer_3plus}.primers.fa \
	-c ${primer_dir}/primer_config.txt \
	-u unclassified.fq -w rescued.fq -A aln_hits.bed \
	-S cdna_classifier_report.tsv -r cdna_classifier_report.pdf \
	-t 4 \
	${cell}_clean.fastq \
	${cell}_full_length.fastq

rm ${cell}_clean.fastq

mv cdna_classifier_report.pdf ${cell}_cdna_classifier_report.pdf
mv cdna_classifier_report.tsv ${cell}_cdna_classifier_report.tsv
rm -f aln_hits.bed rescued.fq unclassified.fq 

# read length < 100 and remove gaps
$seqkit seq -m 108 -g ${cell}_full_length.fastq > ${cell}_full_length_filtered.fastq
rm -f ${cell}_full_length.fastq

# extract umi
$umi_tools extract --bc-pattern=NNNNNNNN \
                   --3prime \
                   --stdin=${cell}_full_length_filtered.fastq \
                   --stdout=${cell}_full_length_filtered.extract.fastq
rm -f ${cell}_full_length_filtered.fastq

# remove extra 8 bp
$cutadapt -u -8 -o ${cell}_full_length_filtered.fastq ${cell}_full_length_filtered.extract.fastq
rm -f ${cell}_full_length_filtered.extract.fastq

# trim poly-A tail
$cutadapt -a "A{10}" -e 0.2 -o ${cell}_full_length_trim_ploya.fastq ${cell}_full_length_filtered.fastq
rm -f ${cell}_full_length_filtered.fastq
QC_Reads=$(cat ${cell}_full_length_trim_ploya.fastq | grep -c runid)

########################################################
# S03. Mapping to cDNA reference and quantification 
########################################################
# mapping to reference cDNA
$minimap2 -t 4 \
	 -ax map-ont \
	 -N 100 \
	 -p 0.99 \
	 ${database_dir}/${ref}/minimap2/*.cdna.all.mmi \
	 ${cell}_full_length_trim_ploya.fastq |\
         $samtools sort - -o ${cell}_cdna.sort.bam

# index bam
$samtools index ${cell}_cdna.sort.bam

# deduplication
$umi_tools dedup \
        --method=directional \
        --edit-distance-threshold=1 \
        --per-gene \
        --per-contig \
        --buffer-whole-contig \
        --stdin=${cell}_cdna.sort.bam \
        --stdout=${cell}_cdna.dedup.bam

# convert deduped bam to fasta
$samtools fastq ${cell}_cdna.dedup.bam > ${cell}.dedup.fastq
rm -f ${cell}_cdna.dedup.bam {cell}_cdna.sort.bam ${cell}_full_length_trim_ploya.fastq

# re-alignment deduped reads
$minimap2 -t 4 \
         -ax map-ont \
         -N 100 \
         -p 0.99 \
         ${database_dir}/${ref}/minimap2/*.cdna.all.mmi \
         ${cell}.dedup.fastq | \
         samtools view -bS > ${cell}_cdna.dedup.bam
gzip ${cell}.dedup.fastq

# quantify umi for each transcript
$salmon quant -p 4 --noErrorModel -l U\
        --minAssignedFrags 1\
        -t ${database_dir}/${ref}/minimap2/*.cdna.all.fa\
        -a ${cell}_cdna.dedup.bam\
        -o ${cell}_cdna_salmon_quant\
        -g ${database_dir}/${ref}/TranscriptID_GeneID_from_cDNA_fasta.txt

rm ${cell}_cdna.dedup.bam

mv ${cell}_cdna_salmon_quant/*.sf ./
mv quant.genes.sf ${cell}_gene_quant.txt
mv quant.sf ${cell}_trans_quant.txt
rm -rf ${cell}_cdna_salmon_quant

########################################################
# S04. Summarize statistics
########################################################
Raw_Reads=$(cat ../../01.demultiplex/${cell}.fastq | grep -c runid)
Q7_Reads=$(cat STAT/${cell}_stat.txt | grep Q7 | awk 'BEGIN{IFS='\t'} {print $2}')
Q7_Percent=$(printf "%.2f" `echo "scale=2;100*${Q7_Reads}/${Raw_Reads}"|bc`)
Q10_Reads=$(cat STAT/${cell}_stat.txt | grep Q10 | awk 'BEGIN{IFS='\t'} {print $2}')
Q10_Percent=$(printf "%.2f" `echo "scale=2;100*${Q10_Reads}/${Raw_Reads}"|bc`)
Q15_Reads=$(cat STAT/${cell}_stat.txt | grep Q15 | awk 'BEGIN{IFS='\t'} {print $2}')
Q15_Percent=$(printf "%.2f" `echo "scale=2;100*${Q15_Reads}/${Raw_Reads}"|bc`)


Mapped_Reads=$(samtools view -F2308 -c ${cell}_cdna.sort.bam)
rm ${cell}_cdna.sort.bam ${cell}_cdna.sort.bam.bai

UMI_count=$(cat *.err | grep "Total # of mapped reads" | awk 'BEGIN{IFS='\t'} {print $7}')
Mapped_percent=$(printf "%.2f" `echo "scale=2;100*${Mapped_Reads}/${QC_Reads}"|bc`)

echo -e ${cell}'\t'${Raw_Reads}'\t'${Q7_Reads}'\t'${Q7_Percent}'\t'${Q10_Reads}'\t'${Q10_Percent}'\t'${Q15_Reads}'\t'${Q15_Percent}'\t'${QC_Reads}'\t'${Mapped_Reads}'\t'${UMI_count}'\t'${Mapped_percent} > ${cell}_ReadCount.txt

mkdir Salmon_Quant
mv ${cell}_gene_quant.txt ${cell}_trans_quant.txt ${cell}_ReadCount.txt Salmon_Quant