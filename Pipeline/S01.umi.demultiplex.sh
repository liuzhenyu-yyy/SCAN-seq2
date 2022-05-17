#!/bin/sh

#######################################################
#
#  Demultiplex barcode for SCAN-Seq2 
#  Author: Liu Zhenyu
#  Latest modification: July 23, 2021
#
#######################################################
fastq=$1
barcode_5=$2
brarode_3=$3

nanoplexer=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/nanoplexer/nanoplexer

mkdir 01.demultiplex

# 5' barcode
$nanoplexer -b ${barcode_5}\
	 -p 01.demultiplex\
	 -t 8\
	 ${fastq}

ls 01.demultiplex | grep Barcode > file.list
cd 01.demultiplex

# 3' barcode
while read fastq
do
cell=${fastq%.*}
Bc=${cell#*e}
mkdir $cell
$nanoplexer -b ../${brarode_3} -p $cell -t 8 ${cell}.fastq
cd $cell
rename Barcode Bc${Bc}_Bc Barcode*
cd ..
done < ../file.list

mv */Bc*fastq ./
cat */unclassified.fastq >> unclassified.fastq
rm -rf Barcode*
rm -rf ../file.list
