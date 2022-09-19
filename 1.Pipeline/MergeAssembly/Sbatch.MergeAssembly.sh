#!/bin/sh

Pipeline_dir=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/Pipeline/SCAN-Seq2

while read CL; \
do \
    mkdir ${CL}; \
    cd ${CL}; \
    my_sbatch.l.sh Merge_Assembly_${CL} 4 "${Pipeline_dir}/MergeAssembly.sh ${CL}"; \
    cd ..; \
done < cell.line
