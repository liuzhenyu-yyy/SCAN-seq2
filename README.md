# SCAN-seq2

`SCAN-seq2` is a high-throughput, high sensitivity full-length single-cell RNA-seq method with Oxford Nanopore Technology (ONT). Our results proved  SCAN-seq2 as a flexible single-cell full-length RNA-seq methods, with potiential implements in isoform quantification, transcriptome assemble, and Immune repertoire profiling.

<p align="center">
    <img  height="600" src="./Protocol.jpg">
</p>

For more detail about `SCAN-seq2`, please check our [publication](https://www.nature.com/articles/s41421-022-00500-4).


This repository provide source code for SCAN-seq2 data processing and downstream analysis.

## Repository Layout: 
The repository are structured as followed:


```bash
.
├── 1.Data_Processing # Pipeline for SCAN-seq2 data processing
│   ├── S01.umi.demultiplex.sh  ## Demultiplex ONT reads
│   ├── S02.umi.process.one.sh  ## QC, deduplication and trancriptome quantification
│   ├── S03.umi.Summarize.sh    ## Summarize statistics, merge UMI count matrix
│   ├── S04.umi.assemble.one.sh ## Reference-guided trancriptome assemble
│   ├── S05.umi.VDJ.BCR.one.sh  ## Immune repertoire profiling of immunoglobulin
│   ├── S05.umi.VDJ.TCR.one.sh  ## Immune repertoire profiling of T cell receptor
│   ├──merge_counts_tsvs.R      ## Merge UMI count matrix
│   └──MergeAssembly
│       ├── MergeAssembly.all.sh ## Merge cell-line assemblies
│       ├── MergeAssembly.sh     ## Merge single-cell assemblies
│       └── Sbatch.MergeAssembly.sh
│
├── 2.Dowstream_Libraries # Downstream analyis of each separate SCAN-seq library
│   ├── 9CL.R      ## 9CL library
│   ├── 4CL.R      ## 4CL library
│   ├── 9CL_Mix.R  ## 9CL_Mix library
│   ├── UMI_100.R  ## UMI_100 library
│   ├── UMI_200.R  ## UMI_200 library
│   └── IGG.R      ## IGG6, IGG24 and IGG48 library
│
├── 3.Downstream_Merge # Downstream analysis of all libraries
│   ├── Part1_Technical_Performance 
│   │   ├── Technical_Performance.R  ## Technical performance of SCAN-seq2
│   │   └── Pseudogene.R             ## Systematic evaluation of pseudogene expression
│   ├── Part2_Transcriptome_Assembly
│   │   ├── BCR_TCR_gtf              ## Annotation for human immunoglobulin and T cell receptor genes
│   │   │   ├── Hg38.IGH.gtf
│   │   │   ├── Hg38.IGK.gtf
│   │   │   ├── Hg38.IGL.gtf
│   │   │   ├── Hg38.TRA.gtf
│   │   │   └── Hg38.TRB.gtf
│   │   ├── ModifGTF.R               ## Modify gtf files for IGV visualization
│   │   ├── Transcriptome_Assembly.R ## Analysis of trancriptome assemblies
│   │   ├── merge_human_6CL.SQANTI3_SQANTI3_report.pdf
│   │   ├── merge_human_6CL.SQANTI3_classification.filtered_lite_SQANTI3_report.pdf
│   │   ├── merge_mouse_3CL.SQANTI3_SQANTI3_report.pdf
│   │   └── merge_mouse_3CL.SQANTI3_classification.filtered_lite_SQANTI3_report.pdf
│   ├── Part3_Isoginkgetin_Treatment
│   │   └── IGGTreatment.R           ## IGG treatment of Hela and HepG2 cell lines1
│   └── Part4_Revise                 ## Additional code during manuscript revision
│
├── primers # Barcode and primer sequences
│
└── README.md
```

## Pipelines for SCAN-seq2 data processing

We provided schematics for the data processing pipelines of SCAN-seq2:

### Demultiplexing and raw read processing:
This step includes identifing barcode sequences in raw ONT reads, trimming adaptors, quality control and deduplication.

```mermaid
%%{init: {"theme": "neutral", 'themeVariables': { "fontSize": "30px","fontFamily": "Arial"}}}%%
graph LR
raw("Raw reads</br>(84,801,426)"):::merge --nanoplexer--> demult("Demultiplexed reads</br>74,861 per cell"):::sc --NanoFilt<br/>Pychopper<br/>cutadapt--> QC(QC reads</br>66,606 per cell):::sc --minimap2<br/>UMI-tools dedup--> 
dedup(deduped reads</br>43,480 per cell):::sc

classDef merge fill:#bebada,stroke:#000000;
classDef sc fill:#8dd3c7,stroke:#000000;
```

### Transcriptome assembly and quantification:
The main processing step for each single cell, where deduped reads were mapped to transcriptome reference for isoform-level quantification, as well as genome reference for transcriptome assembly. Single-cell assemblies were merged in a hierarchical manner, generating a merged transcriptome assembly.

```mermaid
%%{init: {"theme": "neutral", 'themeVariables': { "fontSize": "30px","fontFamily": "Arial"}}}%%
graph TD
dedup("Deduped reads</br>43,480 per cell</br>(fastq)"):::sc
subgraph Transcriptome Mapping
Align_cDNA("Transcriptome alignments</br>(bam)"):::sc --Salmon quant--> Quant("Gene/Isoform Expression Quantification</br>(tsv)"):::sc --merge_counts_tsvs.R--> Matrix("UMI count matrix</br>cell-by-gene/isoform</br>(tsv)"):::merge 
end
subgraph Genome Mapping
Align_Genome("Genome alignments</br>(bam)"):::sc -- StringTie --> sc_Asm("Raw single-cell assembly<br/>4559 genes, 5364 isoforms per cell</br>(gff)"):::sc -- SQANTI3</br>sqanti3_qc.py</br>sqanti3_RulesFilter.py --> sc_Asm_2("Filterd single-cell assembly</br>3685 genes, 4406 isoforms per cell</br>(gff)"):::sc -- TAMA merge --> cl_Asm("Cell-line assembly</br>(bed12)"):::cl -- TAMA merge --> merge_Asm("Raw merged assembly</br>human: 10429 genes, 40712 isoforms</br>mouse: 8923 genes, 23906 isoforms</br>(gff)"):::merge -- SQANTI3</br>sqanti3_qc.py</br>sqanti3_RulesFilter.py -->  merge_Asm_2("Filtered merged assembly</br>human: 10428 genes, 40676 isoforms</br>mouse: 8922 genes, 23894 isoforms</br>(gff)"):::merge 
end
dedup --minimap2 -ax splice</br>Genome mapping--> Align_Genome
dedup --minimap2 -ax map-ont</br>Transcripome mapping--> Align_cDNA
subgraph Downstream Analysis
sc_Asm_2 --> sc_stat("Statistics of single-cell assemblies"):::down
merge_Asm_2 --> merge_stat("Statistics of merged assembly"):::down
Matrix --> NGS("Compare with NGS methods"):::down
Matrix --> DR_cluster("Dimensional Reduction,<br/>Clustering"):::down
Matrix --> DEG("DEG and DTU"):::down
end
classDef merge fill:#bebada,stroke:#000000;
classDef cl fill:#fb8072,stroke:#000000;
classDef sc fill:#8dd3c7,stroke:#00000;
classDef down fill:#b3de69,stroke:#000000;
```

### Immune repertoire profiling:
This step generated polished transcript assembly for IGH/IGL/IGK/TRA/TRB transcripts, and identify V(D)J elements for each single cell.

```mermaid
%%{init: {"theme": "neutral", 'themeVariables': { "fontSize": "30px","fontFamily": "Arial"}}}%%
graph TD

Align_Genome("Genome alignments</br>(bam)"):::sc -- samtools --> reads("IGH/IGL/IGK/TRA/TRB reads<br/>(fastq)"):::sc --usearch--> cluster("Reads clusters</br>(cluster fastq)"):::sc--> large_clustrer("Largest cluster</br>(fastq)"):::sc--centroid--> centroid("centroid reads"):::sc
large_clustrer --> Other("Other reads in the same cluster"):::sc

subgraph Main 
centroid  --racon---> pol1(Polished sequence I):::sc
pol1 --another 3 rounds of racon--->pol4(Polished sequence IV):::sc -- Medaka ----> con("Consensus sequence</br>(fasta)"):::sc
end

subgraph Other reads
Other --racon---> pol1
Other -- another 3 rounds of racon-->pol4
Other -- Medaka --> con("Consensus sequence</br>(fasta)")
end

con --IgBlast--> VDJ("V(D)J calls</br>(tsv)"):::sc

classDef merge fill:#bebada,stroke:#000000;
classDef cl fill:#fb8072,stroke:#000000;
classDef sc fill:#8dd3c7,stroke:#00000;
classDef down fill:#b3de69,stroke:#000000;
```

## Citation 
Please cite the following article if you use SCAN-seq2 in your research:

> Liao, Y., Liu, Z., Zhang, Y. *et al.* High-throughput and high-sensitivity full-length single-cell RNA-seq analysis on third-generation sequencing platform. *Cell Discov* **9**, 5 (2023). https://doi.org/10.1038/s41421-022-00500-4

## Contact
If you have any question please write issue for this repository or contact me directly at liuzhenyu@pku.edu.cn.
