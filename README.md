# SCAN-seq2
SCAN-seq2  is a high-throughput, high sensitivity full-length single-cell RNA-seq method by third-generation sequencing. 

This repository provide code for SCAN-seq2 data processing and downstream analysis.

# Flowchart for SCAN-seq2 data processing:

- Demultiplexing and raw read processing:
```mermaid
graph LR
raw("Raw Reads</br>(84,801,426)"):::merge --nanoplexer--> demult("demultiplexed Reads</br>74,861 per cell"):::sc --NanoFilt<br/>Pychopper<br/>cutadapt--> QC(QC reads</br>66,606 per cell):::sc --minimap2<br/>UMI-tools dedup--> 
dedup(deduped reads</br>43,480 per cell):::sc

classDef merge fill:#bebada,stroke:#000000;
classDef sc fill:#8dd3c7,stroke:#000000;
```

- Transcriptome assembly and quantification
```mermaid
%%{init: {"theme": "default", 'themeVariables': { "fontSize": "30px","fontFamily": "Arial"}}}%%
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
sc_Asm_2 --> sc_stat("Statistics of single-cell assemblies</br>Fig.S5a"):::down
merge_Asm_2 --> merge_stat("Statistics of merged assembly</br>Fig.S5c and Fig.S5e"):::down
Matrix --> NGS("Compare with NGS methods"):::down
Matrix --> DR_cluster("Dimensional Reduction,<br/>Clustering"):::down
Matrix --> DEG("DEG and DTU"):::down
end
classDef merge fill:#bebada,stroke:#000000;
classDef cl fill:#fb8072,stroke:#000000;
classDef sc fill:#8dd3c7,stroke:#00000;
classDef down fill:#b3de69,stroke:#000000;
```
Legend:
```mermaid
graph TD

sc("Single-cell level"):::sc
cl("cell-line level"):::cl
merged("all cells merged"):::merge
down("Dowstream analysis"):::down

classDef merge fill:#bebada,stroke:#000000;
classDef cl fill:#fb8072,stroke:#000000;
classDef sc fill:#8dd3c7,stroke:#000000;
classDef down fill:#b3de69,stroke:#000000;
```