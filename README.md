# The-circulating-metatranscriptome-in-acute-myeloid-leukemia-patients-an-RNA-seq-reanalysis
Specific software and package version\
Bowtie2/2.5.1-GCC-11.3.0\
KrakenUniq v1.0.4\
R version 4.4.0 (2024-04-24 ucrt) -- "Puppy Cup"\
CIBERSORTx version 1.0 https://cibersortx.stanford.edu/ \
edgeR version 4.2.2\
clusterProfiler version 4.12.6\


## Order of script execution\
1. KrakenUniq Profiling\
For every sequencing input bam file, 01_Host_Depletion_And_KrakenUniq_Profiling.sh was to be executed. The expected output is a .report file of KrakenUniq microbial profiling of the samples. All report files for all samples were to be concatenated and filtered into Genus and Phylum level. In the Data/ directory, we provided the report for all genus/phylum, as well as normalized abundance of the top 100 genera for all samples.\

2. Metatranscriptome Profiling\
02_BeatAML_Microbial_Landscape.R is called to generate phylum distribution and PCoA visualization of samples.\
03_BeatAML_Genus_Analysis.R is called to generate Bray-Curtis dissimilarity between AML pairs vs Healthy-AML pairs, samples alpha diversity, differential abundance of genera in Healthy vs AML samples, top 30 genera in term of abundance/proportion in patients samples.\

3. Statistical Analysis\
04_Wrapper_Function.R is called first to initiate the wrapper function\
05_BeatAML_Clinical_Statistics.R is called for statistical analysis between microbial abundance versus patients clinical information/drug response information. Full patient clinical information can be obtained from the BeatAML paper: Tyner JW, Tognon CE, Bottomly D, et al. Functional genomic landscape of acute myeloid leukaemia. Nature. 2018;562(7728):526-531. doi:10.1038/s41586-018-0623-z\

5. Host Expression Analysis
06_BeatAML_EdgeR.R and 07_BeatAML_clusterprofiler.R is meant to be run sequentially, which will process host expression data and perform GSEA pathway enrichment analysis between high-low metatranscriptome load samples. 07_BeatAML_clusterprofiler.R also analyze samples' cell proportion and perform Wilcoxon test between high-low load samples.\
Similarly, 08 and 09 is meant to be run sequentially, which output GSEA hallmark pathways that are overexpressed in samples with low-butyrate producing genera compare to samples with high-butyrate producing genera. 

