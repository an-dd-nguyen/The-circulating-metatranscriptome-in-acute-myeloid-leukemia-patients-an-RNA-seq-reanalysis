library("org.Hs.eg.db", character.only = TRUE)
library(clusterProfiler)
library(enrichplot)
library(fgsea)
library(DOSE)
library(tidyverse)

#https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#pubmed-trend-of-enriched-terms
pathways = read.gmt("../../Data/h.all.v2025.1.Hs.symbols.gmt")

#Faecabacterium

rankings_F = results_batch_F$logFC
names(rankings_F) = df$display_label

rankings_F = rankings_F[!is.na(names(rankings_F))]

rankings_F = sort(rankings_F, decreasing = T)

#rankings = rankings[abs(rankings) >=1]

edo = GSEA(rankings_F, pvalueCutoff = 0.1, TERM2GENE = pathways)


dotplot(edo)
#Enterococcus

rankings_E = results_batch_E$logFC
names(rankings_E) = df$display_label

rankings_E = rankings_E[!is.na(names(rankings_E))]

rankings_E = sort(rankings_E, decreasing = T)

#rankings = rankings[abs(rankings) >=1]

edo1 = GSEA(rankings_E, pvalueCutoff = 0.1, TERM2GENE = pathways)






