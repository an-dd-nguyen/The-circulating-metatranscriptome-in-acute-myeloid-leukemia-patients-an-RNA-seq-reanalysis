library("org.Hs.eg.db", character.only = TRUE)
library(clusterProfiler)
library(enrichplot)
library(fgsea)
library(DOSE)
library(tidyverse)

#https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#pubmed-trend-of-enriched-terms
pathways1 = read.gmt("../Downloads/h.all.v2025.1.Hs.symbols.gmt")
pathways2 = read.gmt("../Downloads/GENTLES_LEUKEMIC_STEM_CELL_UP.v2025.1.Hs.gmt")
pathways3 = read.gmt("../Downloads/GENTLES_LEUKEMIC_STEM_CELL_DN.v2025.1.Hs.gmt")

pathways4 = read.gmt("../Downloads/c2.all.v2025.1.Hs.symbols.gmt")
pathways4 = pathways4[grep("LEUKE", pathways4$term), ]

pathways = rbind(pathways1, pathways2, pathways3)

#Faecabacterium

rankings_F = results_batch_F$logFC
names(rankings_F) = df$display_label

rankings_F = rankings_F[!is.na(names(rankings_F))]

rankings_F = sort(rankings_F, decreasing = T)

#rankings = rankings[abs(rankings) >=1]

edo = GSEA(rankings_F, pvalueCutoff = 0.1, TERM2GENE = pathways)

edo1 = GSEA(rankings, pvalueCutoff = 0.1, TERM2GENE = pathways4)

dotplot(edo)
#Enterococcus

rankings_E = results_batch_E$logFC
names(rankings_E) = df$display_label

rankings_E = rankings_E[!is.na(names(rankings_E))]

rankings_E = sort(rankings_E, decreasing = T)

#rankings = rankings[abs(rankings) >=1]

edo1 = GSEA(rankings_E, pvalueCutoff = 0.1, TERM2GENE = pathways)


rankings2 = results_2$logFC
names(rankings2) = df$display_label

rankings2 = rankings2[!is.na(names(rankings2))]

rankings2 = sort(rankings2, decreasing = T)

edo2 = GSEA(rankings, pvalueCutoff = 0.1, TERM2GENE = pathways)
edo3 = GSEA(rankings, pvalueCutoff = 0.1, TERM2GENE = pathways4)


pathways5 = readRDS("../Downloads/c_Bacteria_pathways.rds")

rankings1 = results_batch$logFC
names(rankings1) = df$display_label

term2gene = read.delim("../Downloads/Term2Gene_Entrez.txt",
                       header = F, sep = " ")
term2gene = term2gene[, -1]
term2gene = term2gene[-1, ]
colnames(term2gene) = c("TermID", "GeneID")


entrezID = mapIds(org.Hs.eg.db, names(rankings1), 
                  column = "ENTREZID", 
                  keytype = "SYMBOL")

#rankings3 = rankings2

names(rankings1) = entrezID

rankings1 = rankings1[!is.na(names(rankings1))]

rankings1 = sort(rankings1, decreasing = T)

edo4 = GSEA(rankings1, pvalueCutoff = 0.1, TERM2GENE = term2gene)
