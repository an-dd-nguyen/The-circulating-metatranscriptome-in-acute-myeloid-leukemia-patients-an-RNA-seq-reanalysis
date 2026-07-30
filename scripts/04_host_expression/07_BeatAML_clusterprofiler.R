library("org.Hs.eg.db", character.only = TRUE)
library(clusterProfiler)
library(enrichplot)
library(fgsea)
library(DOSE)
library(tidyverse)

#https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#pubmed-trend-of-enriched-terms

pathways = readRDS("../../Data/c_Bacteria_pathways.rds")

rankings = results_batch$logFC
names(rankings) = df$display_label

rankings1 = rankings

term2gene = read.delim("../../Data/Term2Gene_Entrez.txt",
                       header = F, sep = " ")
term2gene = term2gene[, -1]
term2gene = term2gene[-1, ]
colnames(term2gene) = c("TermID", "GeneID")
  

entrezID = mapIds(org.Hs.eg.db, names(rankings), 
                  column = "ENTREZID", 
                  keytype = "SYMBOL")

names(rankings) = entrezID

rankings = rankings[!is.na(names(rankings))]

rankings = sort(rankings, decreasing = T)

#rankings = rankings[abs(rankings) >=1]

edo = GSEA(rankings, pvalueCutoff = 0.1, TERM2GENE = term2gene)

rankings1 = rankings1[!is.na(names(rankings1))]

rankings1 = sort(rankings1, decreasing = T)



dotplot(edo)


edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')

#edo = enrichDGN(names(rankings[abs(rankings) >=1]))
#barplot(edo) 

cnetplot(edox, foldChange=rankings)

p1 = heatplot(edox, foldChange=rankings, showCategory=9)

# 1. Identify whether the column name is lowercase 'gene' or uppercase 'Gene' 
# (This varies depending on your version of the enrichplot package)
gene_column <- ifelse("gene" %in% colnames(p1$data), "gene", "Gene")

# 2. Extract the unique genes currently present across your 9 categories in the plot
genes_in_plot <- unique(p1$data[[gene_column]])

# 3. Filter rankings1 for those genes, sort by absolute logFC, and take the top 50
top_50_de_genes <- names(sort(abs(rankings1[names(rankings1) %in% genes_in_plot]), decreasing = TRUE))[1:50]

# Remove NA values in case the 9 pathways contain fewer than 50 unique genes total
top_50_de_genes <- na.omit(top_50_de_genes)

# 4. Filter the data frame inside the ggplot object to keep only those 50 genes
p1$data <- p1$data[p1$data[[gene_column]] %in% top_50_de_genes, ]

# 5. Plot the modified heatmap
p1


edot <- pairwise_termsim(edo)
p1 <- emapplot(edot)
p2 <- emapplot(edot, cex_category=1.5)
p3 <- emapplot(edot, layout="kk")
p4 <- emapplot(edot, cex_category=1.5,layout="kk") 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

ridgeplot(edo)

for(i in 1:nrow(edo@result))
{
  p1 <- gseaplot(edo, geneSetID = i, by = "runningScore", title = edo$Description[i])
  p2 <- gseaplot(edo, geneSetID = i, by = "preranked", title = edo$Description[i])
#p3 <- gseaplot(edo, geneSetID = 1, title = edo$Description[1])
  cowplot::plot_grid(p1, p2, ncol=1)
  title = paste0(edo@result$ID[i], ".pdf")
  ggsave(title, width = 7.2, height = 6, unit = "in", dpi = 600)
}


##Cibersort

cibersort_df = read_tsv("../Data/BeatAML_Cibersortx_Deconvolution.txt")
LSPC = apply(cibersort_df[, 5:7], 1, sum)
immune = cibersort_df$B + cibersort_df$CTL + cibersort_df$Monocyte + cibersort_df$NK + cibersort_df$Plasma + cibersort_df$`T` + cibersort_df$cDC
blast = cibersort_df$`GMP-like` + cibersort_df$`Mono-like` + cibersort_df$`ProMono-like` + cibersort_df$`cDC-like`

cibersort_df$LSPC = LSPC
cibersort_df$immune = immune
cibersort_df$blast = blast

cibersort_df_topbot = cibersort_df |> 
  filter(Mixture %in% c(top_ID_batch, bottom_ID_batch))

cibersort_df_topbot = cibersort_df_topbot[match(c(top_ID_batch, bottom_ID_batch),
                                                cibersort_df_topbot$Mixture), ]

test = c()
for(i in 19 : 21)
{
  top_vec = cibersort_df_topbot[1 : 24, i] |> pull()
  bottom_vec = cibersort_df_topbot[25:48, i] |> pull()
  test = c(test, wilcox.test(top_vec, bottom_vec)$p.value)
}

cibersort_dfx = cibersort_df |>
  filter(Mixture %in% genus_df$SampleID) |>
  arrange(desc(Mixture))
