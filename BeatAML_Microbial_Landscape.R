library(conflicted)
library(readxl)
library(ggfortify)
library(factoextra)
library(plotly)
library(vegan)
library(ggpubr)
library(purrr)
library(broom)
library(forcats)
library(tidyverse)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

x = read_tsv("../Downloads/BeatAML_Results/Data/BeatAML_Normalized_top100_genus.txt") |>
  filter(complete.cases(SampleGroup))
x$centerID = factor(x$centerID)

genus_bray_curtis_dist = vegdist(x |> select_if(is.numeric)) |> as.matrix()

colnames(genus_bray_curtis_dist) = x$Filename
rownames(genus_bray_curtis_dist) = x$Filename


genus_bray_curtis = cmdscale(genus_bray_curtis_dist)

genus_bray_curtis_tibble = tibble(PC1  = genus_bray_curtis[,1],
                                  PC2 = genus_bray_curtis[,2],
                                  SampleGroup = factor(as.numeric(x$SampleGroup)), 
                                  SampleType = x$SampleType)

genus_bray_curtis_tibble = genus_bray_curtis_tibble |>
  mutate(SampleType = case_when(SampleType == "Control Analyte" ~ "12 Technical Replicates",
                                SampleType == "Blood Derived Normal" ~ "Healthy",
                                TRUE ~ "AML"))

(ggplot(genus_bray_curtis_tibble, aes(x = PC1, y = PC2 ,
                                      color = SampleType)) +
    geom_point(alpha = 0.5))
