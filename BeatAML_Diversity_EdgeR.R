library(edgeR)
library(tidyverse)
library(vegan)

df = read_tsv("../Downloads/beataml_waves1to4_counts_dbgap.txt")
genus_df = read_tsv("../Downloads/BeatAML_Results/Data/BeatAML_Normalized_top100_genus.txt")



genus_df = genus_df |>
  filter(SampleType != "Blood Derived Normal" & SampleType != "Control Analyte") |>
  filter(!duplicated(SampleID)) |>
  na.omit(SampleID) |>
  na.omit(SampleGroup)

counts = df[, -c(1 : 4)]
counts = as.matrix(counts)
row.names(counts) = df$stable_id




genus_df$SampleGroup = as.factor(genus_df$SampleGroup)
genus_df$centerID = as.factor(genus_df$centerID)

diversity = diversity(genus_df[, 6 : 105])

genus_df$diversity = diversity

#Faecalibacterium

top_and_bottom_bacteria <- genus_df %>%
  filter(Faecalibacterium > 0)

top_and_bottom_bacteria = top_and_bottom_bacteria %>%
  group_by(SampleGroup) %>%
  
  # Get the Top 2 (highest) diversity values
  # The 'with_ties = FALSE' ensures we only get two rows unless there are ties
  slice_max(order_by = Faecalibacterium, n = 2, with_ties = FALSE) %>%
  
  # Combine the top 2 results with the bottom 2 results
  bind_rows(
    top_and_bottom_bacteria %>%
      group_by(SampleGroup) %>%
      # Get the Bottom 2 (lowest) diversity values
      slice_min(order_by = Faecalibacterium, n = 2, with_ties = FALSE)
  ) %>%
  
  # Remove duplicate rows (which occur when the group size is small, e.g., 3 or 4)
  distinct() %>%
  
  # Sort the final result for a clean display
  arrange(SampleGroup, desc(Faecalibacterium)) 

##Manual inspection
count(top_and_bottom_bacteria, SampleGroup)

top_and_bottom_bacteria = top_and_bottom_bacteria |>
  filter(!(SampleGroup %in% c(7, 8 ,10)))

top_ID_batch = top_and_bottom_bacteria |> 
  group_by(SampleGroup) %>%
  
  # Selects the row(s) with the maximum value in the diversity column
  slice_max(order_by = Faecalibacterium, 
            n = 2, 
            with_ties = FALSE) %>% # 'n = 1' ensures only the single max row is kept
  
  ungroup() |>
  pull(SampleID)

bottom_ID_batch = top_and_bottom_bacteria |> 
  group_by(SampleGroup) %>%
  
  # Selects the row(s) with the maximum value in the diversity column
  slice_min(order_by = Faecalibacterium, 
            n = 2, 
            with_ties = FALSE) %>% # 'n = 1' ensures only the single max row is kept
  
  ungroup() |>
  pull(SampleID)

group_batch = c(rep("Top", 16), rep("Bottom", 16))
group_batch = fct_relevel(group_batch, "Top")

batch = c(top_and_bottom_bacteria |> filter(SampleID %in% top_ID_batch) |> pull(SampleGroup),
          top_and_bottom_bacteria |> filter(SampleID %in% bottom_ID_batch) |> pull(SampleGroup))

batch = batch |> as.numeric() |> as.factor()

y_batch = DGEList(counts = counts[, c(top_ID_batch, bottom_ID_batch)],
                  group = group_batch,
                  genes = row.names(counts))

y_batch <- normLibSizes(y_batch)

design_batch <- model.matrix(~batch + group_batch)
y_batch <- estimateDisp(y_batch, design_batch, robust=TRUE)

fit_batch <- glmQLFit(y_batch, design_batch, robust=TRUE)

qlf_batch <- glmQLFTest(fit_batch)


results_batch_F = qlf_batch[["table"]]

#Enterococcus

top_and_bottom_bacteria <- genus_df %>%
  filter(Enterococcus > 0)

top_and_bottom_bacteria = top_and_bottom_bacteria %>%
  group_by(SampleGroup) %>%
  
  # Get the Top 2 (highest) diversity values
  # The 'with_ties = FALSE' ensures we only get two rows unless there are ties
  slice_max(order_by = Enterococcus, n = 2, with_ties = FALSE) %>%
  
  # Combine the top 2 results with the bottom 2 results
  bind_rows(
    top_and_bottom_bacteria %>%
      group_by(SampleGroup) %>%
      # Get the Bottom 2 (lowest) diversity values
      slice_min(order_by = Enterococcus, n = 2, with_ties = FALSE)
  ) %>%
  
  # Remove duplicate rows (which occur when the group size is small, e.g., 3 or 4)
  distinct() %>%
  
  # Sort the final result for a clean display
  arrange(SampleGroup, desc(Enterococcus)) 

## Manual inspection
count(top_and_bottom_bacteria, SampleGroup)

top_and_bottom_bacteria = top_and_bottom_bacteria |>
  filter(!(SampleGroup %in% c(10)))

top_ID_batch = top_and_bottom_bacteria |> 
  group_by(SampleGroup) %>%
  
  # Selects the row(s) with the maximum value in the diversity column
  slice_max(order_by = Enterococcus, 
            n = 2, 
            with_ties = FALSE) %>% # 'n = 1' ensures only the single max row is kept
  
  ungroup() |>
  pull(SampleID)

bottom_ID_batch = top_and_bottom_bacteria |> 
  group_by(SampleGroup) %>%
  
  # Selects the row(s) with the maximum value in the diversity column
  slice_min(order_by = Enterococcus, 
            n = 2, 
            with_ties = FALSE) %>% # 'n = 1' ensures only the single max row is kept
  
  ungroup() |>
  pull(SampleID)

group_batch = c(rep("Top", 20), rep("Bottom", 20))
group_batch = fct_relevel(group_batch, "Bottom")

batch = c(top_and_bottom_bacteria |> filter(SampleID %in% top_ID_batch) |> pull(SampleGroup),
          top_and_bottom_bacteria |> filter(SampleID %in% bottom_ID_batch) |> pull(SampleGroup))

batch = batch |> as.numeric() |> as.factor()

y_batch = DGEList(counts = counts[, c(top_ID_batch, bottom_ID_batch)],
                  group = group_batch,
                  genes = row.names(counts))

y_batch <- normLibSizes(y_batch)

design_batch <- model.matrix(~batch + group_batch)
y_batch <- estimateDisp(y_batch, design_batch, robust=TRUE)

fit_batch <- glmQLFit(y_batch, design_batch, robust=TRUE)

qlf_batch <- glmQLFTest(fit_batch)


results_batch_E = qlf_batch[["table"]]

clinical = readRDS("../Downloads/BeatAML_clinical_clean.RDS")

clinical = inner_join(genus_df |> select(c(Filename, SampleID)), 
                      clinical,
                      by = c("Filename" = "filename"))

top_clinical = clinical |>
  filter(SampleID %in% top_ID_batch)

bottom_clinical = clinical |>
  filter(SampleID %in% bottom_ID_batch)

cibersort = read_tsv("../Downloads/BeatAML_CIBERSORT_Results.txt")

cibersort = cibersort |> 
  filter(Mixture %in% c(top_ID_batch, bottom_ID_batch))

cibersort_top = cibersort |>
  filter(Mixture %in% top_ID_batch)

cibersort_bottom= cibersort |>
  filter(Mixture %in% bottom_ID_batch)

for(i in 2 : 29)
{
  print(wilcox.test(cibersort_top[, i] |> pull() |> median(), cibersort_bottom[, i] |> pull()))
}


