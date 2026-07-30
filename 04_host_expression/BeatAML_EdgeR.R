library(edgeR)
library(tidyverse)
library(vegan)

df = read_tsv("../Downloads/beataml_waves1to4_counts_dbgap.txt")
genus_df = read_tsv("../Downloads/BeatAML_Results/Data/BeatAML_Normalized_top100_genus.txt")


genus_df = genus_df |>
  filter(SampleType != "Blood Derived Normal" & SampleType != "Control Analyte") |>
  filter(!duplicated(SampleID)) |>
  na.omit(SampleID)

counts = df[, -c(1 : 4)]
counts = as.matrix(counts)
row.names(counts) = df$stable_id




genus_df$SampleGroup = as.factor(genus_df$SampleGroup)
genus_df$centerID = as.factor(genus_df$centerID)

genus_df$total_bacteria = apply(genus_df |> dplyr::select(is.numeric), 1, sum)

genus_df = genus_df |>
  arrange(desc(total_bacteria))


#genus_df = genus_df |>
#  filter(Filename %in%  colnames(genus_bray_curtis_dist))

#### Batch control
top_and_bottom_bacteria <- genus_df %>%
  group_by(SampleGroup) %>%
  
  # Get the Top 2 (highest) total_bacteria values
  # The 'with_ties = FALSE' ensures we only get two rows unless there are ties
  slice_max(order_by = total_bacteria, n = 2, with_ties = FALSE) %>%
  
  # Combine the top 2 results with the bottom 2 results
  bind_rows(
    genus_df %>%
      group_by(SampleGroup) %>%
      # Get the Bottom 2 (lowest) total_bacteria values
      slice_min(order_by = total_bacteria, n = 2, with_ties = FALSE)
  ) %>%
  
  # Remove duplicate rows (which occur when the group size is small, e.g., 3 or 4)
  distinct() %>%
  
  # Sort the final result for a clean display
  arrange(SampleGroup, desc(total_bacteria))

top_ID_batch = top_and_bottom_bacteria |> 
  group_by(SampleGroup) %>%
  
  # Selects the row(s) with the maximum value in the total_bacteria column
  slice_max(order_by = total_bacteria, 
            n = 2, 
            with_ties = FALSE) %>% # 'n = 1' ensures only the single max row is kept
  
  ungroup() |>
  pull(SampleID)

bottom_ID_batch = top_and_bottom_bacteria |> 
  group_by(SampleGroup) %>%
  
  # Selects the row(s) with the minimum value in the total_bacteria column
  slice_min(order_by = total_bacteria, 
            n = 2, 
            with_ties = FALSE) %>% # 'n = 1' ensures only the single max row is kept
  
  ungroup() |>
  pull(SampleID)



group_batch = c(rep("Top", 24), rep("Bottom", 24))
group_batch = fct_relevel(group_batch, "Bottom")

batch = c(top_and_bottom_bacteria |> filter(SampleID %in% top_ID_batch) |> pull(SampleGroup),
          top_and_bottom_bacteria |> filter(SampleID %in% bottom_ID_batch) |> pull(SampleGroup))

y_batch = DGEList(counts = counts[, c(top_ID_batch, bottom_ID_batch)],
            group = group_batch,
            genes = row.names(counts))

y_batch <- normLibSizes(y_batch)

design_batch <- model.matrix(~batch + group_batch)
y_batch <- estimateDisp(y_batch, design_batch, robust=TRUE)

fit_batch <- glmQLFit(y_batch, design_batch, robust=TRUE)

qlf_batch <- glmQLFTest(fit_batch)


results_batch = qlf_batch[["table"]]
#top_ID = top_ID[-c(1, 2)]

group = c(rep("Top", 10), rep("Bottom", 10))
group = fct_relevel(group, "Bottom")

batch = genus_df |> 
  filter(SampleID %in% c(top_ID, bottom_ID)) |>
  select(c(SampleID, SampleGroup)) |>
  pull(SampleGroup)

y = DGEList(counts = counts[, c(top_ID, bottom_ID)],
            group = group,
            genes = row.names(counts))

y <- normLibSizes(y)

design <- model.matrix(~group)
y <- estimateDisp(y, design, robust=TRUE)
plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)

qlf <- glmQLFTest(fit)


results = qlf[["table"]]


##Cibersort
##Cibersort

cibersort_df = read_tsv("../Downloads/CIBERSORTx_Job55_Results (1).txt")
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


top_and_bottom_bacteria = top_and_bottom_bacteria |>
  dplyr::select(c(SampleID)) |>
  mutate(Load_category = (c("High", "High", "Low", "Low"))) |>
  inner_join(cibersort_df, by = c("SampleID" = "Mixture"))

top_and_bottom_bacteria$LSPC = apply(top_and_bottom_bacteria[, 7:9], 1, sum)
top_and_bottom_bacteria$immune = top_and_bottom_bacteria$B + top_and_bottom_bacteria$CTL + top_and_bottom_bacteria$Monocyte + top_and_bottom_bacteria$NK + top_and_bottom_bacteria$Plasma + top_and_bottom_bacteria$`T` + top_and_bottom_bacteria$cDC
top_and_bottom_bacteria$blast = top_and_bottom_bacteria$`GMP-like` + top_and_bottom_bacteria$`Mono-like` + top_and_bottom_bacteria$`ProMono-like` + top_and_bottom_bacteria$`cDC-like`


top_and_bottom_bacteria_long = top_and_bottom_bacteria |>
  ungroup() |>
  dplyr::select(c(Load_category, LSPC, immune, blast)) |>
  pivot_longer(!Load_category, names_to = "CellPopulation", values_to = "Proportion")

library(ggpubr)
library(ggbeeswarm)
library(ggtext)

# 1. Define your colors (keeping them consistent)
blast_col  <- "#E69F00" 
immune_col <- "#0072B2" 
LSPC_col   <- "#CC79A7"

# 2. Re-create the labels and explicitly set the FACTOR LEVELS
# The order in the 'levels' vector determines the order on the plot
top_and_bottom_bacteria_long <- top_and_bottom_bacteria_long %>%
  mutate(CP_Label = case_when(
    CellPopulation == "blast"  ~ paste0("<span style='color:", blast_col,  "'>blast</span>"),
    CellPopulation == "immune" ~ paste0("<span style='color:", immune_col, "'>immune</span>"),
    CellPopulation == "LSPC"   ~ paste0("<span style='color:", LSPC_col,   "'>LSPC</span>")
  )) %>%
  mutate(CP_Label = factor(CP_Label, levels = c(
    paste0("<span style='color:", blast_col,  "'>blast</span>"),
    paste0("<span style='color:", immune_col, "'>immune</span>"),
    paste0("<span style='color:", LSPC_col,   "'>LSPC</span>")
  )))


# 3. Plotting
# Now, ggplot will automatically put blast first, then immune, then LSPC
ggplot(top_and_bottom_bacteria_long, aes(x = CP_Label, y = Proportion)) +
  geom_boxplot(aes(fill = Load_category), outlier.shape = NA, alpha = 0.4, position = position_dodge(width = 0.8)) +
  geom_quasirandom(aes(group = Load_category), dodge.width = 0.8, size = 1.8, alpha = 0.8) +
  scale_fill_manual(values = c("Low" = "#00BFC4", "High" = "#F8766D")) +
  scale_color_manual(values = c("blast" = blast_col, "immune" = immune_col, "LSPC" = LSPC_col)) +
  annotate(
    "text", 
    x = c(1, 2, 3),                        # 1 is your 1st category, 2 is your 2nd, etc.
    y = c(1.02, 1.02, 1.02),               # Height above each boxplot
    label = c("FDR = 0.34", "FDR = 0.06", "FDR = 0.34"),
    size = 4) +
  theme_bw(base_size = 12) +
  labs(y = "Proportion", x = "Cell Population", fill = "Microbial Load", color = NULL) +
  theme(axis.text.x = element_markdown(face = "bold", size = 14), legend.position = "right", legend.text = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 4)))

top_and_bottom_bacteria_long_individual_population = top_and_bottom_bacteria |>
  dplyr::select(!c(SampleGroup, SampleID, LSPC, blast, immune, `P-value`, Correlation, RMSE)) |>
  ungroup() |>
  dplyr::select(-SampleGroup) |>
  pivot_longer(!Load_category, names_to = "CellPopulation", values_to = "Proportion")


top_and_bottom_bacteria_long_individual_population$Cell_category = rep(NA, nrow(top_and_bottom_bacteria_long_individual_population))
top_and_bottom_bacteria_long_individual_population$Cell_category[grep("LSPC", top_and_bottom_bacteria_long_individual_population$CellPopulation, ignore.case = T)] = "LSPC"
top_and_bottom_bacteria_long_individual_population$Cell_category[grep("-like", top_and_bottom_bacteria_long_individual_population$CellPopulation, ignore.case = T)] = "blast"
top_and_bottom_bacteria_long_individual_population$Cell_category[grep("LSPC|-like", top_and_bottom_bacteria_long_individual_population$CellPopulation, ignore.case = T, invert = T)] = "immune"

top_and_bottom_bacteria_long_individual_population <- top_and_bottom_bacteria_long_individual_population %>%
  mutate(CP_Label = case_when(
    Cell_category == "blast"  ~ paste0("<span style='color:", blast_col,  "'>blast</span>"),
    Cell_category == "immune" ~ paste0("<span style='color:", immune_col, "'>immune</span>"),
    Cell_category == "LSPC"   ~ paste0("<span style='color:", LSPC_col,   "'>LSPC</span>")
  )) %>%
  mutate(CP_Label = factor(CP_Label, levels = c(
    paste0("<span style='color:", blast_col,  "'>blast</span>"),
    paste0("<span style='color:", immune_col, "'>immune</span>"),
    paste0("<span style='color:", LSPC_col,   "'>LSPC</span>")
  )))

# 1. Define your category colors
blast_col  <- "#E69F00" 
immune_col <- "#0072B2" 
LSPC_col   <- "#CC79A7"

# 2. Create the label mapping
# This takes your long dataframe and creates a unique list of populations and their colors
label_map <- top_and_bottom_bacteria_long_individual_population %>%
  distinct(CellPopulation, Cell_category) %>%
  mutate(color = case_when(
    Cell_category == "blast"  ~ blast_col,
    Cell_category == "immune" ~ immune_col,
    Cell_category == "LSPC"   ~ LSPC_col
  )) %>%
  # Create the HTML string
  mutate(html_label = paste0("<span style='color:", color, "'>", CellPopulation, "</span>"))

# Create a named vector for scale_x_discrete: c("B" = "<span...>B</span>", ...)
formatted_labels <- setNames(label_map$html_label, label_map$CellPopulation)

# 3. Plot
ggplot(top_and_bottom_bacteria_long_individual_population, 
       aes(x = CellPopulation, y = Proportion)) +
  
  # Boxplots
  geom_boxplot(aes(fill = Load_category), 
               outlier.shape = NA, 
               alpha = 0.4, 
               position = position_dodge(width = 0.8)) +
  
  # Quasirandom points - color matched to Load_category or a neutral grey
  geom_quasirandom(aes(group = Load_category), 
                   color = "grey30", 
                   dodge.width = 0.8, 
                   size = 0.8) +
  
  # Inject the colored HTML labels
  scale_x_discrete(labels = formatted_labels) +
  
  # Colors
  scale_fill_manual(values = c("Low" = "#00BFC4", "High" = "#F8766D")) +
  
  # Center p-values over the dodged boxes
  stat_compare_means(aes(group = Load_category), 
                     method = "wilcox.test", 
                     label = "p.signif", 
                     label.y = 1.05, 
                     size = 3.5) +
  
  # Labels and Theme
  theme_bw(base_size = 14) +
  labs(y = "Proportion", 
       x = "Cell Population",
       fill = "Microbial Load") +
  
  theme(
    # Angle labels and adjust horizontal justification so they align under the ticks
    axis.text.x = element_markdown(face = "bold", angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_text(margin = margin(t = 15)),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )
