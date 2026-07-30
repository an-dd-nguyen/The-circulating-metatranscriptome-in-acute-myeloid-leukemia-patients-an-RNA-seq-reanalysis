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

x = read_tsv("../Input_data/BeatAML_Normalized_top100_genus.txt") |>
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

#Genus level
df_samplegroup = read_tsv("../Input_data/BeatAML_all_genera_report.txt")

#Phylum level
df_phylum = read_tsv("../Input_data/BeatAML_all_phylum_report.txt", col_names = F)

colnames(df_phylum)[2] = "Reads"
colnames(df_phylum)[8] = "Levels"
colnames(df_phylum)[10] = "Name"
colnames(df_phylum)[11] = "Filename"

df_phylum = df_phylum |>
  select(c(Reads, Levels, Name, Filename)) |>
  filter(Levels == "P")

df_phylum_normalized = df_phylum %>%
  # 1. Create the phylum vector with the 5 specified levels
  mutate(
    Phylum = fct_other(
      Name, 
      keep = c("Proteobacteria", "Actinobacteria", "Firmicutes", "Bacteroidetes"),
      other_level = "Others"
    )
  ) %>%
  
  # 2. Group by Filename and the new phylum levels to aggregate reads
  group_by(Filename, Phylum) %>%
  summarise(Reads = sum(Reads), .groups = "drop_last") %>% 
  
  # 3. Calculate the percentage normalization within each file
  mutate(normalize = (Reads / sum(Reads)) * 100) %>%
  ungroup() %>%
  mutate(Filename = str_remove(Filename, "_combined_non_human_non_contaminant_krakenuniq.report"))



df_cuti_reads = df_samplegroup |>
  filter(Name == "Cutibacterium") |>
  filter(SampleType != "Blood Derived Normal" & SampleType != "Control Analyte") |>
  select(c(Filename, Reads, Name)) |>
  arrange(Filename)

df_phylum_cuti = df_phylum_normalized |>
  filter(Phylum == "Actinobacteria") |>
  arrange(Filename) |>
  filter(Filename %in% df_cuti_reads$Filename)

df_phylum_cuti$Reads = df_phylum_cuti$Reads - df_cuti_reads$Reads

df_staphy_reads = df_samplegroup |>
  filter(Name == "Staphylococcus") |>
  filter(SampleType != "Blood Derived Normal" & SampleType != "Control Analyte") |>
  select(c(Filename, Reads, Name)) |>
  arrange(Filename)

df_phylum_staphy = df_phylum_normalized |>
  filter(Phylum == "Firmicutes") |>
  arrange(Filename) |>
  filter(Filename %in% df_staphy_reads$Filename)

df_phylum_staphy$Reads = df_phylum_staphy$Reads - df_staphy_reads$Reads  

cuti_lookup <- df_samplegroup %>%
  filter(Name == "Cutibacterium") %>%
  filter(SampleType != "Blood Derived Normal" & SampleType != "Control Analyte") %>%
  select(Filename, Cuti_Reads = Reads)

staphy_lookup <- df_samplegroup %>%
  filter(Name == "Staphylococcus") %>%
  filter(SampleType != "Blood Derived Normal" & SampleType != "Control Analyte") %>%
  select(Filename, Staphy_Reads = Reads)

# 2. Join them to the main dataframe and update the values conditionally
df_phylum_normalized <- df_phylum_normalized %>%
  # Bring in the genus counts safely mapped by Filename
  left_join(cuti_lookup, by = "Filename") %>%
  left_join(staphy_lookup, by = "Filename") %>%
  
  # Replace any missing (NA) values with 0 so subtraction doesn't result in NA
  mutate(
    Cuti_Reads   = coalesce(Cuti_Reads, 0),
    Staphy_Reads = coalesce(Staphy_Reads, 0)
  ) %>%
  
  # Perform the conditional subtraction based on the Phylum row
  mutate(
    Reads = case_when(
      Phylum == "Actinobacteria" ~ Reads - Cuti_Reads,
      Phylum == "Firmicutes"     ~ Reads - Staphy_Reads,
      TRUE                       ~ Reads # Keep everything else exactly as is
    )
  ) %>%
  
  # 3. Recalculate your normalization percentage since raw reads changed!
  group_by(Filename) %>%
  mutate(normalize = (Reads / sum(Reads)) * 100) %>%
  ungroup() %>%
  
  # Clean up temporary lookup columns
  select(-Cuti_Reads, -Staphy_Reads) %>%
  # Group by each file to calculate the new total reads within that file
  group_by(Filename) %>%
  
  # Recalibrate the percentage based on the updated read counts
  mutate(normalize = (Reads / sum(Reads)) * 100) %>%
  
  # Always ungroup after you are done with grouped calculations
  ungroup()



proteo_sorted_file = df_phylum_normalized |>
  filter(Phylum == "Proteobacteria") |>
  mutate(Percent_Pro = normalize) |>
  arrange((Percent_Pro)) |>
  mutate(Order_Pro = c(1:510))

join_phylum_df = left_join(proteo_sorted_file |> select(Filename, Order_Pro), df_phylum_normalized, by = "Filename")

join_phylum_df$Phylum = factor(join_phylum_df$Phylum, levels = c("Proteobacteria", "Actinobacteria", "Bacteroidetes",
                                                                 "Firmicutes", "Others"))


join_phylum_df = join_phylum_df |>
  filter(Filename %in% df_cuti_reads$Filename)


ggplot(join_phylum_df, aes(x = factor(Order_Pro), y = normalize)) +
  geom_col(aes(fill = Phylum), position = "stack", width = 1) +
  theme_bw(base_size = 12) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(mult = 0)) +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  labs(
    y = "Relative Abundance",
    x = "",
    fill = "Phylum")
