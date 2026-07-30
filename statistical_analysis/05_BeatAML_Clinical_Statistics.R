library(readxl)
library(janitor)
library(performance)
library(caret)
library(randomForest)
library(ggpubr)
library(survival)
library(coin)
library(clinfun)
library(MASS)
library(ggtext)
library(tidyverse)

clinical = readRDS("../Downloads/BeatAML_clinical_clean.RDS")

drug_response = read_excel("../Downloads/BeatAML_DrugResponse.xlsx")
mapping = read_excel("../Downloads/beataml_waves1to4_sample_mapping.xlsx")

bug_norm = read_tsv("../Downloads/BeatAML_Results/Data/BeatAML_Normalized_top100_genus.txt")

drug_response = left_join(drug_response, mapping |> select(c(labId, dbgap_rnaseq_sample)), by = c("lab_id" = "labId"))

drug_response = drug_response |>
  filter(complete.cases(dbgap_rnaseq_sample))



ic50_counts <- drug_response |>
  group_by(inhibitor) |>
  summarise(
    ic50_greater_8_count = sum(ic50 > 8, na.rm = TRUE),
    ic50_lesser_2_count = sum(ic50 < 2, na.rm = TRUE),
    ratio = ic50_greater_8_count/ic50_lesser_2_count
  )

kept_drug = ic50_counts |>
  filter(ratio > 0.2 & ratio < 5) |>
  pull(inhibitor)


colnames(clinical)[81] = "Filename"




bug_norm = bug_norm |>
  dplyr::filter(SampleType != "Blood Derived Normal" & SampleType != "Control Analyte") |>
  filter(complete.cases(SampleGroup))

clinical = clinical |>
  filter(Filename %in% bug_norm$Filename)

clinical = clinical |> arrange(Filename)
bug_norm = bug_norm |> arrange(Filename)

drug_response = drug_response |> 
  filter(dbgap_rnaseq_sample %in% bug_norm$SampleID)


drug_response_wide = drug_response |>
  filter((inhibitor %in% kept_drug)) |>
  dplyr::select(inhibitor, ic50, dbgap_rnaseq_sample) |>
  filter(ic50 < 2 | ic50 > 8) |>
  pivot_wider(names_from = inhibitor, values_from =  ic50) |>
  arrange(dbgap_rnaseq_sample)

drug_response_wide_logical = drug_response_wide |>
  mutate(
    # Use across() to apply the transformation to all columns
    # except 'Filename'.
    across(
      .cols = -dbgap_rnaseq_sample, # Apply to all columns except 'Filename'
      .fns = ~case_when(
        .x > 8 ~ TRUE,
        .x < 2 ~ FALSE,
        TRUE ~ NA_real_ # Handles values between 5 and 9, and existing NAs
      )
    )
  )


##Clinical

age = bug_norm$SampleGroup |> as.factor()

bug_norm = bug_norm |> dplyr::select(-c(SampleID, SampleType, Filename, SampleGroup, centerID))

clinical <- clinical |> 
  mutate(
    # Fix: Ensure the parenthesis closes AFTER the ifelse arguments
    specificdxatinclusion_clean = factor(ifelse(specificdxatinclusion %in% c(
      "AML with mutated NPM1", "AML with myelodysplasia-related changes",
      "Acute myeloid leukaemia, NOS", "Therapy-related myeloid neoplasms"), 
      specificdxatinclusion, NA_character_)),
    
    specificdxatacquisition_clean = factor(ifelse(specificdxatacquisition %in% c(
      "AML with mutated NPM1", "AML with myelodysplasia-related changes",
      "Acute myeloid leukaemia, NOS", "Therapy-related myeloid neoplasms"), 
      specificdxatacquisition, NA_character_)),
    
    specimentype_clean = case_when(
      specimentype == "Bone Marrow Aspirate" ~ TRUE,
      specimentype == "Peripheral Blood"    ~ FALSE,
      TRUE                                   ~ NA),
    
    # Numeric thresholds to Boolean
    cumulativetreatmenttypecount_clean = cumulativetreatmenttypecount < 2,
    cumulativetreatmentregimencount_clean = cumulativetreatmentregimencount < 3,
    cumulativetreatmentstagecount_clean = cumulativetreatmentstagecount < 3,
    
    # Clinical outcomes to Boolean
    responsetoinductiontx_clean = case_when(
      responsetoinductiontx == "Complete Response" ~ TRUE,
      responsetoinductiontx == "Refractory" ~ FALSE,
      TRUE ~ NA # Handles UNKNOWN or other
    ),
    
    # Your prioritized treatment logic (Consolidation = TRUE)
    cumulativetreatmentstages_clean = case_when(
      str_detect(cumulativetreatmentstages, "Consolidation") ~ TRUE,
      str_detect(cumulativetreatmentstages, "Induction")     ~ FALSE,
      TRUE ~ NA
    ),
    
    # Risk categories
    eln2017_clean = case_when(
      str_detect(eln2017, "Adverse")      ~ "Adverse",
      str_detect(eln2017, "Intermediate") ~ "Intermediate",
      eln2017 == "Favorable"              ~ "Favorable",
      TRUE                                ~ NA_character_
    ) |> factor(levels = c("Favorable", "Intermediate", "Adverse")),
    
    fabblast_morphology = as.factor(fabblastmorphology),
    
    # Mutation status (Negative = FALSE, Mutation = TRUE)
    tp53_clean = case_when(
      tp53 == "negative" ~ FALSE,
      is.na(tp53)        ~ NA,
      TRUE               ~ TRUE
    ),
    
    runx1_clean = case_when(
      runx1 == "negative" ~ FALSE,
      is.na(runx1)        ~ NA,
      TRUE                ~ TRUE
    ),
    specimengroups_clean = ifelse(specimengroups_clean == "Initial AML Diagnosis",
                                  TRUE, FALSE),
    mostrecenttreatmenttype_clean = ifelse(mostrecenttreatmenttype_clean == "Standard Chemotherapy",
                                            TRUE, FALSE)
  )

clinical$eln2017_clean = factor(clinical$eln2017_clean, 
                                levels = c("Favorable", "Intermediate", "Adverse"))

clinical_clean = clinical |> dplyr::select(-c(Filename, rnaseq, exomeseq, 
                                              analysisexomeseq, dbgap_rnaseq_sample, 
                                              dbgap_subject_id, GDC_download,
                                              specificdxatinclusion, 
                                              specificdxatacquisition,
                                              specimentype, cumulativetreatmenttypecount,
                                              responsetoinductiontx, cumulativetreatmentstages,
                                              eln2017, consensus_sex, inferred_sex, inferred_ethnicity,
                                              centerid, cebpa_biallelic,
                                              specificdxatacquisition_mdsmpn,
                                              nonaml_mdsmpn_specificdxatacquisition,
                                              priormalignancytype, 
                                              priormalignancyradiationtx, 
                                              priormdsmpnmorethantwomths,
                                              priormdsmpn, priormpn, priormpnmorethantwomths,
                                              dxatinclusion, dxatspecimenacquisition,
                                              ageatspecimenacquisition, ageatdiagnosis,
                                              timeofsamplecollectionrelativetoinclusion,
                                              specimengroups, totaldrug, analysisdrug,
                                              cumulativetreatmenttypes, cumulativetreatmentregimens,
                                              responsedurationtoinductiontx, mostrecenttreatmenttype,
                                              currentregimen, currentstage, mostrecenttreatmentduration,
                                              karyotype, othercytogenetics, surfaceantigensimmunohistochemicalstains,
                                              asxl1, tp53, runx1, fabblastmorphology, sampletype))

vital_status = ifelse(clinical_clean$vitalstatus %in% c("Dead", "Alive"), clinical_clean$vitalstatus, NA_character_)
vital_status[vital_status == "Dead"] = 1
vital_status[vital_status == "Alive"] = 0
vital_status = as.numeric(vital_status)
surv_obj = Surv(time = clinical_clean$overallsurvival, event = vital_status)


test_vec = c()
results_vec = c()
year = clinical$ageatdiagnosis
gender = clinical$consensus_sex



for(i in 1:ncol(clinical_clean)) {
  # Pull clinical variable (length N)
  clin_vec = clinical_clean[[i]] 
  
  for(j in c(1:26, 28: 30)) {
    # Pull bug variable (length N)
    raw_bug_vec = log(bug_norm[[j]] + 1)
    
    # 1. Calculate Outlier Thresholds
    bug_mean <- mean(raw_bug_vec, na.rm = TRUE)
    bug_sd   <- sd(raw_bug_vec, na.rm = TRUE)
    
    # 2. Identify outliers (Values > 2SD from mean)
    # We use TRUE/FALSE to find them
    outlier_idx <- (raw_bug_vec > (bug_mean + 2 * bug_sd)) | 
      (raw_bug_vec < (bug_mean - 2 * bug_sd))
    
    # 3. REPLACE outliers with NA (Vector length stays exactly N)
    bug_vec_masked <- raw_bug_vec
    bug_vec_masked[outlier_idx] <- NA
    
    # 4. Run function - All vectors (bug, clin, age) are now same length N
    # Internal to getPvecA, R will ignore rows where bug is NA.
    results = getPvecAGY(var = clin_vec, bug = bug_vec_masked, age = age, gender = gender, year = year)
    
    # Handle NULL results
    if(is.null(results) || length(results) == 0) { results <- NA }
    
    # 5. Store results
    test = paste0(colnames(clinical_clean)[i], ".", colnames(bug_norm)[j])
    test_vec = c(test_vec, test)
    results_vec = c(results_vec, results)
  }
}


t = tibble(test = test_vec, result = results_vec, FDR = p.adjust(abs(results_vec), "BH"))
t = t[!grepl("vitalstatus|overallsurvival|causeofdeath", t$test), ]
t$FDR = p.adjust(abs(t$result), "BH")

t_prevotella_priormds = tibble(bug = log(bug_norm$Prevotella + 1),
                               var = clinical$priormdsmorethantwomths)
t_prevotella_priormds = t_prevotella_priormds |>
  filter(bug < mean(bug) + 2 * sd(bug) & bug > mean(bug) - 2 * sd(bug))
##Get rid of status in the legend + remove *

#Prevotella-Priormds
ggplot(t_prevotella_priormds, aes(x = var, y = (bug), fill = var)) +
  geom_boxplot() +
  xlab("") +
  ylab("log(*Prevotella* abundance + 1)") +
  theme_bw() +
  theme(axis.title.y = element_markdown(),
        axis.text.x = element_blank()) +
  labs(fill = "Prior MDS") +
  geom_text(
    data = data.frame(x = 1.5, y = 5.7, label = "Adjusted P = 0.08"),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE, # Cleans out the 'fill = var' mapping heritage
    size = 4
  )

##Drug response


#drug_response_wide = drug_response_wide[, -c(31, 42)]

#One for presence of microbes, one for abundance using the same loop
#but different variable name

bug_norm_drug = bug_norm |>
  filter(SampleID %in% drug_response_wide$dbgap_rnaseq_sample) |>
  filter(complete.cases(SampleGroup)) |>
  arrange(SampleID)



drug_response_wide_logical = drug_response_wide_logical |>
  filter(dbgap_rnaseq_sample %in% bug_norm_drug$SampleID) |>
  arrange(dbgap_rnaseq_sample)



test_vec = c()
results_vec = c()
spearman_vec = c()
ratio_vec = c()

clinical_drug = clinical |>
  filter(Filename %in% bug_norm_drug$Filename) |>
  arrange(dbgap_rnaseq_sample)


age = as.factor(bug_norm_drug$SampleGroup)
gender = as.factor(clinical_drug$consensus_sex)
year = clinical_drug$ageatdiagnosis

drug_response_auc_wide = drug_response_auc_wide |>
  filter(Filename %in% bug_norm_drug$Filename) |>
  arrange(Filename)


#Logical drug
drug_response_wide_logical = drug_response_wide_logical |>
  filter(dbgap_rnaseq_sample %in% bug_norm_drug$SampleID) |>
  arrange(dbgap_rnaseq_sample)

test_vec = c()
results_vec = c()

for(i in 2:57)
{
  var = as.logical((pull(drug_response_wide_logical[, i])))
  for(j in c(6:31,33 : 35))
  {
    bug = (pull(bug_norm_drug[, j]))
    #t = tibble(drug = var, bug = bug, age = age) 
    results = getPvecAGY(var = var, bug = bug, age = age, year = year, gender = gender)
    test = paste0(colnames(drug_response_wide_logical)[i], ".", colnames(bug_norm[, j]))
    results_vec = c(results_vec, results)
    test_vec = c(test_vec, test)
  }
}

t_bug_abundance_dug_logical = tibble(test = test_vec,
                                     result = results_vec,
                                     fdr = p.adjust(abs(results_vec), "BH"))



t_Moraxella_BEZ235 = tibble(bug = bug_norm_drug$Moraxella,
                            var = drug_response_wide_logical$BEZ235,
                            SampleGroup = age) |>
  na.omit() |>
  mutate(var = factor(case_when(var == 1 ~ "Resistant",
                         var == 0 ~ "Sensitive")))

## Moraxella - BEZ235
ggplot(t_Moraxella_BEZ235, aes(x = var, y = log(bug + 1))) +
  geom_boxplot() +
  geom_violin(aes(fill = var)) +
  geom_jitter() +
  theme_bw(base_size = 12) +
  ylab("log(*Moraxella* + 1)") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown()) +
  annotate("text", x = 1, y = 6, label = "Adjust P = 0.002")
#drug_norm_auc50 = data.frame(test = test_vec, results = results_vec)
drug_norm_auc50 = data.frame(test = test_vec, results = results_vec, spearman = spearman_vec)
drug_norm_auc50 = drug_norm_auc50 |>
  mutate(fdr = p.adjust(abs(results), method = "fdr")) 


