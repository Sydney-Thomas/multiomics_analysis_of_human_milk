library(tidyverse)
library(stringr)
library(tibble)
library(dplyr)
library(tidyr)
library(corrplot)
library(RColorBrewer)
library(Hmisc)
library(rstatix)
library(ggpubr)


HMO <- read.csv("./Data/HMO_scaled_centered.csv")
questionnaire <- read.csv("./Data/US/DFU_QuestionnaireData_03212022.csv")
metadata <- read.csv("./Data/US/Metadata_SEED.csv")

names <- c("Opioids", "Secretor", "Child_1_Delivery", "Child_1_Preterm_Birth", "Antibiotic", "SSRI", "Allergy", "Pain_reliver", "Caffeine", "Alcohol",
           "Breastfeeding_Experience", "Marijuana", "Depression", "Anxiety", "Asthma", "Diabetes", "Migraines", "Anemia", "HPV", "Breastfeeding_frequency_factor", "EBF")
metadata <- metadata %>% mutate(across(any_of(names), as.factor))
rm(names)

data <- questionnaire %>% 
  gather(Sample, Result, 4:ncol(questionnaire))
data$Sample <- str_replace(data$Sample, "X", "")
data <- data %>% 
  separate(Sample, c("Sample", "Child"), "_")
data$Test <- str_replace_all(data$Test, " ", "_")
data$Test <- str_replace_all(data$Test, "-", "_")
data$Test <- str_replace_all(data$Test, "_.mts.", "")
data$Test <- str_replace_all(data$Test, ":", "")
data <- left_join(data, metadata, by = "Sample")
data$Result <- as.numeric(data$Result)
data <- data %>% filter(!is.na(Result)) %>% filter(!is.na(Amigos))
data$Amigos <- factor(data$Amigos, levels = c("Low", "Medium", "High"))

# Just look at ITSEA scores
# ITSEA tests were performed for the most children in our study compared to all other tests
data <- data %>% filter(str_detect(Category, "ITSEA")) %>% filter(!str_detect(Test, "_Test|_Age"))

data$Test <- str_replace_all(data$Test, "dysreg", "Dysregulation")
data$Test <- str_replace_all(data$Test, "compet", "Competence")
data$Test <- str_replace_all(data$Test, "extern", "Externalizing")
data$Test <- str_replace_all(data$Test, "intern", "Internalizing")

# Run anova on ITSEA scores and 6DDL
fit <- data %>% group_by(Test, Months) %>% tukey_hsd(Result ~ Amigos)

# Plot results
# Here, we filter the data just to look at externalizing subsection scores at 12 months (and one score which was clearly a typo)
sig_fit <- data %>% filter(Months == 12) %>% filter(str_detect(Test, "Extern")) %>% filter(Test != "Externalizing_agress_def"|Result<5) %>% filter(!str_detect(Test, "_t"))

comparisons <- list(c("Low", "High"), c("Medium", "High"), c("Low", "Medium"))
ggplot(sig_fit, aes(x = Amigos, y = Result, color = Amigos)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  stat_compare_means(comparisons = comparisons) +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~Test, ncol = 3)

# Run correlations on ITSEA scores and 6DDL
ITSEA_6DDL <- data %>% group_by(Test, Months) %>% rstatix::cor_test(vars=c("Result"), vars2=c("z6DDL"), method="spearman")

#### Now let's look at Rothia/Veillonella and ITSEA scores #####################################
# Read in microbiome data
microbiome <- read.csv("./Data/US/Microbiome_table_SEED_genus_rclr.csv")

# Filter out Rothia and Veillonella (and also create a single column with the abundance of both)
microbiome <- microbiome %>% dplyr::select(Sample, contains("othia"), contains("eillon")) %>% rowwise() %>% mutate(RV = sum(Rothia + Veillonella))
microbiome[microbiome == 0] <- NA

# Merge datasets
microbiome <- inner_join(data, microbiome, by = "Sample", multiple = "all")
microbiome <- microbiome %>% dplyr::select(Sample, Test, Months, Result, Rothia, Veillonella, RV, Amigos) %>% unique()
microbiome <- microbiome %>% pivot_longer(Rothia:RV, names_to = "Microbe", values_to = "Abundance")

# Run correlation
ITSEA_microbes <- microbiome %>% group_by(Test, Months, Microbe) %>% rstatix::cor_test(vars=c("Result"), vars2=c("Abundance"), method="spearman")

# Plot results
# Here, we filter the data just to look at externalizing subsection scores at 12 months (and one score which was clearly a typo)
sig_fit <- microbiome %>% filter(Months == 12) %>% filter(str_detect(Test, "Extern")) %>% filter(Test != "Externalizing_agress_def"|Result<5) %>% filter(Microbe == "RV")
ggscatter(sig_fit, x = "Result", y = "Abundance", add = "reg.line",  conf.int = TRUE, 
          cor.method = "spearman") + stat_cor() + facet_wrap(~Test, scales = "free_x", ncol = 3)
