library(tidyverse)
library(stringr)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(car)
library(gridExtra)
library(vegan)
library(Rcpm)
library(corrplot)
library(RColorBrewer)
library(Hmisc)
library(ggExtra)
library(scales)
library(lemon)
library(igraph)
library(ggraph)
library(rstatix)
library(ggpubr)
library(factoextra)
library(FactoMineR)

########### Make CANOPUS barplots ###############################################################
### Read in Canopus data
coeff_ES <- read.csv("./Data/ES/canopus_HMO_correlations.csv")
coeff_US <- read.csv("./Data/US/canopus_HMO_correlations.csv")
coeff_ES$Dataset <- "Collado"
coeff_US$Dataset <- "SEED"
coeff_ES <- coeff_ES %>% filter(Metabolite %in% coeff_US$Metabolite)
coeff <- rbind(coeff_ES, coeff_US)
rm(coeff_ES, coeff_US)

sig <- coeff %>% filter(p<0.05)
coeff <- coeff %>% filter(Metabolite %in% sig$Metabolite)

canopus_all <- read.csv("./Data/US/canopus_SEED.csv")

# Correlation plots
data <- canopus_all %>%
  dplyr::select(canopus, Type)
data$Type <- factor(data$Type, levels = c("NPC.pathway", "NPC.superclass", "NPC.class", "ClassyFire.superclass", "ClassyFire.class", "ClassyFire.subclass", "ClassyFire.level.5"))
names(data)[1] <- "Metabolite"
coeff <- left_join(coeff, data, by = "Metabolite")
coeff$Metabolite <- str_replace_all(coeff$Metabolite, "_", " ")
coeff$Metabolite <- str_replace_all(coeff$Metabolite, " Classy", "")
coeff$Metabolite <- str_replace_all(coeff$Metabolite, " NPC", "")
coeff$sig <- ifelse(coeff$padj < 0.05, 1, ifelse(coeff$p<0.05, 2, 3))
coeff$z6DDL <- coeff$z6DDL*-1

# Bar plot ###############################
ggplot(coeff, aes(x=z6DDL, y=reorder(Metabolite, z6DDL))) +
  geom_bar(stat="identity", aes(fill=z6DDL)) +
  scale_fill_distiller(palette="RdBu") +
  facet_rep_grid(Type~Dataset, scales = "free", space = "free", repeat.tick.labels = T)

rm(sig, coeff, data, canopus_all)

########### Look at similar metabolites ##########################################################################################
tukey_US <- read.csv("./Data/US/Metabolite_Changes.csv")
tukey_ES <- read.csv("./Data/ES/Metabolite_Changes.csv")
corr_US <- read.csv("./Data/US/HMO_Metabolite_correlations.csv")
corr_ES <- read.csv("./Data/ES/HMO_Metabolite_correlations.csv")

tukey_US <- tukey_US %>% filter(Secretor == 1)
tukey_ES <- tukey_ES %>% filter(Secretor == 1)
corr_US <- corr_US %>% dplyr::select(Metabolite:padj)
corr_ES <- corr_ES %>% dplyr::select(Metabolite:padj)
all_US <- right_join(corr_US, tukey_US, by = "Metabolite", multiple="all")
all_ES <- right_join(corr_ES, tukey_ES, by = "Metabolite", multiple="all")
all_US <- all_US %>% filter(p<0.05 & p.adj<0.05)
all_ES <- all_ES %>% filter(p<0.1 & p.adj<0.1)

all_US <- all_US %>% dplyr::select(-(null.value:conf.high)) %>% dplyr::select(-(Network_Link)) %>% dplyr::select(Metabolite, Network_Number, everything())
all_ES <- all_ES %>% dplyr::select(-(null.value:conf.high)) %>% dplyr::select(-(Network_Link)) %>% dplyr::select(Metabolite, Network_Number, everything())

# Choose a metabolite to look at
i <- "X1766" 

# Difference between HMO groups in US samples
norm_SEED <- read.csv("./Data/US/Metabolites_normalized.csv")
metadata_SEED <- read.csv("./Data/US/Metadata_SEED.csv")
norm_SEED <- left_join(norm_SEED, metadata_SEED, by = "Sample")
norm_SEED <- norm_SEED %>% filter(!is.na(Amigos)) %>% filter(!is.na(Secretor))
norm_SEED$Amigos <- factor(norm_SEED$Amigos, levels = c("Low", "Medium", "High"))
rm(metadata_SEED)

# Look at individual metadata 
ggplot(norm_SEED, aes(x = Amigos, y = (!!sym(i)), color = Amigos)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=0.9, shape=16) +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~Secretor)

# Choose a metabolite to look at
j <- "X13415" 

# Difference between HMO groups in ES samples
norm_Collado <- read.csv("./Data/ES/Metabolites_normalized.csv")
metadata_Collado <- read.csv("./Data/ES/Metadata_Collado.csv")
norm_Collado <- left_join(norm_Collado, metadata_Collado, by = "Sample")
norm_Collado <- norm_Collado %>% filter(!is.na(Amigos)) %>% filter(!is.na(Secretor))
norm_Collado$Amigos <- factor(norm_Collado$Amigos, levels = c("Low", "Medium", "High"))
rm(metadata_Collado)

# Look at individual metadata 
ggplot(norm_Collado, aes(x = Amigos, y = (!!sym(j)), color = Amigos)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=0.9, shape=16) +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~Secretor)

# Combine metabolites
norm_ES <- norm_Collado %>% dplyr::select(Sample, j, z6DDL, Secretor, Amigos)
colnames(norm_ES)[2] <- "Metabolite"
norm_ES$Dataset <- "ES"
norm_US <- norm_SEED %>% dplyr::select(Sample, i, z6DDL, Secretor, Amigos)
colnames(norm_US)[2] <- "Metabolite"
norm_US$Dataset <- "US"
norm <- rbind(norm_ES, norm_US)
rm(norm_ES, norm_US)

comparisons <- list(c("Low", "High"), c("Medium", "High"), c("Low", "Medium"))
ggplot(norm, aes(x = Amigos, y = Metabolite, color = Amigos)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  stat_compare_means(comparisons = comparisons) +
  #theme(legend.position = "none") +
  geom_jitter(aes(shape=Dataset)) +
  facet_wrap(~Secretor)

ggscatter(norm, x = "z6DDL", y = "Metabolite", add = "reg.line",  conf.int = TRUE, 
          cor.method = "spearman", shape = "Dataset") + stat_cor() + facet_wrap(~Secretor)

########### Metabolite correlation to bacterial genera ####################################################
########### Read in normalized metabolites ################################################################
norm_SEED <- read.csv("./Data/US/Metabolites_normalized.csv")
norm_Collado <- read.csv("./Data/ES/Metabolites_normalized.csv")
metadata_Collado <- read.csv("./Data/ES/Metadata_Collado.csv")
metadata_Collado <- metadata_Collado %>% dplyr::select(Sample, SampleID)
norm_Collado <- inner_join(metadata_Collado, norm_Collado, by = "Sample")
norm_Collado <- norm_Collado %>% dplyr::select(-Sample)
norm_Collado <- norm_Collado %>% filter(SampleID != "")
rm(metadata_Collado)

Microbiome_SEED <- read.csv("./Data/US/Microbiome_table_SEED_genus_rclr.csv")
Microbiome_Coll <- read.csv("./Data/ES/Microbiome_table_Collado_genus_rclr.csv")
Microbiome_SEED[Microbiome_SEED == 0] <- NA
Microbiome_Coll[Microbiome_Coll == 0] <- NA

norm_Collado <- inner_join(norm_Collado, Microbiome_Coll, by = "SampleID")
norm_SEED <- inner_join(norm_SEED, Microbiome_SEED, by = "Sample")
rm(Microbiome_Coll, Microbiome_SEED)

tukey_US <- read.csv("./Data/US/Metabolite_Changes.csv")
tukey_ES <- read.csv("./Data/ES/Metabolite_Changes.csv")
corr_US <- read.csv("./Data/US/HMO_Metabolite_correlations.csv")
corr_ES <- read.csv("./Data/ES/HMO_Metabolite_correlations.csv")

tukey_US <- tukey_US %>% filter(Secretor == 1)
tukey_ES <- tukey_ES %>% filter(Secretor == 1)
corr_US <- corr_US %>% dplyr::select(Metabolite:padj)
corr_ES <- corr_ES %>% dplyr::select(Metabolite:padj)
all_US <- right_join(corr_US, tukey_US, by = "Metabolite", multiple="all")
all_ES <- right_join(corr_ES, tukey_ES, by = "Metabolite", multiple="all")
all_US <- all_US %>% filter(p<0.05 & p.adj<0.05)
all_ES <- all_ES %>% filter(p<0.1 & p.adj<0.1)
all_US$Metabolite <- paste0("X",all_US$Metabolite)
all_ES$Metabolite <- paste0("X",all_ES$Metabolite)
rm(tukey_US, tukey_ES, corr_US, corr_ES)

# Function - flatten correlation matrix (from sthda website)
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

# For US samples
data <- norm_SEED %>% column_to_rownames(var = "Sample")
data <- data %>% dplyr::select(where(~mean(is.na(.)) < 0.7)) %>% as.matrix()
corr <- rcorr(data, type = "spearman")
corr_flat <- flattenCorrMatrix(corr$r,corr$P)
corr_flat <- corr_flat %>% filter(str_detect(row, "^X1|^X2|^X3|^X4|^X5|^X6|^X7|^X8|^X9")) %>% filter(!str_detect(column, "^X1|^X2|^X3|^X4|^X5|^X6|^X7|^X8|^X9"))

# Check for metabolites that are significantly associated with 6DDL
corr_flat <- corr_flat %>% filter(row %in% all_US$Metabolite)
corr_flat$padj <- p.adjust(corr_flat$p, method = "BH")
corr_flat <- corr_flat %>% filter(p<0.05)

Mic_6DDL <- read.csv("C:/Users/sthomas/OneDrive - University of California, San Diego Health/SEED Grant/Microbiome/Microbiome_Comparison/Microbiome_HMO_correlation_network.csv")
Mic_6DDL <- Mic_6DDL %>% filter(from == "US") %>% dplyr::select(-from)
Mic_6DDL <- Mic_6DDL %>% dplyr::select(to, correlation) %>% dplyr::rename(Microbe_HMO_cor = correlation, Genera = to)

corr_flat <- corr_flat %>% rename(Metabolite = row, Genera = column, Microbe_Metabolite_cor = cor)
corr_flat <- inner_join(corr_flat, Mic_6DDL, by = "Genera")
corr_flat <- corr_flat %>% unique()
corr_flat$Direction <- ifelse(corr_flat$Microbe_Metabolite_cor<0, -1, 1)

write_csv(corr_flat, "Metabolite_Microbe_Correlations_SEED.csv")

# For ES samples
data <- norm_Collado %>% column_to_rownames(var = "SampleID")
data <- data %>% dplyr::select(where(~mean(is.na(.)) < 0.7)) %>% as.matrix()
corr <- rcorr(data, type = "spearman")
corr_flat <- flattenCorrMatrix(corr$r,corr$P)
corr_flat <- corr_flat %>% filter(str_detect(row, "^X1|^X2|^X3|^X4|^X5|^X6|^X7|^X8|^X9")) %>% filter(!str_detect(column, "^X1|^X2|^X3|^X4|^X5|^X6|^X7|^X8|^X9"))

# Check for metabolites that are significantly associated with 6DDL
corr_flat <- corr_flat %>% filter(row %in% all_ES$Metabolite)
corr_flat$padj <- p.adjust(corr_flat$p, method = "BH")
corr_flat <- corr_flat %>% filter(p<0.1)

Mic_6DDL <- read.csv("C:/Users/sthomas/OneDrive - University of California, San Diego Health/SEED Grant/Microbiome/Microbiome_Comparison/Microbiome_HMO_correlation_network.csv")
Mic_6DDL <- Mic_6DDL %>% filter(from == "ES") %>% dplyr::select(-from)
Mic_6DDL <- Mic_6DDL %>% dplyr::select(to, correlation) %>% dplyr::rename(Microbe_HMO_cor = correlation, Genera = to)

corr_flat <- corr_flat %>% rename(Metabolite = row, Genera = column, Microbe_Metabolite_cor = cor)
corr_flat <- inner_join(corr_flat, Mic_6DDL, by = "Genera")
corr_flat <- corr_flat %>% unique()
corr_flat$Direction <- ifelse(corr_flat$Microbe_Metabolite_cor<0, -1, 1)

write_csv(corr_flat, "Metabolite_Microbe_Correlations_Collado.csv")


#### PCA plots ###########################################
### Read in data and rclr transform
norm_SEED <- read.csv("./Data/US/Metabolites_normalized.csv")
norm_SEED[is.na(norm_SEED)] <- 0
norm_SEED <- column_to_rownames(norm_SEED, var = "Sample")
norm_SEED <- norm_SEED %>% decostand(method = "rclr")
norm_SEED <- norm_SEED %>% rownames_to_column(var = "Sample")
metadata_SEED <- read.csv("./Data/US/Metadata_SEED.csv")
metadata_SEED <- metadata_SEED %>% dplyr::select(Sample, Amigos)
norm_SEED <- inner_join(metadata_SEED, norm_SEED, by = "Sample")
norm_SEED <- norm_SEED %>% filter(!is.na(Amigos))

norm_Collado <- read.csv("./Data/ES/Metabolites_normalized.csv")
norm_Collado[is.na(norm_Collado)] <- 0
norm_Collado <- column_to_rownames(norm_Collado, var = "Sample")
norm_Collado <- norm_Collado %>% decostand(method = "rclr")
norm_Collado <- norm_Collado %>% rownames_to_column(var = "Sample")
metadata_Collado <- read.csv("./Data/ES/Metadata_Collado.csv")
metadata_Collado <- metadata_Collado %>% dplyr::select(Sample, Amigos)
norm_Collado <- inner_join(metadata_Collado, norm_Collado, by = "Sample")
norm_Collado <- norm_Collado %>% filter(!is.na(Amigos))
rm(metadata_Collado, metadata_SEED)

pca_Collado <- PCA(norm_Collado[,3:ncol(norm_Collado)], scale.unit = FALSE, ncp = 5, graph = FALSE)
pca_SEED <- PCA(norm_SEED[,3:ncol(norm_SEED)], scale.unit = FALSE, ncp = 5, graph = FALSE)
fviz_pca_ind(pca_Collado, axes = c(1,2), habillage = as.factor(norm_Collado$Amigos), label = "none", addEllipses = TRUE, ellipse.level = 0.95, invisible="quali")
fviz_pca_ind(pca_SEED, axes = c(1,2), habillage = as.factor(norm_SEED$Amigos), label = "none", addEllipses = TRUE, ellipse.level = 0.95, invisible="quali")
fviz_pca_biplot(pca_Collado, axes = c(1,2), habillage = as.factor(norm_Collado$Amigos), addEllipses = TRUE, select.var = list(contrib = 5), label = "var", invisible="quali", col.var = "black")
fviz_pca_biplot(pca_SEED, axes = c(1,2), habillage = as.factor(norm_SEED$Amigos), addEllipses = TRUE, select.var = list(contrib = 5), label = "var", invisible="quali", col.var = "black")

pca_matrix_Collado <- norm_Collado %>% dplyr::select(-Sample, -Amigos) %>% as.matrix()
dist_Collado <- vegdist(pca_matrix_Collado, method = "euclidean")
permanova_genus_Collado <- adonis2(dist_Collado ~ as.factor(norm_Collado$Amigos), norm_Collado, na.action = na.omit, method = "euclidean")

pca_matrix_SEED <- norm_SEED %>% dplyr::select(-Sample, -Amigos) %>% as.matrix()
dist_SEED <- vegdist(pca_matrix_SEED, method = "euclidean")
permanova_genus_SEED <- adonis2(dist_SEED ~ as.factor(norm_SEED$Amigos), norm_SEED, na.action = na.omit, method = "euclidean")
