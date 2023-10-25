library(tidyverse)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(car)
library(gridExtra)
library(corrplot)
library(RColorBrewer)
library(ggbreak)
library(Hmisc)
library(multcomp)
library(dplyr)
library(nnet)
library(mixOmics)
library(factoextra)
library(NbClust)
library(ggpubr)
library(FactoMineR)
library(vegan)
library(cluster)
library(ggfortify)
library(dendextend)
library(rstatix)


########### Read, normalize, and clean datasets ######################################################################################
HMO_US <- read.csv("./Data/US/HMO_SEED.csv")
HMO_US <- HMO_US %>% dplyr::select(-Secretor)
HMO_ES <- read.csv("./Data/ES/HMO_Collado.csv")
HMO_ES <- HMO_ES %>% dplyr::select(-Secretor)
HMO_AU <- read.csv("./Data/AU/HMO_AU.csv")
HMO_AU <- HMO_AU %>% dplyr::select(-Secretor, -Infant_age_days)
HMO_US$Dataset <- "US"
HMO_ES$Dataset <- "ES"
HMO_AU$Dataset <- "AU"

metadata_US <- read.csv("./Data/US/Metadata_US.csv")
metadata_ES <- read.csv("./Data/ES/Metadata_Collado.csv")
metadata_AU <- read.csv("./Data/AU/Metadata_Geddes.csv")

rnm <- c(Infant_age_days = "Child_1_Age", Infant_sex = "Child_1_Gender", Infant_sex = "infant_sex", C_section = "Child_1_Delivery", Preterm = "Child_1_Preterm_Birth", Parent_BMI = "Mother_BMI",
         Parent_Age = "Age", Parent_Age = "Mother_age", Ethnicity = "Race", Ethnicity = "ethnicity", Parity = "parity", Parent_BMI = "Pre_BMI", Breastfed_months = "Child_1_Previous_Breastfed_Length_months",
         Breastfed_months = "Breastfeeding_length_months")
metadata_US <- metadata_US %>% rename(any_of(rnm)) %>% dplyr::select(-z6DDL, -Amigos)
metadata_ES <- metadata_ES %>% rename(any_of(rnm)) %>% dplyr::select(-z6DDL, -Amigos)
metadata_AU <- metadata_AU %>% rename(any_of(rnm)) %>% dplyr::select(-z6DDL, -Amigos)

HMO_ES <- full_join(metadata_ES, HMO_ES, by = "Sample")
HMO_ES$Sample <- as.character(HMO_ES$Sample)
HMO_ES$Sample <- str_replace_all(HMO_ES$Sample, "w", "W")
HMO_US <- full_join(HMO_US, metadata_US, by = "Sample")
HMO_AU <- full_join(HMO_AU, metadata_AU, by = "Sample")

Microbiome_US <- read.csv("./Data/US/SEED_microbiome_shannon_genus.csv")
Microbiome_ES <- read.csv("./Data/ES/Collado_microbiome_shannon_genus.csv")
Microbiome_AU <- read.csv("./Data/AU/Geddes_microbiome_shannon_genus.csv")
Microbiome_US$Shannon[Microbiome_US$Shannon == 0] <- NA
Microbiome_ES$Shannon[Microbiome_ES$Shannon == 0] <- NA
Microbiome_AU$Shannon[Microbiome_AU$Shannon == 0] <- NA
Microbiome_US <- Microbiome_US %>% mutate(Shannon_eq = Shannon/(max(Shannon, na.rm=TRUE)))
Microbiome_ES <- Microbiome_ES %>% mutate(Shannon_eq = Shannon/(max(Shannon, na.rm=TRUE)))
Microbiome_AU <- Microbiome_AU %>% mutate(Shannon_eq = Shannon/(max(Shannon, na.rm=TRUE)))

HMO_ES <- full_join(Microbiome_ES, HMO_ES, by = "SampleID") %>% dplyr::select(-Sample) %>% rename(Sample = SampleID) %>% filter(Sample != "")
HMO_US <- full_join(Microbiome_US, HMO_US, by = "Sample") %>% filter(Sample != "")
HMO_AU <- full_join(Microbiome_AU, HMO_AU, by = "Study_ID") %>% dplyr::select(-Sample) %>% rename(Sample = Study_ID) %>% filter(Sample != "")

rm(metadata_US, metadata_ES, metadata_AU, Microbiome_US, Microbiome_ES, Microbiome_AU, rnm)

HMO_US <- HMO_US %>% 
  gather(HMO, value, 2:ncol(HMO_US)) %>%
  spread(Sample, value)
HMO_ES <- HMO_ES %>% 
  gather(HMO, value, 2:ncol(HMO_ES)) %>%
  spread(Sample, value)
HMO_AU <- HMO_AU %>% 
  gather(HMO, value, 2:ncol(HMO_AU)) %>%
  spread(Sample, value)

HMO <- full_join(HMO_US, HMO_AU, by = "HMO")
HMO <- full_join(HMO, HMO_ES, by = "HMO")

HMO <- HMO %>% 
  gather(Sample, value, 2:ncol(HMO)) %>%
  spread(HMO, value)

HMO <- HMO %>%
  filter(!is.na(Secretor)) %>% 
  filter(!is.na(Infant_age_days))


HMO$Secretor <- as.factor(HMO$Secretor)
HMO$Shannon <- as.numeric(HMO$Shannon)
HMO$Parent_Age <- as.numeric(HMO$Parent_Age)
HMO$Parent_BMI <- as.numeric(HMO$Parent_BMI)
HMO$milk_production <- as.numeric(HMO$milk_production)
HMO$Breastfeeding_frequency <- as.numeric(HMO$Breastfeeding_frequency)
HMO$Infant_age_days <- as.numeric(HMO$Infant_age_days)
HMO$Ethnicity <- str_replace_all(HMO$Ethnicity, "African", "Black")
HMO$Income <- str_replace_all(HMO$Income, "\\$", "")
HMO$Income <- str_replace_all(HMO$Income, "10,001 - 49,999", ">10,000")
HMO$Income <- str_replace_all(HMO$Income, "50,000 - 59,999", ">10,000")
HMO$Income <- str_replace_all(HMO$Income, "10,000 - 49,999", ">10,000")

rm(HMO_ES, HMO_AU, HMO_US)

HMO <- HMO %>% dplyr::select(Sample, X6_SL, LSTc, DSLNH, DFLNH, X2_FL, X3_FL, X3_SL, DFLac, DFLNT, DSLNT, FDSLNH, FLNH, LNFP_I, LNFP_II, LNFP_III, LNH, LNnT, LNT, LSTb, Shannon, Shannon_eq, SUM, everything()) %>% mutate(across(X6_SL:SUM, as.numeric))

# Normalize data by sample
HMO_s <- HMO %>% mutate(across(X6_SL:LSTb, ~ ./SUM))

# Normalized by dataset uncentered
HMO_u <- HMO %>% group_by(Dataset, Secretor) %>% mutate(across(X6_SL:LSTb, ~ c(scale(., center=FALSE))))

# Normalize by dataset centered
# This normalization does the best job of removing dataset-dependent differences and is used for most of the rest of the script
HMO <- HMO %>% group_by(Dataset, Secretor) %>% mutate(across(X6_SL:LSTb, ~ c(scale(., center=TRUE))))


### Run IPCA for each dataset ########################################################################
ipca_hmos <- list()  # new empty list
for (n in c("AU", "US", "ES", "US|AU|ES")) {
  data <- HMO_s %>% dplyr::select(Sample, Secretor, Dataset, X6_SL, LSTc, DSLNH, DFLNH, X2_FL, X3_FL, X3_SL, DFLac, DFLNT, DSLNT, FDSLNH, FLNH, LNFP_I, LNFP_II, LNFP_III, LNH, LNnT, LNT, LSTb)
  data <- data %>% column_to_rownames(var = "Sample") %>% filter(str_detect(Dataset, n)) %>% filter(Secretor == 0) %>% dplyr::select(-Dataset, -Secretor) %>% as.matrix()
  dist <- ipca(data, scale = TRUE)
  data1 <- data.frame(dist$loadings$X)
  data1$Dataset <- n
  data <- t(data)
  km <- kmeans(data, 5, iter.max = 200, nstart = 50)
  data1$Cluster <- km$cluster
  ipca_hmos[[n]] <- data1
}

ipca_hmos <- do.call(rbind, ipca_hmos)
ipca_hmos$Dataset <- str_replace_all(ipca_hmos$Dataset, c("\\|" = ""))
ipca_hmos$Dataset <- str_replace_all(ipca_hmos$Dataset, c("USAUES" = "Total"))
ipca_hmos <- rownames_to_column(ipca_hmos, var = "HMO")
ipca_hmos$HMO <- str_replace_all(ipca_hmos$HMO, c(AU.="", ES.="", US.="", "\\AU|US|ES." = ""))
ipca_hmos$Datset <- factor(ipca_hmos$Dataset, levels = c("AU", "ES", "US", "Total"))

ggplot(ipca_hmos, aes(x = IPC1, y = IPC2, label=HMO, color=as.factor(Cluster))) +
  geom_point() +
  geom_text(position = position_nudge(y=0.01, x=0.01)) +
  theme_light() +
  theme(legend.position = "none") +
  facet_wrap(~Dataset)

rm(n, dist, km, data, data1, ipca_hmos)

########## Run PCA #########################################################################################################
i <- "Dataset"

data <- HMO_s %>% dplyr::select(Sample, Secretor, (!!sym(i)), X6_SL, LSTc, DSLNH, DFLNH, X2_FL, X3_FL, X3_SL, DFLac, DFLNT, DSLNT, FDSLNH, FLNH, LNFP_I, LNFP_II, LNFP_III, LNH, LNnT, LNT, LSTb) %>% filter(!is.na((!!sym(i))))
data[i] <- as.factor(data[[i]])
data <- data %>% column_to_rownames(var = "Sample") %>% filter(Secretor == 0)
data1 <- data %>% dplyr::select(-Secretor, -(!!sym(i))) %>% as.matrix()

# PCA
pca <- PCA(data[,3:ncol(data)], scale.unit = FALSE, ncp = 5, graph = FALSE)
fviz_screeplot(pca)
pca_var <- get_pca_var(pca)
fviz_pca_ind(pca, axes = c(1,2), habillage = data[[i]], label = "none", addEllipses = TRUE, ellipse.level = 0.95, invisible="quali")
fviz_pca_var(pca, axes = c(1,2))
fviz_pca_biplot(pca, axes = c(1,2), habillage = data[[i]], addEllipses = TRUE, select.var = list(contrib = 5), label = "var", invisible="quali", col.var = "black")

dist_hmos <- vegdist(data1, method = "euclidean")
dist_hmos <- as.matrix(dist_hmos)
permanova <- adonis2(dist_hmos ~ as.factor(data[[i]]), data, na.action = na.omit, method = "euclidean")

rm(permanova, dist_hmos, pca, pca_var, data1, data)

########## 6DDL Calculations ######################################################################################
data <- HMO %>% ungroup() %>% dplyr::select(Sample, X6_SL, LSTc, DSLNH, DFLNH, everything()) 
data <- data %>% mutate(z6DDL = rowSums(data[,2:5], na.rm = TRUE)) %>% dplyr::select(Sample, z6DDL, everything())
data$z6DDL <- as.numeric(data$z6DDL)

# Make HMOs categorical variable based on child age
fit <- lm(z6DDL ~ log2(Infant_age_days) + Secretor, data=data)
pred <- as.data.frame(predict(fit, interval = "confidence", level = .999999999999))
data = cbind(data, pred)
data$Amigos <- ifelse(data$z6DDL<pred[,2], "Low", ifelse(data$z6DDL>pred[,3], "High", "Medium"))
rm(fit)
rm(pred)

data1 <- data %>% filter(Secretor == 1)

ggplot(data, aes(x = log2(Infant_age_days), y = z6DDL)) +
  geom_point(aes(color=Amigos, shape=Dataset)) +
  geom_ribbon(data=data1, aes(ymin = lwr, ymax = upr), alpha = 0.15) +
  theme(legend.position = "none") +
  geom_line(data=data1, aes(y=fit)) 

Amigos <- data %>%
  dplyr::select(Sample, Secretor, z6DDL, Amigos)
write_tsv(Amigos, "Amigos.tsv")

Amigos <- Amigos %>% dplyr::select(Sample, z6DDL, Amigos)
HMO_c <- full_join(HMO, Amigos, by = "Sample")
HMO_s1 <- full_join(HMO_s, Amigos, by = "Sample")
HMO_u1 <- full_join(HMO_u, Amigos, by = "Sample")
write_csv(HMO_c, "HMO_scaled_centered.csv")
write_csv(HMO_u1, "HMO_unscaled_centered.csv")
write_csv(HMO_s1, "HMO_sample_normalized.csv")
rm(HMO_c, HMO_s1, HMO_u1)

# If you want to check how many points are in each group
data %>% filter(str_detect(Dataset, "US")) %>% count(Amigos)

########## Correlation of all HMOs ########################################################################
# Use spearman correlation since HMOs aren't evenly distributed ###########################################
corr <- HMO_u %>% ungroup() %>% dplyr::select(Sample:LSTb)
corr <- left_join(corr, Amigos, by = "Sample")
corr <- corr %>%
  filter(Secretor == 1) %>%
  dplyr::select(-Secretor, -Amigos)
colnames(corr) <- gsub("X", "", colnames(corr))
corr <- column_to_rownames(corr, var = "Sample")
cor <- rcorr(as.matrix(corr), type=c("spearman"))

# Extract coefficients
coeff <- as.data.frame(cor$r)
coeff <- rownames_to_column(coeff, var="HMO")

# Extract p values
sig <- as.data.frame(cor$P)
sig <- lapply(sig[], \(x) p.adjust (x, method = "BH", n = length(x)))
sig <- as.data.frame(sig)
colnames(sig) <- gsub("X", "", colnames(sig))
HMO_labels <- coeff$HMO
sig <- cbind(HMO_labels, sig)
rm(HMO_labels)

coeff <- column_to_rownames(coeff, var="HMO")
sig <- column_to_rownames(sig, var = "HMO_labels")
coeff <- coeff*-1
corrplot(as.matrix(coeff),diag = FALSE, type = "lower", order = "alphabet", p.mat = as.matrix(sig),  sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9, insig = 'label_sig', pch.col = 'grey20', cl.pos = 'r', cl.ratio = 0.20, addgrid.col = NA)

rm(corr, cor, sig, coeff)

# Calculate correlation between HMOs and microbiome shannon diversity #####################################
# Use spearman correlation since HMOs aren't evenly distributed ###########################################
data <- Amigos %>% dplyr::select(-Secretor)
data <- full_join(HMO_s, data, by = "Sample")
corr <- data %>%
  filter(Secretor==1) %>%
  filter(Dataset == "AU") %>%
  dplyr::select(Sample:Shannon, z6DDL)
rownames(corr) <- NULL
corr <- column_to_rownames(corr, var = "Sample")
cor <- rcorr(as.matrix(corr), type=c("spearman"))

# Extract coefficients
coeff <- as.data.frame(cor$r)
coeff <- rownames_to_column(coeff, var="HMO")
coeff <- coeff %>%
  dplyr::select(HMO, Shannon) %>%
  filter(!str_detect(HMO, "Shannon"))

# Extract p values
sig <- as.data.frame(cor$P)
sig <- rownames_to_column(sig, var="HMO")
sig <- sig %>%
  dplyr::select(HMO, Shannon) %>%
  filter(!str_detect(HMO, "Shannon"))
sig <- lapply(sig[-1], \(x) p.adjust (x, method = "BH", n = length(x)))
sig <- as.data.frame(sig)
sig <- round(sig, digits = 5)
colnames(sig) <- paste("padj", colnames(sig), sep = "_")
coeff <- cbind(coeff, sig)

# Make plot
coeff$HMO <- str_replace(coeff$HMO, "X", "")
coeff <- coeff %>%
  arrange(HMO)
coeff1 <- coeff %>%
  dplyr::select(HMO, Shannon)
sig1 <- coeff %>%
  dplyr::select(HMO, padj_Shannon)
coeff1 <- column_to_rownames(coeff1, var="HMO")
sig1 <- column_to_rownames(sig1, var = "HMO")
corrplot(as.matrix(coeff1), is.corr = FALSE, cl.pos = 'r', cl.ratio = 1, addgrid.col = NA, p.mat = as.matrix(sig1),  sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 0.9, insig = 'label_sig', pch.col = 'grey20', col.lim = c(-0.4,0.4), cl.length = 9, tl.pos = "l")
rm(coeff, coeff1, corr, cor, sig, sig1)

######## OTHER POTENTIAL CONFOUNDERS #####################################################################################
data <- Amigos %>% dplyr::select(-Secretor)
data <- full_join(HMO, data, by = "Sample")
data$Amigos <- factor(data$Amigos, levels = c("Low", "Medium", "High"))
# Plot Parity vs Four Amigos ###################
data1 <- data %>% 
  filter(!is.na(Parity)) %>% 
  filter(Parity != "12") 
data1 <- data1 %>% 
  group_by(Amigos, Secretor) %>% 
  count(Parity)
ggplot(data1, aes(x= Amigos, y = n, fill = as.factor(Parity))) + 
  geom_bar(stat = "identity", position = "fill") + 
  theme_minimal() +
 # theme(legend.position = "none") +
  facet_wrap(~Secretor) 

fit <- data %>% group_by(Secretor) %>% tukey_hsd(z6DDL ~ Parity)

# Plot EBF vs. Four amigos ##################################################################################
data1 <- data %>% 
  filter(!is.na(EBF))
data1 <- data1 %>% 
  group_by(Amigos, Secretor) %>% 
  count(EBF)
ggplot(data1, aes(x= Amigos, y = n, fill = EBF)) + 
  geom_bar(stat = "identity", position = "fill") + 
  theme_minimal() +
  facet_wrap(~Secretor) +
  theme(legend.position = "none")

data2 <- data1 %>% ungroup() %>% pivot_wider(names_from = "EBF", values_from = "n") %>% filter(Secretor == 0) %>% dplyr::select(-Secretor)
data2 <- data2 %>% column_to_rownames(var="Amigos") %>% as.matrix(labels=TRUE)
fit <- pairwise_fisher_test(data2)
fit <- data %>% group_by(Secretor) %>% tukey_hsd(z6DDL ~ EBF)
rm(data2)

# Plot BMI vs Four Amigos ####################
data1 <- data %>% 
  filter(!is.na(Parent_BMI))

ggplot(data1, aes(x = Amigos, y = Parent_BMI, fill = Amigos)) +
  geom_violin() +
  theme_minimal() +
  geom_boxplot(width = 0.1) +
  theme(legend.position = "none") +
  facet_wrap(~Secretor)

fit <- data1 %>% group_by(Secretor) %>% tukey_hsd(Parent_BMI ~ Amigos)

# Plot Parent Age vs Four Amigos ####################
data1 <- data %>% 
  filter(!is.na(Parent_Age))

ggplot(data1, aes(x = Amigos, y = Parent_Age, fill = Amigos)) +
  geom_violin() +
  theme_minimal() +
  geom_boxplot(width = 0.1) +
  theme(legend.position = "none") +
  facet_wrap(~Secretor)

fit <- data1 %>% group_by(Secretor) %>% tukey_hsd(Parent_Age ~ Amigos)

# Ethnicity #############################################
data1 <- data %>% 
  filter(!is.na(Ethnicity))
data1 <- data1 %>% 
  group_by(Amigos, Secretor) %>% 
  count(Ethnicity)
ggplot(data1, aes(x= Amigos, y = n, fill = Ethnicity)) + 
  geom_bar(stat = "identity") +
  theme_minimal() +
  #theme(legend.position = "none") +
  facet_wrap(~Secretor) 

fit <- data %>% group_by(Secretor) %>% tukey_hsd(z6DDL ~ Parity)

# Infant sex #######################################################
data1 <- data %>% 
  filter(!is.na(Infant_sex))
data1 <- data1 %>% 
  group_by(Amigos, Secretor) %>% 
  count(Infant_sex)

ggplot(data1, aes(x= Amigos, y = n, fill = Infant_sex)) + 
  geom_bar(stat = "identity", position = "fill") + 
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~Secretor) 

data2 <- data1 %>% ungroup() %>% pivot_wider(names_from = "Infant_sex", values_from = "n") %>% filter(Secretor == 1) %>% dplyr::select(-Secretor)
data2 <- data2 %>% column_to_rownames(var="Amigos") %>% as.matrix(labels=TRUE)
fit <- pairwise_fisher_test(data2)
fit <- data %>% group_by(Secretor) %>% tukey_hsd(z6DDL ~ Infant_sex)

rm(data2)

# Breastfed Length ####################################
data1 <- data %>% 
  filter(!is.na(Breastfed_months))
data1$Breastfed_months <- as.numeric(data1$Breastfed_months)
data1 <- data1 %>% filter(Breastfed_months < 37)

ggplot(data1, aes(x = Amigos, y = Breastfed_months, fill = Amigos)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=0.9, shape=16) +
  theme_minimal() +
  facet_wrap(~Secretor) +
  theme(legend.position = "none") 

fit <- data1 %>% group_by(Secretor) %>% tukey_hsd(Breastfed_months ~ Amigos)
average <- data1 %>% group_by(Secretor, Amigos) %>% summarise(mn = median(Breastfed_months), std = sd(Breastfed_months))
rm(average)

# Milk Production ####################################
data1 <- data %>% 
  filter(!is.na(milk_production))

ggplot(data1, aes(x = Amigos, y = milk_production, fill = Amigos)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=0.9, shape=16) +
  theme_minimal() +
  facet_wrap(~Secretor) +
  theme(legend.position = "none")

fit <- data1 %>% group_by(Secretor) %>% tukey_hsd(milk_production ~ Amigos)

# Breastfeeding Frequency ####################################
data1 <- data %>% 
  filter(!is.na(Breastfeeding_frequency))

ggplot(data1, aes(x = Amigos, y = Breastfeeding_frequency, fill = Amigos)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme_minimal() +
  facet_wrap(~Secretor) +
  theme(legend.position = "none")

fit <- data1 %>% group_by(Secretor) %>% tukey_hsd(Breastfeeding_frequency ~ Amigos)

rm(p, fit, data1)

# Income ############################################
data1 <- data %>% 
  filter(!is.na(Income))

data1 <- data1 %>% 
  group_by(Amigos, Secretor) %>% 
  count(Income)
ggplot(data1, aes(x= Amigos, y = n, fill = as.factor(Income))) + 
  geom_bar(stat = "identity", position = "fill") + 
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~Secretor) 

######## Plot differences in all HMOs between groups #######################################################################
data <- Amigos %>% dplyr::select(-Secretor)
data <- full_join(HMO_s, data, by = "Sample")
data <- data %>% ungroup() %>%
  dplyr::select(Sample, Secretor, Amigos, X6_SL, LSTc, DSLNH, DFLNH, DFLac, DFLNT, DSLNT, FDSLNH, FLNH, LNFP_I, LNFP_II, LNFP_III, LNH, LNnT, LNT, LSTb, X2_FL, X3_FL, X3_SL)

data <- data %>% pivot_longer(X6_SL:X3_SL, names_to = "HMO", values_to = "value")
data <- data %>%
  group_by(Secretor, Amigos, HMO) %>%
  summarise_at(c("value"), mean)
# Change levels to match levels of fold change between low and high groups
data$HMO <- factor(data$HMO, levels = c("X6_SL", "DSLNH", "LSTc", "DFLNH", "FLNH", "LNFP_I", "LNT", "DSLNT", "LNH", "LNFP_III", "X3_FL", "FDSLNH", "LNFP_II", "LNnT", "X2_FL", "DFLNT", "DFLac", "LSTb", "X3_SL"))
data$Amigos <- factor(data$Amigos, levels = c("Low", "Medium", "High"))

colorcount <- length(unique(data$HMO))
getPalette <- colorRampPalette(brewer.pal(11, "RdBu"))
ggplot(data, aes(x = Amigos, y = value, fill = HMO)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = getPalette(colorcount)) +
  facet_wrap(~Secretor)

#### Do the stats #######################
data <- Amigos %>% dplyr::select(-Secretor)
data <- full_join(HMO_s, data, by = "Sample")
groups <- data %>%
  dplyr::select(Sample, Secretor, Amigos, Infant_age_days, X6_SL, LSTc, DSLNH, DFLNH, DFLac, DFLNT, DSLNT, FDSLNH, FLNH, LNFP_I, LNFP_II, LNFP_III, LNH, LNnT, LNT, LSTb, X2_FL, X3_FL, X3_SL) %>%
  filter(Secretor == 1) %>%
  dplyr::select(-Secretor)
groups <- groups %>% pivot_longer(X6_SL:X3_SL, names_to = "HMO", values_to = "value")

# Now you can start doing calculations. First, average all the replicates in each sample set
average <- groups %>%
  group_by(Amigos, HMO) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)
# You'll get a warning on this command, it just means that R couldn't average your sample names, which is probably good since they're letter and not numbers

# Next calculate standard deviations for each set of replicates
stdev <- groups %>%
  group_by(Amigos, HMO) %>%
  summarise_if(is.numeric, sd, na.rm = TRUE)

# Now calculate fold changes (TRUE/FALSE)
fold_change <- average %>%
  group_by(HMO) %>%
  mutate_if(is.numeric, ~ log2(./.[Amigos == 'Low']))
fold_change <- fold_change %>%
  filter(Amigos == "High") %>%
  dplyr::select(-Infant_age_days)


# Do glm for all HMOs
long <- groups

# Convert to factors
long <- split(long, paste(long$HMO))

tmpfn <- function(dat) {
  fit <- multinom(Amigos ~ value + Infant_age_days, data=dat)
  p <-  (1 - pnorm(abs(summary(fit)$coefficients/summary(fit)$standard.errors), 0, 1)) * 2
  as.data.frame(p)
}

pval <- long %>%
  purrr::map(tmpfn) %>%
  bind_rows(.id = "HMO")
pval <- rownames_to_column(pval, var = "Type")
pval$Type <- str_replace(pval$Type, "Medium.*", "Medium")
pval$Type <- str_replace(pval$Type, "High.*", "High")
pval <- pval %>%
  dplyr::select(HMO, Type, value)
pval <- pval %>%
  mutate_if(is.numeric, ~ p.adjust(.x, method = "fdr", n = length(.x)))

rm(long)
rm(tmpfn)

########### Calculate average values for each HMO using metadata categories ############################
########### Need to do this separately for secretor/non-secretor ####################################
data <- Amigos %>% dplyr::select(-Secretor)
data <- full_join(HMO, data, by = "Sample")
groups <- data %>% ungroup %>% 
  dplyr::select(Sample:LSTb, C_section, Migraines, Pain_reliever, SSRI, Alcohol, Anxiety, Depression, EBF, Opioids, Preterm, z6DDL, Secretor, Infant_age_days) %>% 
  filter(Secretor == 1) %>% 
  dplyr::select(-Secretor)

groups <- groups %>% 
  pivot_longer(C_section:Preterm, names_to="Condition", values_to = "Count")
groups <- groups %>% 
  filter(!is.na(Count))

# Now you can start doing calculations. First, average all the replicates in each sample set 
average <- groups %>% 
  group_by(Condition, Count) %>% 
  summarise_if(is.numeric, mean, na.rm = TRUE)

# Next calculate standard deviations for each set of replicates
stdev <- groups %>% 
  group_by(Condition, Count) %>% 
  summarise_if(is.numeric, sd, na.rm = TRUE)

# Now calculate fold changes (TRUE/FALSE)
fold_change <- average %>% 
  group_by(Condition) %>% 
  mutate_if(is.numeric, ~ .-.[Count == 'FALSE'])
fold_change <- fold_change %>% 
  filter(Count == "TRUE") %>% 
  dplyr::select(-Count)

# Do glm for all HMOs
long <- groups
long$Sample_Type <- paste(long$Condition,"-",long$Count,sep="")
long <- long %>% 
  dplyr::select(Sample_Type, Infant_age_days, where(is.numeric))
long <- long %>%   
  gather(HMO, value, 3:ncol(long)) 
long <-separate(long, Sample_Type, c("Factor1", "Factor2"), sep = "-")

# Convert to factors
long$Factor1 <- as.factor(long$Factor1)
long$Factor2 <- as.factor(long$Factor2)
long <- split(long, paste(long$HMO, long$Factor1))

tmpfn <- function(dat) {
  fit <- glm(Factor2 ~ value + Infant_age_days, data=dat, family = "binomial")
  sumfit <- summary(fit)
  pval <- sumfit[["coefficients"]][11]
  as.data.frame(pval)
}

pval <- long %>%
  purrr::map(tmpfn) %>%
  bind_rows(.id = "HMO")

pval <- pval %>% separate(HMO, into = c("HMO", "Condition"), sep = " ")
pval <- pval %>% 
  spread(Condition, pval)
pval <- pval %>% 
  mutate_if(is.numeric, ~ p.adjust(.x, method = "fdr", n = length(.x)))
colnames(pval)[-1] <- paste("padj", colnames(pval)[-1], sep = "_")

######### Make plot ###################################### 
# Make correlation plot
fold_change <- fold_change %>% 
  gather(HMO, value, 2:ncol(fold_change)) %>%
  spread(Condition, value)
fold_change <- fold_change %>% 
  filter(!str_detect(HMO, "Infant"))
fold_change <- left_join(fold_change, pval, by = "HMO")
fold_change <- fold_change %>% 
  dplyr::select(HMO, C_section, Migraines, Pain_reliever,  Preterm, SSRI, Alcohol, Anxiety, Depression, EBF, Opioids, 
                padj_C_section, padj_Migraines, padj_Pain_reliever, padj_Preterm, padj_SSRI, padj_Alcohol, padj_Anxiety, padj_Depression, padj_EBF, padj_Opioids)
fold_change$HMO <- str_replace_all(fold_change$HMO, "X", "")
fold_change <- fold_change %>% arrange(HMO)

fc1 <- fold_change %>% dplyr::select(HMO, !contains("padj"))
fc2 <- fold_change %>% dplyr::select(HMO, contains("padj"))
fc1 <- column_to_rownames(fc1, var="HMO")
fc2 <- column_to_rownames(fc2, var = "HMO")
#fc1[fc1>2] <- 1.9
fc1 <- fc1*-1
#fc1[fc1>2] <- 1.9
corrplot(as.matrix(fc1), is.corr = FALSE, cl.pos = 'r', cl.ratio = .2, addgrid.col = NA, p.mat = as.matrix(fc2),  
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9, insig = 'label_sig', pch.col = 'grey20', col.lim = c(-4,4), cl.length = 5)


rm(long, tmpfn, pval, stdev, average, groups, fc1, fc2)

data1 <- data %>% 
  group_by(Secretor) %>% 
  summarise_all(~sum(!is.na(.)))

############## Look at Four amigos in metadata over time ###############################################
# Make HMOs categorical variable based on child age
data <- Amigos %>% dplyr::select(-Secretor)
data <- full_join(HMO, data, by = "Sample")
data1 <- data %>% filter(!is.na(Preterm)) 

data1 <- data1 %>% filter(Secretor ==1)
data1$Infant_age_days <- log2(data1$Infant_age_days)

# Need to change color here!
ggscatter(data1, y = "z6DDL", x = "Infant_age_days", add = "reg.line",  conf.int = TRUE,
          cor.method = "spearman", shape = "Dataset", color = "C_section", ggtheme = theme_bw()) + 
  stat_cor(aes(color = C_section)) + facet_wrap(~Secretor)

fit <- data1 %>% group_by(Secretor, C_section) %>% rstatix::cor_test(vars=c("Infant_age_days"), vars2=c("z6DDL"), method="spearman")
fit$padj <- p.adjust(fit$p, method = "BH")

