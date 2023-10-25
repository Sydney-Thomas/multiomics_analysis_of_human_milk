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
library(compositions)
library(rstatix)
library(FactoMineR)

metadata <- read.csv("./Data/US/Metadata_SEED.csv")

names <- c("Opioids", "Secretor", "Child_1_Delivery", "Child_1_Preterm_Birth", "Antibiotic", "SSRI", "Allergy", "Pain_reliver", "Caffeine", "Alcohol",
           "Breastfeeding_Experience", "Marijuana", "Depression", "Anxiety", "Asthma", "Diabetes", "Migraines", "Anemia", "HPV", "Breastfeeding_frequency_factor", "EBF")
metadata <- metadata %>% mutate(across(any_of(names), as.factor))
rm(names)

## Normalization ######################################################################################################
## This step performs rclr transformation and filters metabolites that are found in less than 10% of samples
## These thresholds can be changed depending on the dataset
## Clean up dataframe from Mzmine
norm <- read.csv("C:/Users/sthomas/OneDrive - University of California, San Diego Health/SEED Grant/Original_Data/Original_Data_SEED/SEED_Grant_iimn_gnps_quant.csv")
norm <- norm %>% dplyr::select(-X)
norm <- norm %>% rename_with(~str_replace_all(., ".mzML.Peak.area", ".M"))
norm <- norm %>% dplyr::select(-(row.m.z:correlation.group.ID), -(best.ion:neutral.M.mass))
norm <- norm %>% mutate(across(row.ID:annotation.network.number, as.character))

## Merge ions from the same molecule (iimn)
iimn <- norm %>%
  group_by(annotation.network.number) %>%
  summarise(across(where(is.numeric), ~sum(., na.rm = TRUE))) %>%
  filter(!is.na(annotation.network.number)) %>%
  rename(row.ID = annotation.network.number)
labels <- norm %>% dplyr::select(row.ID, annotation.network.number)
iimn$row.ID <- paste0(iimn$row.ID, "_i")
norm_iimn <- norm %>% filter(is.na(annotation.network.number)) %>% dplyr::select(-annotation.network.number)
norm_iimn <- rbind(norm_iimn, iimn)

## Read in FBMN library IDs from GNPS
library_ID <- read.delim("./Data/US/FBMN_IDs.tsv")
library_ID <- library_ID %>% dplyr::select(cluster.index, LibraryID, precursor.mass, componentindex)
colnames(library_ID) <- c("row.ID", "LibraryID", "Precursor_Mass", "Network_Number")
library_ID$row.ID <- as.character(library_ID$row.ID)
library_ID <- right_join(labels, library_ID, by = "row.ID", multiple = "all")
library_ID[library_ID == "N/A"] <- NA
library_ID$Network_Number[library_ID$Network_Number == -1] <- NA

## Collapse ions (keeps all IDs, precursor masses, and network numbers of features in the iimn network separated by OR)
result_list <- list()
for (colname in c("LibraryID", "Precursor_Mass", "Network_Number")) {
  result <- library_ID %>%
    filter(!is.na(annotation.network.number)) %>%
    group_by(annotation.network.number, !!sym(colname)) %>%
    summarise(n = n()) %>%
    group_by(annotation.network.number) %>%
    summarise(concat = paste(na.omit(!!sym(colname)), collapse = " OR ")) %>%
    rename(!!colname := "concat")
  result_list[[colname]] <- result
}

library_iimn <- result_list %>% reduce(full_join, by = "annotation.network.number") %>% rename("row.ID" = "annotation.network.number")
library_iimn$row.ID <- paste0(library_iimn$row.ID, "_i")
library_f <- library_ID %>% filter(is.na(annotation.network.number)) %>% dplyr::select(-annotation.network.number)
library_iimn <- rbind(library_f, library_iimn)
library_iimn <- library_iimn %>% rename(Metabolite = row.ID)
library_iimn[library_iimn == ""] <- NA
write_csv(library_iimn, "library_matches.csv")

## Read in canopus from sirius workflow
canopus <- read.delim("canopus_compound_summary.tsv")
canopus$name <- str_replace_all(canopus$id, ".*iimn_", "")
canopus <- canopus %>% dplyr::select(name, NPC.superclass, NPC.class, NPC.pathway, ClassyFire.superclass, ClassyFire.class, ClassyFire.subclass, ClassyFire.level.5)
canopus[canopus == ""] <- NA
## Remove duplicates
canopus <- canopus %>% distinct(name, .keep_all = TRUE)
canopus <- right_join(labels, canopus, by = c("row.ID" = "name"))

## Collapse ions (keeps all canopus IDs in the iimn network separated by OR)
result_list <- list()
for (colname in c("NPC.superclass", "NPC.class", "NPC.pathway", "ClassyFire.superclass", "ClassyFire.class", "ClassyFire.subclass", "ClassyFire.level.5")) {
  result <- canopus %>%
    filter(!is.na(annotation.network.number)) %>%
    group_by(annotation.network.number, !!sym(colname)) %>%
    summarise(n = n()) %>%
    group_by(annotation.network.number) %>%
    summarise(concat = paste(na.omit(!!sym(colname)), collapse = " OR ")) %>%
    rename(!!colname := "concat")
  result_list[[colname]] <- result
}
canopus_iimn <- result_list %>% reduce(full_join, by = "annotation.network.number") %>% rename("row.ID" = "annotation.network.number")
canopus_iimn$row.ID <- paste0(canopus_iimn$row.ID, "_i")
canopus_f <- canopus %>% filter(is.na(annotation.network.number)) %>% dplyr::select(-annotation.network.number)
canopus_iimn <- rbind(canopus_f, canopus_iimn)
canopus_iimn <- canopus_iimn %>% rename(Metabolite = row.ID)
write_csv(canopus_iimn, "canopus_matches.csv")

norm_iimn <- norm_iimn %>%
  gather(Sample, value, 2:ncol(norm_iimn)) %>%
  spread(row.ID, value)
norm_iimn$Sample <- str_replace_all(norm_iimn$Sample, "^X", "")
norm_iimn <- column_to_rownames(norm_iimn, var = "Sample")
## Perform rclr preprocessing
norm_iimn <- norm_iimn %>% decostand(method = "rclr")
## Filter metabolites that are found in at least 10% of samples
norm_iimn <- norm_iimn %>% dplyr::select(where(~sum(. != 0) >= (0.1*nrow(norm_iimn))))
## Missing value imputation
norm_iimn <- rownames_to_column(norm_iimn, var = "Sample")
set.seed(141)
norm_imp <- norm_iimn
for (i in 2:ncol(norm_imp)) {
  norm_imp[,i] <- ifelse(norm_imp[,i] == 0,
                         round(runif(nrow(norm_imp), min = min(norm_imp[,i])-1, max = min(norm_imp[,i])), 1),
                         norm_imp[,i])
}
write_csv(norm_iimn, "Metabolites_normalized.csv")
write_csv(norm_imp, "Metabolites_normalized_imputed.csv")
write_csv(labels, "iimn_groups.csv")
rm(i, norm, norm_imp, labels, canopus, canopus_iimn, library_ID, library_iimn, result, result_list, iimn, colname, canopus_f, library_f, norm_iimn)

#############################################################################
# Read in normalized data (if you've already performed normalization above)
normalized <- read.csv("Metabolites_normalized_imputed.csv")
library_ID <- read.csv("library_matches.csv")
canopus <- read.csv("canopus_matches.csv")
data <- left_join(normalized, metadata, by = "Sample")

######## Do PCA plot of Secretor status ###############################
pca_sec <- data %>% dplyr::select(Sample, Secretor, starts_with("X")) %>% filter(!is.na(Secretor))
pca <- PCA(pca_sec[,3:ncol(pca_sec)], scale.unit = FALSE, ncp = 5, graph = FALSE)
fviz_pca_ind(pca, axes = c(1,2), habillage = as.factor(pca_sec$Secretor), label = "none", addEllipses = TRUE, ellipse.level = 0.95, invisible="quali")
fviz_pca_biplot(pca, axes = c(1,2), habillage = as.factor(pca_sec$Secretor), addEllipses = TRUE, select.var = list(contrib = 5), label = "var", invisible="quali", col.var = "black")
pca_mat <- pca_sec %>% dplyr::select(-Sample, -Secretor) %>% as.matrix()
dist <- vegdist(pca_mat, method = "euclidean")
permanova <- adonis2(dist ~ as.factor(pca_sec$Secretor), pca_sec, na.action = na.omit, method = "euclidean")
rm(pca_sec, pca, pca_mat, dist, permanova)

########### ANOVAS ################################################################################################################
########### Calculate average values for each metabolite with HMO production ############################
groups <- data %>% 
  dplyr::select(Sample, Amigos, starts_with("X")) %>% 
  filter(!is.na(Amigos))

# Now you can start doing calculations. First, average all the replicates in each sample set 
average <- groups %>% 
  group_by(Amigos) %>% 
  summarise_if(is.numeric, mean, na.rm = TRUE)
# You'll get a warning on this command, it just means that R couldn't average your sample names, which is probably good since they're letter and not numbers

# Next calculate standard deviations for each set of replicates
stdev <- groups %>% 
  group_by(Amigos) %>% 
  summarise_if(is.numeric, sd, na.rm = TRUE)

# Run tukey_HSD
pval <- groups %>% pivot_longer(4:ncol(groups), names_to = "Metabolite", values_to = "Abundance")
pval <- pval %>% dplyr::select(-Sample) %>% group_by(Metabolite) %>% tukey_hsd(Abundance~Amigos)
pval$p.adj.BH <- p.adjust(pval$p.adj, method = "BH")

pval$Metabolite <- str_replace(pval$Metabolite, "X", "")
pval <- left_join(pval, library_ID, by = "Metabolite", multiple = "all")
pval <- left_join(pval, canopus, by = "Metabolite", multiple = "all")
write.csv(pval, "Metabolite_Changes.csv", row.names = F)

# Look at individual metadata 
ggplot(groups, aes(x = Amigos, y = X276)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=0.9, shape=16) +
  facet_wrap(~Secretor)

########## END ANOVAS #####################################################################################

########## CORRELATIONS ###################################################################################
# Compare to HMOs #########################################################################################
# Use spearman correlation since HMOs aren't evenly distributed ###########################################
data <- left_join(normalized, metadata, by = "Sample")
corr <- data %>% dplyr::select(Sample, z6DDL, starts_with("X")) 
corr <- column_to_rownames(corr, var = "Sample")
cor <- rcorr(as.matrix(corr), type=c("spearman"))

# Extract coefficients
coeff <- as.data.frame(cor$r)
coeff <- rownames_to_column(coeff, var="Metabolite")
coeff <- coeff %>% 
  dplyr::select(Metabolite, z6DDL) %>% 
  filter(!str_detect(Metabolite, "z6DDL"))
coeff <- coeff %>% 
  filter(abs(z6DDL) > 0.0)

# Extract p values
sig <- as.data.frame(cor$P)
sig <- rownames_to_column(sig, var="Metabolite")
sig <- sig %>% 
  dplyr::select(Metabolite, z6DDL) %>% 
  filter(!str_detect(Metabolite, "z6DDL")) 
sig <- sig %>% 
  filter(Metabolite %in% coeff$Metabolite)
sig$padj <- p.adjust(sig$z6DDL, method = "BH")
sig <- as.data.frame(sig)
colnames(sig) <- c("Metabolite", "p", "padj")
coeff <- left_join(coeff, sig, by = "Metabolite")

coeff$Metabolite <- str_replace(coeff$Metabolite, "X", "")
coeff <- left_join(coeff, library_ID, by = "Metabolite", multiple = "all")
coeff <- left_join(coeff, canopus, by = "Metabolite", multiple = "all")

write.csv(coeff, "HMO_Metabolite_correlations.csv", row.names = F)

rm(corr)
rm(cor)
rm(sig)

##########################################################################################################################################################
####################### CANOPUS ##########################################################################################################################
data <- read.csv("Metabolites_normalized.csv")
data <- data %>% 
  gather(Metabolite, value, 2:ncol(data)) %>% 
  spread(Sample, value)
data[data == 0] <- NA
data$Metabolite <- str_replace(data$Metabolite, "X", "")
data <- left_join(data, canopus, by = "Metabolite", multiple = "all")

# Calculate mean value of each class, superclass, subclass, etc
ClassyFire.superclass <- data %>%
  dplyr::select(ClassyFire.superclass, where(is.numeric)) %>%
  group_by(ClassyFire.superclass) %>%
  filter(ClassyFire.superclass != "") %>%
  summarise_all(mean, na.rm = TRUE)
names(ClassyFire.superclass)[1] <- "canopus"
ClassyFire.superclass$Type <- "ClassyFire.superclass"
ClassyFire.superclass$canopus <- paste(as.character(ClassyFire.superclass$canopus), "_Classy", sep = "")
ClassyFire.class <- data %>% 
  dplyr::select(ClassyFire.class, where(is.numeric)) %>% 
  group_by(ClassyFire.class) %>% 
  filter(ClassyFire.class != "") %>% 
  summarise_all(mean, na.rm = TRUE) 
names(ClassyFire.class)[1] <- "canopus"
ClassyFire.class$Type <- "ClassyFire.class"
ClassyFire.class$canopus <- paste(as.character(ClassyFire.class$canopus), "_Classy", sep = "")
ClassyFire.subclass <- data %>% 
  dplyr::select(ClassyFire.subclass, where(is.numeric)) %>% 
  group_by(ClassyFire.subclass) %>% 
  filter(ClassyFire.subclass != "") %>% 
  summarise_all(mean, na.rm = TRUE) 
names(ClassyFire.subclass)[1] <- "canopus"
ClassyFire.subclass$Type <- "ClassyFire.subclass"
ClassyFire.subclass$canopus <- paste(as.character(ClassyFire.subclass$canopus), "_Classy", sep = "")
ClassyFire.level.5 <- data %>% 
  dplyr::select(ClassyFire.level.5, where(is.numeric)) %>% 
  group_by(ClassyFire.level.5) %>% 
  filter(ClassyFire.level.5 != "") %>% 
  summarise_all(mean, na.rm = TRUE) 
names(ClassyFire.level.5)[1] <- "canopus"
ClassyFire.level.5$Type <- "ClassyFire.level.5"
ClassyFire.level.5$canopus <- paste(as.character(ClassyFire.level.5$canopus), "_Classy", sep = "")

canopus_all <- rbind(ClassyFire.superclass, ClassyFire.class, ClassyFire.subclass, ClassyFire.level.5)
canopus_all$canopus <- str_replace_all(canopus_all$canopus, " |,|-", "_")
write_csv(canopus_all, "canopus_SEED.csv")

rm(ClassyFire.superclass, ClassyFire.class, ClassyFire.subclass, ClassyFire.level.5)

########## Look at subclass correlations to HMOs ####################################################
data <- canopus_all %>% 
  dplyr::select(-Type)
data <- data %>% 
  gather(Sample, value, 2:ncol(data)) %>% 
  spread(canopus, value)
data <- left_join(data, metadata, by = "Sample")
corr <- data %>% 
  dplyr::select(Sample, z6DDL, 2:Tertiary_amines_Classy) 
corr <- column_to_rownames(corr, var = "Sample")
cor <- rcorr(as.matrix(corr), type=c("spearman"))

# Extract coefficients
coeff <- as.data.frame(cor$r)
coeff <- rownames_to_column(coeff, var="Metabolite")
coeff <- coeff %>% 
  dplyr::select(Metabolite, z6DDL) %>% 
  filter(!str_detect(Metabolite, "z6DDL")) 
coeff <- coeff %>% 
  filter(abs(z6DDL) > 0.0)


# Extract p values
sig <- as.data.frame(cor$P)
sig <- rownames_to_column(sig, var="Metabolite")
sig <- sig %>% 
  dplyr::select(Metabolite, z6DDL) %>% 
  filter(!str_detect(Metabolite, "z6DDL")) %>% 
  arrange(z6DDL)
sig <- sig %>% 
  filter(Metabolite %in% coeff$Metabolite)
sig$padj <- p.adjust(sig$z6DDL, method = "BH")
sig <- as.data.frame(sig)
colnames(sig) <- c("Metabolite", "p", "padj")
coeff <- left_join(coeff, sig, by = "Metabolite")
write.csv(coeff, "canopus_HMO_correlations.csv", row.names = F)

rm(corr)
rm(cor)
rm(sig)

# Correlation plots
data <- canopus_all %>% 
  dplyr::select(canopus, Type)
data$Type <- factor(data$Type, levels = c("ClassyFire.superclass", "ClassyFire.class", "ClassyFire.subclass", "ClassyFire.level.5"))
names(data)[1] <- "Metabolite"
coeff <- left_join(coeff, data, by = "Metabolite")
ggplot(coeff, aes(x=z6DDL, y=reorder(Metabolite, z6DDL))) +
  geom_bar(stat="identity", aes(fill=z6DDL)) + 
  scale_fill_distiller(palette="RdBu") +
  facet_rep_grid(Type~., scales = "free", space = "free", repeat.tick.labels = T)



  