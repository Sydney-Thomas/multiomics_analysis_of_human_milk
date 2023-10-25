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
library(ggExtra)
library(ggpubr)
library(rstatix)
library(igraph)
library(ggraph)
library(factoextra)
library(FactoMineR)
library(vegan)

##### Compare phylum-level data from three datasets ###################################
Phylum_SEED <- read.csv("./Data/US/phylum_table_SEED.csv")
colnames(Phylum_SEED) <- c("Phylum", "US")
Phylum_SEED$Phylum <- str_replace_all(Phylum_SEED$Phylum, " ", "")
Phylum_Collado <- read.csv("./Data/ES/phylum_table_Collado.csv")
colnames(Phylum_Collado) <- c("Phylum", "ES")
Phylum_Collado$Phylum <- str_replace_all(Phylum_Collado$Phylum, "Campylobacterota", "Campilobacterota")
Phylum_Geddes <- read.csv("./Data/AU/phylum_table_Geddes.csv")
colnames(Phylum_Geddes) <- c("Phylum", "AU")

Phylum <- full_join(Phylum_SEED, Phylum_Collado, by = "Phylum")
Phylum <- full_join(Phylum, Phylum_Geddes, by = "Phylum")

# Remove any phyla that were only found twice in each dataset
Phylum <- Phylum %>%
  mutate(Sum = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>%
  filter(Sum > 2) %>%
  dplyr::select(-Sum)
Phylum <- pivot_longer(Phylum, 2:ncol(Phylum), names_to = "Dataset", values_to = "Count")
colorcount <- length(unique(Phylum$Phylum))
getPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
ggplot(Phylum, aes(x = Dataset, y = Count, fill = Phylum)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = getPalette(colorcount))

rm(Phylum, Phylum_Collado, Phylum_Geddes, Phylum_SEED)
rm(colorcount, getPalette)

### Look at changes in Phylum abundance over time #####################################################################
### Only use US data here, the rest don't have enough time points #####################################################
US_phyla <- read.csv("./Data/US/Phylum_abundance_SEED.csv")
metadata_SEED <- read.csv("./Data/US/Metadata_SEED.csv")
rnm <- c(Infant_age_days = "Child_1_Age", Infant_sex = "Child_1_Gender", Infant_sex = "infant_sex", C_section = "Child_1_Delivery", Preterm = "Child_1_Preterm_Birth", Parent_BM<- = "Mother_BMI",
         Parent_Age = "Age", Parent_Age = "Mother_age", Ethnicity = "Race", Ethnicity = "ethnicity", Parity = "parity", Parent_BMI = "Pre_BMI")
metadata_SEED <- metadata_SEED %>% rename(any_of(rnm))

US_phyla <- US_phyla %>% pivot_longer(2:ncol(US_phyla), names_to = "Phylum", values_to = "Abundance")
US_phyla <- left_join(US_phyla, metadata_SEED, by = "Sample")
rm(metadata_SEED, rnm)

US_phyla$Amigos <- factor(US_phyla$Amigos, levels = c("Low", "Medium", "High"))
US_phyla$Dataset <- "US"
US_phyla <- US_phyla %>% dplyr::select(Sample, Phylum, Abundance, Infant_age_days, Amigos, Dataset)
Phyla <- US_phyla
rm(US_phyla)

Phyla <- Phyla %>% filter(!str_detect(Phylum, "Bacteria_unclassified"))
top <- Phyla %>% filter(Abundance != 0) %>% count(Phylum) %>% slice_max(n, n = 8)
Phyla <- Phyla %>% filter(Phylum %in% top$Phylum) %>% filter(Infant_age_days < 1000)

ggplot(Phyla, aes(x=Infant_age_days, y= Abundance, color=Phylum)) +
  geom_smooth(method = "loess") +
  theme_minimal()

####### Compare genus-level data from three datasets ######################################
Genus_SEED <- read.csv("./Data/US/Microbiome_table_SEED_genus_rclr.csv")
Genus_Collado <- read.csv("./Data/ES/Microbiome_table_Collado_genus_rclr.csv")
Genus_Geddes <- read.csv("./Data/AU/Microbiome_table_Geddes_genus_rclr.csv")

metadata_SEED <- read.csv("./Data/US/Metadata_SEED.csv")
metadata_Coll <- read.csv("./Data/ES/Metadata_Collado.csv")
metadata_Geddes <- read.csv("./Data/AU/Metadata_Geddes.csv")

rnm <- c(Infant_age_days = "Child_1_Age", Infant_sex = "Child_1_Gender", Infant_sex = "infant_sex", C_section = "Child_1_Delivery", Preterm = "Child_1_Preterm_Birth", Parent_BMI = "Mother_BMI",
         Parent_Age = "Age", Parent_Age = "Mother_age", Ethnicity = "Race", Ethnicity = "ethnicity", Parity = "parity", Parent_BMI = "Pre_BMI")
metadata_SEED <- metadata_SEED %>% rename(any_of(rnm))
metadata_Coll <- metadata_Coll %>% rename(any_of(rnm))
metadata_Geddes <- metadata_Geddes %>% rename(any_of(rnm))

Genus_SEED <- Genus_SEED %>% pivot_longer(2:ncol(Genus_SEED), names_to = "Genus", values_to = "Abundance")
Genus_Collado <- Genus_Collado %>% pivot_longer(2:ncol(Genus_Collado), names_to = "Genus", values_to = "Abundance")
Genus_Geddes <- Genus_Geddes %>% pivot_longer(2:ncol(Genus_Geddes), names_to = "Genus", values_to = "Abundance")

Genus_Collado <- left_join(Genus_Collado, metadata_Coll, by = "SampleID") %>% dplyr::select(-Sample) %>% rename(Sample = SampleID)
Genus_SEED <- left_join(Genus_SEED, metadata_SEED, by = "Sample")
Genus_Geddes <- left_join(Genus_Geddes, metadata_Geddes, by = "Study_ID") %>% dplyr::select(-Sample) %>% rename(Sample = Study_ID)
rm(metadata_SEED, metadata_Coll, metadata_Geddes, rnm)

Genus_Collado$Amigos <- factor(Genus_Collado$Amigos, levels = c("Low", "Medium", "High"))
Genus_SEED$Amigos <- factor(Genus_SEED$Amigos, levels = c("Low", "Medium", "High"))
Genus_Geddes$Amigos <- factor(Genus_Geddes$Amigos, levels = c("Low", "Medium", "High"))

Genus_SEED$Dataset <- "US"
Genus_Collado$Dataset <- "ES"
Genus_Geddes$Dataset <- "AU"
Genus_Geddes <- Genus_Geddes %>% filter(!str_detect(Genus, "_"))

# Look at specific genera ##########################################################################################
j = "othia|eillone"

Genus_1 <- Genus_SEED %>% filter(str_detect(Genus, j)) %>% dplyr::select(Sample, Genus, Abundance, Secretor, Amigos, Dataset, Infant_age_days, Preterm, C_section)
Genus_2 <- Genus_Collado %>% filter(str_detect(Genus, j)) %>% dplyr::select(Sample, Genus, Abundance, Secretor, Amigos, Dataset, Infant_age_days, Preterm, C_section)
Genus_3 <- Genus_Geddes %>% filter(str_detect(Genus, j)) %>% dplyr::select(Sample, Genus, Abundance, Secretor, Amigos, Dataset, Infant_age_days, Preterm, C_section)
Genus_3 <- Genus_3 %>% filter(!str_detect(Genus, "_unclass"))
Genus_plot <- rbind(Genus_1, Genus_2, Genus_3)
Genus_plot$Abundance[Genus_plot$Abundance == 0] <- NA
Perc_zero <- Genus_plot %>% pivot_wider(names_from = Genus, values_from = Abundance)
Perc_zero$Category <- ifelse((is.na(Perc_zero$Rothia)&is.na(Perc_zero$Veillonella)),"zNeither", ifelse((is.na(Perc_zero$Rothia)&!is.na(Perc_zero$Veillonella)), "Veillonella",
                             ifelse((!is.na(Perc_zero$Rothia)&is.na(Perc_zero$Veillonella)), "Rothia", "Both")))
Perc_zero <- Perc_zero %>% group_by(Amigos) %>% filter(!is.na(Amigos)) %>% count(Category)
ggplot(Perc_zero, aes(x=Amigos, y=n, fill=Category)) +
  geom_bar(stat="identity", position = "fill")

Genus_plot <- Genus_plot %>% filter(!is.na(Abundance)) %>% filter(!is.na(Secretor))
rm(Genus_1, Genus_2, Genus_3)

comparisons <- list(c("Low", "High"), c("Medium", "High"), c("Low", "Medium"))
ggplot(Genus_plot, aes(x = Amigos, y = Abundance, color = Amigos)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape=Dataset)) +
  stat_compare_means(comparisons = comparisons) +
  theme_minimal() +
  #theme(legend.position = "none") +
  facet_wrap(~Genus)

comparisons <- list(c("AU", "ES"), c("ES", "US"), c("AU", "US"))
ggplot(Genus_plot, aes(x = Dataset, y = Abundance, color = Dataset)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape=Dataset)) +
  stat_compare_means(comparisons = comparisons) +
  theme_minimal() +
  #theme(legend.position = "none") +
  facet_wrap(~Genus)

Genus_plot <- Genus_plot %>% pivot_wider(names_from = Genus, values_from = Abundance)
ggplot(Genus_plot, aes(x=Veillonella, y= Rothia, color=Dataset)) +
  geom_point(aes(shape = Dataset)) +
  geom_smooth(method = "glm") +
  stat_cor(aes(color = Dataset)) +
  theme_minimal()

# Make HMOs categorical variable based on child age
Genus_plot <- Genus_plot %>% rowwise() %>% mutate(Total = mean(c_across(c('Rothia', 'Veillonella')), na.rm = TRUE))
Genus_plot_1 <- Genus_plot %>% filter(Infant_age_days < 1000)

# # Choose metadata you want to compare ########################################################################################################################
i <- "C_section"

Genus_Collado <- Genus_Collado %>% filter(!is.na((!!sym(i)))) %>% dplyr::select(Sample, Genus, Abundance, (!!sym(i)))
Genus_SEED <- Genus_SEED %>% filter(!is.na((!!sym(i)))) %>% dplyr::select(Sample, Genus, Abundance, (!!sym(i)))
Genus_Geddes <- Genus_Geddes %>% filter(!is.na((!!sym(i)))) %>% dplyr::select(Sample, Genus, Abundance, (!!sym(i)))

Genus_Collado$Abundance[Genus_Collado$Abundance == 0] <- NA
Ct <- Genus_Collado %>% group_by_at(vars(i, Genus)) %>% summarise(sm = sum(!is.na(Abundance)))  %>% filter(sm >= 5) %>% ungroup() %>% dplyr::select(-(!!sym(i)), -sm)
Ct <- subset(Ct, duplicated(Ct))
Genus_Collado <- Genus_Collado %>% filter(Genus %in% Ct$Genus)

Genus_SEED$Abundance[Genus_SEED$Abundance == 0] <- NA
St <- Genus_SEED %>% group_by_at(vars(i, Genus)) %>% summarise(sm = sum(!is.na(Abundance)))  %>% filter(sm >= 5) %>% ungroup() %>% dplyr::select(-(!!sym(i)), -sm)
St <- subset(St, duplicated(St))
Genus_SEED <- Genus_SEED %>% filter(Genus %in% St$Genus)

Genus_Geddes$Abundance[Genus_Geddes$Abundance == 0] <- NA
Gt <- Genus_Geddes %>% group_by_at(vars(i, Genus)) %>% summarise(sm = sum(!is.na(Abundance)))  %>% filter(sm >= 5) %>% ungroup() %>% dplyr::select(-(!!sym(i)), -sm)
Gt <- subset(Gt, duplicated(Gt))
Genus_Geddes <- Genus_Geddes %>% filter(Genus %in% Gt$Genus)
rm(Gt, St, Ct)

Genus_SEED$Dataset <- "US"
Genus_Collado$Dataset <- "ES"
Genus_Geddes$Dataset <- "AU"

Genus_C <- Genus_Collado %>% group_by_at(vars(i, Genus, Dataset)) %>% summarise_at(c("Abundance"), mean, na.rm = TRUE)
Genus_S <- Genus_SEED %>% group_by_at(vars(i, Genus, Dataset)) %>% summarise_at(c("Abundance"), mean, na.rm = TRUE)
Genus_G <- Genus_Geddes %>% group_by_at(vars(i, Genus, Dataset)) %>% summarise_at(c("Abundance"), mean, na.rm = TRUE)

Genus <- rbind(Genus_C, Genus_S, Genus_G)
Genus_raw <- rbind(Genus_Collado, Genus_SEED, Genus_Geddes)

Genus_raw[i] <- as.factor(Genus_raw[[i]])

rm(Genus_C, Genus_G, Genus_S, Genus_Collado, Genus_SEED, Genus_Geddes)

# Run PCA #########################################
pca_data <- Genus_raw %>% pivot_wider(names_from = "Genus", values_from = "Abundance")
pca_data[is.na(pca_data)] <- 0


pca_microbes <- list()  # new empty list
for (n in c("AU", "ES", "US")) {
  # Filter data
  pca_filtered <- pca_data %>% filter(str_detect(Dataset, n))
  pca <- PCA(pca_filtered[,4:ncol(pca_filtered)], scale.unit = FALSE, ncp = 5, graph = FALSE)
  plot <- fviz_pca_biplot(pca, axes = c(1,2), habillage = pca_filtered[[i]], addEllipses = TRUE, select.var = list(contrib = 5), label = "var", invisible="quali", col.var = "black")
  pca_microbes[[n]] <- plot
}

patchwork::wrap_plots(pca_microbes, ncol = 3)

permanova <- list()  # new empty list
for (n in c("AU", "ES", "US")) {
  # Filter data
  pca_filtered <- pca_data %>% filter(str_detect(Dataset, n))
  pca_matrix <- pca_filtered %>% dplyr::select(-Dataset, -Sample, -Amigos) %>% as.matrix()
  dist <- vegdist(pca_matrix, method = "euclidean")
  permanova_genus <- adonis2(dist ~ pca_filtered[[i]], pca_filtered, na.action = na.omit, method = "euclidean")
  permanova[[n]] <- permanova_genus
}

permanova_microbes <- do.call(rbind, permanova)
rm(pca_filtered, pca_microbes, pca_var, plot, pca, pca_data, permanova_genus, permanova, permanova_microbes)

# Calculate p values ############################################################################
fm <- as.formula(paste("Abundance", "~", i))
pval <- Genus_raw %>% group_by(Dataset, Genus) %>% tukey_hsd(fm)
rm(fm)

# Now calculate fold changes (TRUE/FALSE)
fold_change <- Genus %>% dplyr::select(Genus, Dataset, (!!sym(i)), Abundance) %>% pivot_wider(names_from = (!!sym(i)), values_from = Abundance)
colnames(fold_change) <- c("Genus", "Dataset", "false", "true")
fold_change <- fold_change %>% dplyr::mutate(fc = false - true)

# Look at all datasets together
Total <- full_join(fold_change, pval, by = c("Genus", "Dataset"))
Total <- Total %>% filter(p.adj < 0.1)
Total$sig <- ifelse(Total$p.adj < 0.001, 1, ifelse(Total$p.adj<0.05, 2, 3))
# If the number of rows in network_table is below the number of nodes, then the graph will give an error
# You can artificially add rows to your network table by specifying n value here
Total_max <- Total %>% ungroup() %>% slice_max(abs(fc), n=10)
Total_max$fc <- Total_max$fc*-1
Total_max$Dataset <- "Extra"
Total_max$sig <- 1
Total <- rbind(Total, Total_max)
rm(Total_max)

# Make graph #############
# Generate nodes ###########################################
nodes_HMO <- Total %>% ungroup() %>% dplyr::select(Dataset) %>% unique() 
nodes_HMO <- data.frame(node = nodes_HMO, category = rep("Dataset",length(nodes_HMO)))
names(nodes_HMO)[1] <- "node"
nodes_bac <- Total %>% ungroup() %>% dplyr::select(Genus) %>% unique()
nodes_bac <- data.frame(node = nodes_bac, category = rep("Bacteria",length(nodes_bac)))
names(nodes_bac)[1] <- "node"
nodes <- rbind(nodes_HMO, nodes_bac)
node_colors <- c("#83A4C9", "#5F8AB9")
names(node_colors) <- c("Dataset", "Bacteria")
rm(nodes_bac, nodes_HMO)

# Generate network table
network_table <- Total %>% ungroup() %>% dplyr::select(Dataset, Genus, fc, everything())
#network_table <- network_table[complete.cases(network_table),]
colnames(network_table)[1:3] <- c("from", "to", "correlation")
network_table$direction <- ifelse(network_table$correlation < 0, "neg", "pos")
network_table$weight <- abs(network_table$correlation)
# Generate network
network_seed <- graph_from_data_frame(d = network_table, vertices = nodes, directed = FALSE)
SEED_Network <- ggraph(network_seed) +
  geom_edge_link(aes(color = correlation, linetype = as.factor(sig))) +
  scale_edge_color_gradientn(colors = brewer.pal(n=8, name="RdBu")) +
  geom_node_point(aes(shape = category), size = 3) +
  theme_void() +
  geom_node_text(aes(label = name), size = 3, repel = TRUE) +
  labs(shape = "Category", edge_color = "Fold Change")
SEED_Network

###############################################################################################################################################################
###################################### Look at correlations between HMOs and microbiome between datasets ########################################################
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

######### Read in HMO data ################################################
HMO <- read.csv("./Data/HMO_scaled_centered.csv")
HMO <- HMO %>% dplyr::select(Sample:LSTb, z6DDL, Secretor, Amigos, Dataset)

Microbiome_SEED <- read.csv("./Data/US/Microbiome_table_SEED_genus_rclr.csv")
Microbiome_Coll <- read.csv("./Data/ES/Microbiome_table_Collado_genus_rclr.csv")
Microbiome_Geddes <- read.csv("./Data/AU/Microbiome_table_Geddes_genus_rclr.csv")
Microbiome_SEED[Microbiome_SEED == 0] <- NA
Microbiome_Coll[Microbiome_Coll == 0] <- NA
Microbiome_Geddes[Microbiome_Geddes == 0] <- NA

Microbiome_SEED <- Microbiome_SEED %>% 
  gather(Microbiome, value, 2:ncol(Microbiome_SEED)) %>%
  spread(Sample, value)
Microbiome_Coll <- Microbiome_Coll %>% 
  gather(Microbiome, value, 2:ncol(Microbiome_Coll)) %>%
  spread(SampleID, value)
Microbiome_Geddes <- Microbiome_Geddes %>% 
  gather(Microbiome, value, 2:ncol(Microbiome_Geddes)) %>%
  spread(Study_ID, value)
Microbiome_Geddes <- Microbiome_Geddes %>% filter(!str_detect(Microbiome, "_"))

Microbiome <- full_join(Microbiome_SEED, Microbiome_Geddes, by = "Microbiome")
Microbiome <- full_join(Microbiome, Microbiome_Coll, by = "Microbiome")

Microbiome <- Microbiome %>% 
  gather(Sample, value, 2:ncol(Microbiome)) %>%
  spread(Microbiome, value)

Microbiome <- full_join(HMO, Microbiome, by = "Sample")
Microbiome <- Microbiome %>% filter(!is.na(z6DDL))
Microbiome <- Microbiome %>% dplyr::select(Sample, Dataset, Secretor, Amigos, everything()) %>% mutate(across(5:ncol(Microbiome), as.numeric))

rm(Microbiome_Coll, Microbiome_Geddes, Microbiome_SEED, HMO)

# Individual correlation matrices
corr_hmos <- list()  # new empty list
for (n in c("AU", "US", "ES")) {
  data <- Microbiome %>% column_to_rownames(var = "Sample") %>% filter(str_detect(Dataset, n)) %>% dplyr::select(-Dataset, -Secretor, -Amigos)
  data <- data %>% dplyr::select(where(~mean(is.na(.)) < 0.9)) %>% as.matrix()
  corr <- rcorr(data, type = "spearman")
  corr_flat <- flattenCorrMatrix(corr$r,corr$P)
  corr_flat1 <- corr_flat 
  colnames(corr_flat1) <- c("column", "row", "cor", "p") 
  corr_flat <- corr_flat1 %>% dplyr::select(row, column, cor, p) %>% rbind(corr_flat1, corr_flat)
  corr_flat <- corr_flat %>% dplyr::filter(str_detect(column, "6_SL|^DSLNH|LSTc|DFLNH|z6DDL")) %>% 
               dplyr::filter(!str_detect(row, "X6_SL|LSTc|DSLNH|DFLNH|X2_FL|X3_FL|X3_SL|DFLac|DFLNT|DSLNT|FDSLNH|FLNH|LNFP_I|LNFP_II|LNFP_III|LNH|LNnT|LNT|LSTb|z6DDL")) %>% unique()
  corr_flat$Direction <- ifelse(corr_flat$cor < 0, -1, 1)
  #corr_flat <- corr_flat %>% filter(row %in% Dir$row) %>% dplyr::select(-Direction) 
  corr_flat$padj <- p.adjust(corr_flat$p, method = "BH")
  corr_flat$Dataset <- n
  corr_hmos[[n]] <- corr_flat
}

corr_hmos <- do.call(rbind, corr_hmos)
rownames(corr_hmos) <- NULL
corr_sig <- corr_hmos %>% filter(p < 0.1) 
corr_sig$sig <- ifelse(corr_sig$p < 0.05, ifelse(corr_sig$padj < 0.05, 1, 2), 3)

rm(data, corr, corr_flat, n, corr_flat1)

corr_max <- corr_sig %>% slice_min(cor, n=1) 
corr_max1 <- corr_sig %>% slice_min(cor,n=1)
corr_max$cor <- corr_max$cor
corr_max$column <- "Extra"
corr_max$sig <- 1
corr_max1$cor <- corr_max1$cor*-1
corr_max1$column <- "Extra"
corr_max1$sig <- 3
corr_max <- rbind(corr_max, corr_max1)
corr_max1 <- corr_max
corr_max1$Dataset <- "US"
corr_max <- rbind(corr_max, corr_max1)
corcorr_max1 <- corr_max
corr_max1$Dataset <- "AU"
corr_max <- rbind(corr_max, corr_max1)
rm(corr_max1)

# Plot individual datasets
plot_hmos <- list()  # new empty list
for (n in c("AU", "ES", "US")) {
  # Filter data 
  corr_filtered <- corr_sig %>% filter(str_detect(Dataset, n)) %>% filter(!str_detect(column, "z6DDL"))
  corr_filtered <- rbind(corr_filtered, corr_max)
  # Generate nodes 
  nodes_HMO <- corr_filtered %>% dplyr::select(column) %>% unique() 
  nodes_HMO <- data.frame(node = nodes_HMO, category = rep("HMO",length(nodes_HMO)))
  names(nodes_HMO)[1] <- "node"
  nodes_bac <- corr_filtered %>% dplyr::select(row) %>% unique()
  nodes_bac <- data.frame(node = nodes_bac, category = rep("Bacteria",length(nodes_bac)))
  names(nodes_bac)[1] <- "node"
  nodes <- rbind(nodes_HMO, nodes_bac)
  node_colors <- c("#83A4C9", "#5F8AB9")
  names(node_colors) <- c("HMO", "Bacteria")
  rm(nodes_bac, nodes_HMO)
  # Generate network table
  network_table <- corr_filtered %>% dplyr::select(column, row, cor, everything())
  network_table <- network_table[complete.cases(network_table),]
  colnames(network_table)[1:3] <- c("from", "to", "correlation")
  network_table$direction <- ifelse(network_table$correlation < 0, "neg", "pos")
  network_table$weight <- abs(network_table$correlation)
  # Generate network
  network_seed <- graph_from_data_frame(d = network_table, vertices = nodes, directed = FALSE)
  SEED_Network <- ggraph(network_seed) +
    geom_edge_link(aes(color = correlation, linetype = as.factor(sig)), n=2) +
    scale_edge_color_gradientn(colors = brewer.pal(n=8, name="RdBu")) +
    geom_node_point(aes(shape = category), size = 1.5) +
    theme_void() +
    geom_node_text(aes(label = name), size = 3, repel = TRUE) +
    labs(shape = "Category", edge_color = "Spearman")
  plot_hmos[[n]] <- SEED_Network
}

patchwork::wrap_plots(plot_hmos, ncol = 3)
rm(corr_filtered, network_seed, network_table, nodes, n, node_colors, SEED_Network)

# Look at all datasets together
corr_all <- corr_sig %>% filter(str_detect(column, "z6DDL"))
corr_max <- corr_all %>% ungroup() %>% slice_max(abs(cor),n=3)
corr_max$cor <- corr_max$cor*-1
corr_max$Dataset <- "Extra"
corr_all <- rbind(corr_all, corr_max)
rm(corr_max)

# Make graph #############
# Generate nodes ###########################################
nodes_HMO <- corr_all %>% ungroup() %>% dplyr::select(Dataset) %>% unique() 
nodes_HMO <- data.frame(node = nodes_HMO, category = rep("Dataset",length(nodes_HMO)))
names(nodes_HMO)[1] <- "node"
nodes_bac <- corr_all %>% ungroup() %>% dplyr::select(row) %>% unique()
nodes_bac <- data.frame(node = nodes_bac, category = rep("Bacteria",length(nodes_bac)))
names(nodes_bac)[1] <- "node"
nodes <- rbind(nodes_HMO, nodes_bac)
node_colors <- c("#83A4C9", "#5F8AB9")
names(node_colors) <- c("Dataset", "Bacteria")
rm(nodes_bac, nodes_HMO)

# Generate network table
network_table <- corr_all %>% ungroup() %>% dplyr::select(Dataset, row, cor, everything())
network_table <- network_table[complete.cases(network_table),]
colnames(network_table)[1:3] <- c("from", "to", "correlation")
network_table$direction <- ifelse(network_table$correlation < 0, "neg", "pos")
network_table$weight <- abs(network_table$correlation)
# Generate network
network_seed <- graph_from_data_frame(d = network_table, vertices = nodes, directed = FALSE)
SEED_Network <- ggraph(network_seed) +
  geom_edge_link(aes(color = correlation, linetype = as.factor(sig))) +
  scale_edge_color_gradientn(colors = brewer.pal(n=8, name="RdBu")) +
  geom_node_point(aes(shape = category), size = 3) +
  theme_void() +
  geom_node_text(aes(label = name), size = 3, repel = TRUE) +
  labs(shape = "Category", edge_color = "Spearman")
SEED_Network
write_csv(network_table, "Microbiome_HMO_correlation_network.csv")

################### Look at alpha diversity in all three datasets #############################################################
#################################################################################################################################################
Microbiome_SEED <- read.csv("./Data/US/SEED_microbiome_shannon_genus.csv")
Microbiome_Coll <- read.csv("./Data/ES/Collado_microbiome_shannon_genus.csv")
Microbiome_Geddes <- read.csv("./Data/AU/Geddes_microbiome_shannon_genus.csv")
Microbiome_SEED$Dataset <- "US"
Microbiome_Coll$Dataset <- "ES"
Microbiome_Geddes$Dataset <- "AU"

metadata_SEED <- read.csv("./Data/US/Metadata_SEED.csv")
metadata_Coll <- read.csv("./Data/ES/Metadata_Collado.csv")
metadata_Geddes <- read.csv("./Data/AU/Metadata_Geddes.csv")
rnm <- c(Infant_age_days = "Child_1_Age", Infant_sex = "Child_1_Gender", Infant_sex = "infant_sex", C_section = "Child_1_Delivery", Preterm = "Child_1_Preterm_Birth", Parent_BMI = "Mother_BMI",
         Parent_Age = "Age", Parent_Age = "Mother_age", Ethnicity = "Race", Ethnicity = "ethnicity", Parity = "parity", Parent_BMI = "Pre_BMI")
metadata_SEED <- metadata_SEED %>% rename(any_of(rnm))
metadata_Coll <- metadata_Coll %>% rename(any_of(rnm))
metadata_Geddes <- metadata_Geddes %>% rename(any_of(rnm))

Microbiome_Coll <- left_join(Microbiome_Coll, metadata_Coll, by = "SampleID") %>% dplyr::select(-Sample) %>% rename(Sample = SampleID)
Microbiome_SEED <- left_join(Microbiome_SEED, metadata_SEED, by = "Sample")
Microbiome_Geddes <- left_join(Microbiome_Geddes, metadata_Geddes, by = "Study_ID") %>% dplyr::select(-Sample) %>% rename(Sample = Study_ID)
rm(metadata_SEED, metadata_Coll, metadata_Geddes, rnm)

Microbiome_SEED <- Microbiome_SEED %>%
  gather(Microbiome, value, 2:ncol(Microbiome_SEED)) %>%
  spread(Sample, value)
Microbiome_Coll <- Microbiome_Coll %>%
  gather(Microbiome, value, 2:ncol(Microbiome_Coll)) %>%
  spread(Sample, value)
Microbiome_Geddes <- Microbiome_Geddes %>%
  gather(Microbiome, value, 2:ncol(Microbiome_Geddes)) %>%
  spread(Sample, value)

Microbiome <- full_join(Microbiome_SEED, Microbiome_Geddes, by = "Microbiome")
Microbiome <- full_join(Microbiome, Microbiome_Coll, by = "Microbiome")

Microbiome <- Microbiome %>%
  gather(Sample, value, 2:ncol(Microbiome)) %>%
  spread(Microbiome, value)

Microbiome <- Microbiome %>%
  filter(!is.na(Secretor)) %>%
  filter(!is.na(Infant_age_days))

Microbiome$Shannon <- as.numeric(Microbiome$Shannon)
Microbiome$Parent_Age <- as.numeric(Microbiome$Parent_Age)
Microbiome$Parent_BMI <- as.numeric(Microbiome$Parent_BMI)
Microbiome$milk_production <- as.numeric(Microbiome$milk_production)
Microbiome$Breastfeeding_frequency <- as.numeric(Microbiome$Breastfeeding_frequency)
Microbiome$Infant_age_days <- as.numeric(Microbiome$Infant_age_days)

rm(Microbiome_Coll, Microbiome_Geddes, Microbiome_SEED)
Microbiome$Amigos <- str_replace(Microbiome$Amigos, "FALSE", "Low")
Microbiome$Amigos <- str_replace(Microbiome$Amigos, "TRUE", "High")
Microbiome$Amigos[is.na(Microbiome$Amigos)] <- "Medium"
Microbiome$Amigos <- factor(Microbiome$Amigos, levels = c("Low", "Medium", "High"))


# If you want to check how many points are in each group
data %>% filter(Dataset=="US") %>% count(Amigos)

data <- Microbiome

# ######## OTHER POTENTIAL CONFOUNDERS #####################################################################################
# Plot Secretor vs. Shannon Diversity ###################################
data1 <- data %>%
  filter(!is.na(Secretor))

ggplot(data1, aes(x = Secretor, y = Shannon, color = Secretor)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=0.9, shape=16, alpha = 0.5) +
  stat_compare_means(aes(label = after_stat(p.signif))) +
  facet_wrap(~Dataset, nrow = 1)

# Differences in datasets
ggplot(data, aes(x = Dataset, y = Shannon, color = Dataset)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=0.9, shape=16, alpha = 0.5) +
  stat_compare_means(aes(label = after_stat(p.signif)), comparisons = list(c("AU", "ES"), c("ES", "US"), c("AU", "US")))
#
# # Plot EBF vs. Shannon Diversity ###############
data1 <- data %>%
  filter(!is.na(EBF))

ggplot(data1, aes(x = EBF, y = Shannon, color = EBF)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=0.9, shape=16, alpha = 0.5) +
  stat_compare_means() +
  facet_wrap(~Dataset, nrow = 1)

# Plot Parent BMI vs Shannon ####################
data1 <- data %>%
  filter(!is.na(Parent_BMI))

ggscatter(data1, x = "Parent_BMI", y = "Shannon", add = "reg.line",  conf.int = TRUE,
          cor.method = "spearman", color = "Dataset") + stat_cor(aes(color = Dataset))

# Plot Parent Age vs Shannon ####################
data1 <- data %>%
  filter(!is.na(Parent_Age))

ggscatter(data1, x = "Parent_Age", y = "Shannon", add = "reg.line",  conf.int = TRUE,
          cor.method = "spearman", color = "Dataset") + stat_cor(aes(color = Dataset))

# Infant sex #######################################################
data1 <- data %>%
  filter(!is.na(Infant_sex))
ggplot(data1, aes(x = Infant_sex, y = Shannon, color = Infant_sex)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=0.9, shape=16, alpha = 0.5) +
  stat_compare_means() +
  facet_wrap(Dataset~Secretor, nrow = 1)

# Child_Age ####################################
data1 <- data %>%
  filter(!is.na(Infant_age_days)) %>%
  filter(Dataset != "AU") %>%  #Not enough difference in age to plot
  mutate(Infant_age_days = log2(Infant_age_days))

ggscatter(data1, x = "Infant_age_days", y = "Shannon", add = "reg.line",  conf.int = TRUE,
          cor.method = "spearman", color = "Dataset") + stat_cor(aes(color = Dataset))

# C-Section ###########################################################################
data1 <- data %>%
  filter(!is.na(C_section))

ggplot(data1, aes(x = C_section, y = Shannon, color = C_section)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=0.9, shape=16, alpha = 0.5) +
  stat_compare_means() +
  facet_wrap(~Dataset, nrow = 1)

# Preterm ###########################################################################
data1 <- data %>%
  filter(!is.na(Preterm))

ggplot(data1, aes(x = Preterm, y = Shannon, color = Preterm)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=0.9, shape=16, alpha = 0.5) +
  stat_compare_means() +
  facet_wrap(~Dataset, nrow = 1)

# Check whether C-section is still signifcant when accounting for preterm
data1 <- data %>%
  filter(Dataset == "US")

fit <- multinom(C_section ~ Shannon + Preterm, data=data1)
p <-  (1 - pnorm(abs(summary(fit)$coefficients/summary(fit)$standard.errors), 0, 1)) * 2
p

data1 <- data1 %>% group_by(Preterm) %>% count(C_section)
rm(fit, Microbiome, data, data1, p)
