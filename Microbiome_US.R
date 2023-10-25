# Load packages
library(tidyverse)
library(lubridate)
library(phyloseq)
library(CoDaSeq)
library(viridis)
library(ggpubr)
library(mixOmics)
library(vegan)
library(viridis)
library(effsize)
library(gridExtra)
library(patchwork)
library(janitor)
library(Hmisc)

# Function to perform pre-filtering (from mixOmics website) - remove OTUs accounting for less 1% of abundance
low_count_removal <- function(data, percent = 0.01) {
  
  keep_otu <- which(colSums(data)*100/(sum(colSums(data))) > percent)
  data_filter <- data[,keep_otu]
  return(list(data_filter = data_filter, keep_otu = keep_otu))
  
}

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

# Load data
metadata <- read.csv("./Data/US/Metadata_SEED.csv", header = TRUE)
microbiome <- read.csv("./Data/US/ASV_table.csv", header = TRUE)
colnames(microbiome) <- str_replace_all(colnames(microbiome), "X", "")
hmos <- read.csv("./Data/US/HMO_SEED.csv")
metadata_hmos <- metadata %>% dplyr::select(-Secretor) %>% dplyr::left_join(hmos, by = "Sample")

# Taxonomy table
taxonomy_table <- microbiome %>% dplyr::select(ASV, taxonomy)
ASV_sequence <- microbiome %>% dplyr::select(ASV, OTU_ID)

taxonomy_table_split <- taxonomy_table %>% separate(taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")

taxonomy_table_split$Kingdom <- gsub("d__", "", taxonomy_table_split$Kingdom)
taxonomy_table_split$Phylum <- gsub("p__", "", taxonomy_table_split$Phylum)
taxonomy_table_split$Class <- gsub("c__", "", taxonomy_table_split$Class)
taxonomy_table_split$Order <- gsub("o__", "", taxonomy_table_split$Order)
taxonomy_table_split$Family <- gsub("f__", "", taxonomy_table_split$Family)
taxonomy_table_split$Genus <- gsub("g__", "", taxonomy_table_split$Genus)
taxonomy_table_split$Species <- gsub("s__", "", taxonomy_table_split$Species)

taxonomy_table_split[taxonomy_table_split == "" | taxonomy_table_split == " "] <- NA

# Fill unknown taxa
taxonomy_table_fill <- taxonomy_table_split
taxonomy_table_fill$Class[is.na(taxonomy_table_fill$Class)] <- as.character(taxonomy_table_fill$Phylum[is.na(taxonomy_table_fill$Class)])
taxonomy_table_fill$Order[is.na(taxonomy_table_fill$Order)] <- as.character(taxonomy_table_fill$Class[is.na(taxonomy_table_fill$Order)])
taxonomy_table_fill$Family[is.na(taxonomy_table_fill$Family)] <- as.character(taxonomy_table_fill$Order[is.na(taxonomy_table_fill$Family)])
taxonomy_table_fill$Genus[is.na(taxonomy_table_fill$Genus)] <- as.character(taxonomy_table_fill$Family[is.na(taxonomy_table_fill$Genus)])

# Prepare for phyloseq
taxonomy_table_fill <- taxonomy_table_fill %>% column_to_rownames("ASV")
taxa_final <- as.matrix(taxonomy_table_fill)

# OTU table
OTUs <- microbiome %>% dplyr::select(-c("OTU_ID", "taxonomy")) %>% column_to_rownames("ASV") 
OTUs <- OTUs[colnames(OTUs) %in% metadata_hmos$Sample]
OTUs_sort <- OTUs %>% dplyr::select(order(colnames(OTUs), decreasing = TRUE))
OTU_final <- as.matrix(OTUs_sort)
metadata_final <- metadata_hmos %>% column_to_rownames("Sample")

rm(OTUs_sort, taxonomy_table_split, hmos, metadata)

############################
# Generate phyloseq object #
############################
ps <- phyloseq(otu_table(OTU_final, taxa_are_rows = TRUE), sample_data(metadata_final), tax_table(taxa_final))

ps <- filter_taxa(ps, function(x) sum(x)>10, TRUE)

# Remove archaea and eukaryota
table(tax_table(ps)[, "Kingdom"], exclude = NULL)
psF <- subset_taxa(ps, Kingdom != "Archaea")
psF <- subset_taxa(psF, Kingdom != "Eukaryota")
psF <- subset_taxa(psF, Kingdom != "Unassigned")

# Check phyla and remove unclassified Phyla
table(tax_table(psF)[, "Phylum"], exclude = NULL)
psF <- subset_taxa(psF, Phylum != "NA")

# Get rid of DNA from chloroplasts or mitochondria
table(tax_table(psF)[, "Family"], exclude = NULL)
psF <- subset_taxa(psF, Family != " Chloroplast")
psF <- subset_taxa(psF, Family != " mitochondria")

# Remove samples with less than 10000 reads
ps_filter <- prune_samples(sample_sums(psF) > 1000, psF)

########################
# CHECK AT GENUS LEVEL #
########################
psT <- tax_glom(ps_filter, taxrank = "Genus")

# Generate density plot of sequencing depth
dplot <- ggdensity(data.frame(Sequencing_Depth = sample_sums(psT)),
                   x = "Sequencing_Depth",
                   xlab = "Sequencing Depth",
                   y = "..count..", alpha = 0.6,
                   ylab = "Number of samples (%)",title = "Density Plot - Microbiome",
                   legend = "right") +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 5), legend.title = element_text(size = 8),
        legend.text = element_text(size = 8))
dplot

# # Visualize rarefaction curves
rarcurve <- psT %>%
  otu_table() %>%     # extract abundance table
  as("matrix") %>%    # convert phyloseq object to matrix
  t()                 # transpose - samples have to be in rows

rarecurve(rarcurve,
          step = 100, ylab = "Unique Genera",
          xlab = "Sequencing Depth", label = FALSE)

phylum_table <- table(tax_table(psT)[, "Phylum"], exclude = NULL) %>% as.data.frame() %>% arrange(desc(Freq))
write.csv(phylum_table, "phylum_table_SEED.csv", row.names = F)

# Export table of phylum abundance by sample
phylum_2 <- psT %>% tax_glom(taxrank = "Phylum") %>% psmelt() %>%
  dplyr::select(Phylum, Sample, Abundance) %>% spread(Sample, Abundance)
# Filter out phyla that weren't found in at least 5 samples
phylum_2[phylum_2 == 0] <- NA
phylum_2 <- phylum_2 %>% mutate(sum_na = rowSums(is.na(.)))
phylum_2 <- phylum_2 %>% filter(sum_na<(ncol(phylum_2)-6)) %>% dplyr::select(-sum_na)
phylum_2[is.na(phylum_2)] <- 0
# Rclr transform
phylum_2 <- phylum_2 %>%
  gather(Sample, value, 2:ncol(phylum_2)) %>%
  spread(Phylum, value)
phylum_2 <- column_to_rownames(phylum_2, var = "Sample")
phylum_2 <- phylum_2 %>% decostand(method = "rclr")
phylum_2 <- phylum_2 %>% rownames_to_column(var = "Sample")
write_csv(phylum_2, "Phylum_abundance_SEED.csv")

phylum_table$Var1 <- factor(phylum_table$Var1, levels = phylum_table$Var1)
phylum_pie <- phylum_table %>% ggpie("Freq", label = "Var1", fill = "Var1", legend = "none")
phylum_pie

rm(ps, psF)

################################################
# ALPHA DIVERSITY ## ASV LEVEL ####################
################################################
ps_asv_rar <- rarefy_even_depth(ps_filter, 
                                sample.size = min(sample_sums(ps_filter)), 
                                rngseed = 1234, # set this so that subsampling is reproducible
                                replace = FALSE, 
                                trimOTUs = TRUE)


# Calculate Shannon diversity at asv level 
shannon_asv <- ps_asv_rar %>%
  estimate_richness(measures = "Shannon") %>% rownames_to_column("Sample")

write_csv(shannon_asv, "SEED_microbiome_shannon_asv.csv")

i <- "Child_1_Preterm_Birth"

shannon_asv %>% left_join(metadata_hmos, by = "Sample") %>%    
  dplyr::filter(Secretor == 1) %>%
  dplyr::filter(Child_1_Age < 2000) %>% 
  dplyr::filter(!is.na(!!sym(i))) %>% 
  ggboxplot(x = i, y = "Shannon", add = "jitter",
            add.params = list(color = "Child_1_Age"), legend = "right",
            xlab = FALSE, ylab = FALSE, nrow = 1, title = "Shannon Diversity - ASV") +
  scale_color_distiller(palette = "RdBu") +
  stat_compare_means() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

shannon_asv %>% left_join(metadata_hmos, by = "Sample") %>%    
  dplyr::filter(Secretor == 1) %>%
  dplyr::filter(Child_1_Age < 2000) %>% 
  dplyr::filter(!is.na(Child_1_Age)) %>% 
  ggplot(aes(x = Child_1_Age, y = Shannon)) +
  geom_point() +
  geom_smooth() 


##################################
# ALPHA DIVERSITY at Genus level #
##################################
ps_asv_rar <- rarefy_even_depth(psT, 
                                sample.size = min(sample_sums(psT)), 
                                rngseed = 1234, # set this so that subsampling is reproducible
                                replace = FALSE, 
                                trimOTUs = TRUE)


# Calculate Shannon diversity at genus level 
shannon_genus <- ps_asv_rar %>%
  estimate_richness(measures = "Shannon") %>% rownames_to_column("Sample")

write_csv(shannon_genus, "SEED_microbiome_shannon_genus.csv")

i <- "Child_1_Delivery"

shannon_genus %>% left_join(metadata_hmos, by = "Sample") %>%    
  dplyr::filter(Secretor == 1) %>%
  dplyr::filter(!is.na(!!sym(i))) %>% 
  ggboxplot(x = i, y = "Shannon", add = "jitter",
            add.params = list(alpha = 0.3, color = i), legend = "right",
            xlab = FALSE, ylab = FALSE, nrow = 1, title = "Shannon Diversity - Genus") +
  stat_compare_means() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

rm(ps_asv_rar)

##############
# ORDINATION #
##############

# Filter out genera not detected in at least 3% of the samples (> 12 samples)
prev <- apply(X = otu_table(psT),
              MARGIN = ifelse(taxa_are_rows(psT), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf <- data.frame(Prevalence = prev,
                     TotalAbundance = taxa_sums(psT),
                     rownames(otu_table(psT)))
prevalenceThreshold <- 0.03 * nsamples(psT)
psA <- prune_taxa((prev > prevalenceThreshold), psT)

# Extract count table and metadata
OTUG <- as.data.frame(psA@otu_table) %>% t() %>% as.data.frame()
metaG <- psA@sam_data %>% as.matrix() %>% as.data.frame() %>% rownames_to_column("Sample")

# Remove samples with less of 1000 reads
c <- data.frame(Seq = rowSums(OTUG)) %>% filter(Seq > 1000)
psS <- prune_samples(rownames(c), psA)

# Extract count table and metadata
OTUG <- as.data.frame(psS@otu_table) %>% t() %>% as.data.frame()
metaG <- psS@sam_data %>% as.matrix() %>% as.data.frame() %>% rownames_to_column("Sample") 

# CLR Transformation
OTUzeros <- cmultRepl(OTUG, label = 0, method = "CZM", output = "p-counts")
OTUclr <- codaSeq.clr(OTUzeros, samples.by.row = TRUE) %>% as.data.frame()

# Robust CLR transformation
OTU_rclr <- OTUG %>% decostand(method = "rclr")

# Add taxonomy information
taxa_clr_filter <- taxa_final %>% as.data.frame() %>% rownames_to_column("ASV") %>% dplyr::filter(ASV %in% colnames(OTUclr)) %>% dplyr::select(Genus)

all(taxa_clr_filter$rowname == colnames(OTUclr))
all(taxa_clr_filter$rowname == colnames(OTUG))
all(taxa_clr_filter$rowname == colnames(OTU_rclr))

colnames(OTUG) <- taxa_clr_filter$Genus
colnames(OTUclr) <- taxa_clr_filter$Genus
colnames(OTU_rclr) <- taxa_clr_filter$Genus

OTUG <- rownames_to_column(OTUG, var = "Sample")
OTUclr <- rownames_to_column(OTUclr, var = "Sample")
OTU_rclr <- rownames_to_column(OTU_rclr, var = "Sample")

write_csv(OTUG, "Microbiome_table_SEED_genus.csv")
write_csv(OTUclr, "Microbiome_table_SEED_genus_clr.csv")
write_csv(OTU_rclr, "Microbiome_table_SEED_genus_rclr.csv")

##########################################################
# Remove samples for which there is not HMOs information #
##########################################################

######## Choose whether to look at secretors or non-secretors! #######################################
metadata_filter <- metaG %>% dplyr::filter(!is.na(Secretor)) %>% as.data.frame() %>% mutate(across(c("Secretor":"Fuc"), as.numeric)) %>% dplyr::filter(Secretor == 1)

OTU_filter <- OTU_rclr %>% dplyr::filter(Sample %in% metadata_filter$Sample) %>% column_to_rownames("Sample") %>% as.data.frame()

# PCA
PCA_all <- mixOmics::pca(OTU_filter, ncomp = 2, center = TRUE, scale = FALSE)
PCA_all_scores <- data.frame(PCA_all$variates$X) %>% rownames_to_column("Sample") %>% dplyr::left_join(metadata_filter)

i <- "Amigos"

PCA_plot <- PCA_all_scores %>% 
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6, 
            title = paste("PCA - Microbiome", i, sep = " "),
            xlab = paste("PC1 (", round(PCA_all$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_all$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_all_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()
PCA_plot

dist_genus <- vegdist(OTU_filter, method = "euclidean")
permanova_genus <- adonis2(dist_genus ~ PCA_all_scores[[i]], PCA_all_scores, na.action = na.omit, method = "euclidean")


# Plot HMOs abundance
pca_hmos <- list()  # new empty list

for (n in c("X6_SL", "LSTc", "DSLNH", "DFLNH")) {
  
  PCA_plot <- PCA_all_scores %>% 
    mutate_at(c("X6_SL", "LSTc", "DSLNH", "DFLNH"), log2) %>%
    ggscatter(x = "PC1", y = "PC2", color = n, alpha = 0.6, 
              title = paste("PCA - Microbiome", n, sep = " "),
              xlab = paste("PC1 (", round(PCA_all$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_all$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) + scale_color_viridis_c() +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  pca_hmos[[n]] <- PCA_plot  # add each plot into plot list
  
}

microbe_hmos <- wrap_plots(pca_hmos, ncol = 2)
microbe_hmos

############
# PCA hmos #
############

HMOs_log <- metadata_filter %>% dplyr::select(Sample, X2_FL:DSLNH) %>% mutate(across(X2_FL:DSLNH, log2)) %>% column_to_rownames(var = "Sample")

# PCA
PCA_all <- mixOmics::pca(HMOs_log, ncomp = 2, center = TRUE, scale = FALSE)
PCA_all_scores <- data.frame(PCA_all$variates$X) %>% rownames_to_column("Sample") %>% dplyr::left_join(metadata_filter)

PCA_plot <- PCA_all_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6, 
            title = paste("PCA - HMOs", i, sep = " "),
            xlab = paste("PC1 (", round(PCA_all$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_all$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_all_scores %>% group_by((!!sym(i))) %>% 
             summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()
PCA_plot

dist_hmos <- vegdist(HMOs_log, method = "euclidean")
permanova_hmos <- adonis2(dist_hmos ~ PCA_all_scores[[i]], PCA_all_scores, na.action = na.omit, method = "euclidean")

mixOmics::plotVar(PCA_all, cex = 3)


