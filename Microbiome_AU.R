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


library(tidyverse)
library(phyloseq)
library(vegan)
library(Hmisc)
library(ggpubr)
library(CoDaSeq)
library(ggplot2)
library(patchwork)
library(mixOmics)

data <- read.delim("./Data/AU/1humanmilk.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons (1).shared") %>% 
  dplyr::select(-c("label", "numOtus")) 
taxa <- read.delim("./Data/AU/1humanmilk.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons (1).taxonomy") %>% 
  dplyr::select(-Size) %>% column_to_rownames("OTU")
taxa$Taxonomy <- gsub("\\([^()]*\\)", "", taxa$Taxonomy)
hmos <- read_csv("./Data/AU/HMO_Geddes.csv") %>% dplyr::select(-Infant_age_days)
metadata <- read_csv("./Data/AU/Metadata_Geddes.csv")

# Geddes data is PACBIO. They have generated OTUs and annotated them with SILVA.
# Even looking at ASVs will not be exactly the same as the short reads Illumina.
# Work with OTUs and collapse them at Genus level to make them comparable

# Fix taxa table to generate phyloseq object
taxa_split <- taxa %>% separate(Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";")

# Prepare for phyloseq
otu_table <- data %>% column_to_rownames("Group") %>% t()
taxa_table <- taxa_split %>% as.matrix()

metadata_filter <- metadata %>% dplyr::filter(Study_ID %in% data$Group) %>% left_join(hmos) %>% column_to_rownames("Study_ID")

############################
# Generate phyloseq object #
############################
ps <- phyloseq(otu_table(otu_table, taxa_are_rows = TRUE), tax_table(taxa_table), sample_data(metadata_filter))

# Remove archaea
table(tax_table(ps)[, "Class"], exclude = NULL)
psF <- subset_taxa(ps, Kingdom != "Archaea")

# Remove samples with less than 10000 reads
ps_filter <- prune_samples(sample_sums(psF) > 10000, psF)

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

# Visualize rarefaction curves
rarcurve <- psT %>%
  otu_table() %>%     # extract abundance table
  as("matrix") %>%    # convert phyloseq object to matrix
  t()                 # transpose - samples have to be in rows

rarecurve(rarcurve,
          step = 100, ylab = "Unique Genera", 
          xlab = "Sequencing Depth", label = FALSE)

phylum_table <- table(tax_table(psT)[, "Phylum"], exclude = NULL) %>% as.data.frame() %>% arrange(desc(Freq))
write.csv(phylum_table, "phylum_table_Geddes.csv", row.names = F)

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
write_csv(phylum_2, "Phylum_abundance_Geddes.csv")

phylum_table$Var1 <- factor(phylum_table$Var1, levels = phylum_table$Var1)
phylum_pie <- phylum_table %>% ggpie("Freq", label = "Var1", fill = "Var1", legend = "none")
phylum_pie

##################################
# ALPHA DIVERSITY at Genus level #
##################################
ps_asv_rar <- rarefy_even_depth(psT, 
                                sample.size = min(sample_sums(psT)), 
                                rngseed = 1234, # set this so that subsampling is reproducible
                                replace = FALSE, 
                                trimOTUs = TRUE)

# Calculate Shannon diversity at asv level 
shannon_genus <- ps_asv_rar %>%
  estimate_richness(measures = "Shannon") %>% rownames_to_column("Study_ID")

write.csv(shannon_genus, "Geddes_microbiome_shannon_genus.csv", row.names = F)

i <- "C_section"

shannon_genus %>% left_join(metadata, by = "Study_ID") %>%    
  dplyr::filter(Secretor == 1) %>%
  dplyr::filter(!is.na(!!sym(i))) %>% 
  ggboxplot(x = i, y = "Shannon", add = "jitter",
            add.params = list(alpha = 0.3, color = i), legend = "right",
            xlab = FALSE, ylab = FALSE, nrow = 1, title = "Shannon Diversity - Genus") +
  stat_compare_means() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())



# Filter out genera not detected in at least 3 samples
prev <- apply(X = otu_table(psT),
              MARGIN = ifelse(taxa_are_rows(psT), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf <- data.frame(Prevalence = prev,
                     TotalAbundance = taxa_sums(psT),
                     rownames(otu_table(psT)))
psA <- prune_taxa((prev > 2), psT)

# Extract count table and metadata
OTUG <- as.data.frame(psA@otu_table) %>% t() %>% as.data.frame()
metaG <- psA@sam_data %>% as.matrix() %>% as.data.frame() %>% rownames_to_column("SampleID")

# CLR Transformation
OTUzeros <- cmultRepl(OTUG, label = 0, method = "CZM", output = "p-counts")
OTUclr <- codaSeq.clr(OTUzeros, samples.by.row = TRUE) %>% as.data.frame()

# Robust CLR transformation
OTU_rclr <- OTUG %>% decostand(method = "rclr")

# Add taxonomy information
taxa_clr_filter <- taxa_table %>% as.data.frame() %>% rownames_to_column("OTU") %>%
  dplyr::filter(OTU %in% colnames(OTUclr)) %>% dplyr::select(Genus)

all(taxa_clr_filter$rowname == colnames(OTUclr))
all(taxa_clr_filter$rowname == colnames(OTUG))
all(taxa_clr_filter$rowname == colnames(OTU_rclr))

colnames(OTUG) <- taxa_clr_filter$Genus
colnames(OTUclr) <- taxa_clr_filter$Genus
colnames(OTU_rclr) <- taxa_clr_filter$Genus

OTUG <- rownames_to_column(OTUG, var = "Study_ID")
OTUclr <- rownames_to_column(OTUclr, var = "Study_ID")
OTU_rclr <- rownames_to_column(OTU_rclr, var = "Study_ID")
write_csv(OTUG, "Microbiome_table_Geddes_genus.csv")
write_csv(OTUclr, "Microbiome_table_Geddes_genus_clr.csv")
write_csv(OTU_rclr, "Microbiome_table_Geddes_genus_rclr.csv")

##########################################################
# Remove samples for which there is not HMOs information #
##########################################################

OTU_filter <- OTUclr %>% rownames_to_column("SampleID") %>% 
  dplyr::filter(SampleID %in% (metaG %>% dplyr::filter(!is.na(Secretor)))$SampleID) %>% 
  column_to_rownames("SampleID")

metadata_filter <- metaG %>% dplyr::filter(!is.na(Secretor))

# PCA
PCA_all <- mixOmics::pca(OTU_filter, ncomp = 2, center = TRUE, scale = FALSE)
PCA_all_scores <- data.frame(PCA_all$variates$X) %>% rownames_to_column("SampleID") %>% dplyr::left_join(metadata_filter)

i <- "Secretor"

PCA_plot <- PCA_all_scores %>% dplyr::filter(!(is.na(Secretor))) %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6, 
            title = paste("PCA - Microbiome", i, sep = " "),
            xlab = paste("PC1 (", round(PCA_all$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_all$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_all_scores %>% dplyr::filter(!(is.na(Secretor))) %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

dist_genus <- vegdist(OTU_filter, method = "euclidean")
permanova_genus <- adonis2(dist_genus ~ PCA_all_scores$Secretor, PCA_all_scores, na.action = na.omit, method = "euclidean")

# Plot HMOs abundance
pca_hmos <- list()  # new empty list

for (i in c("X6_SL", "LSTc", "DSLNH", "DFLNH")) {
  
  PCA_plot <- PCA_all_scores %>% 
    mutate_all(as.numeric) %>%
    mutate_at(c("X6_SL", "LSTc", "DSLNH", "DFLNH"), log10) %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6, 
              title = paste("PCA - Microbiome", i, sep = " "),
              xlab = paste("PC1 (", round(PCA_all$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_all$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) + scale_color_viridis_c() +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  pca_hmos[[i]] <- PCA_plot  # add each plot into plot list
  
}

microbe_hmos <- wrap_plots(pca_hmos, ncol = 2)

############
# PCA hmos #
############

HMOs_log <- metadata_filter %>% dplyr::select("SampleID", 27:48) %>% column_to_rownames("SampleID") %>%
  mutate_all(as.numeric) %>% mutate_all(log10) %>% dplyr::select(-c("SUM", "Sia", "Fuc"))

# PCA
PCA_all <- mixOmics::pca(HMOs_log, ncomp = 2, center = TRUE, scale = FALSE)
PCA_all_scores <- data.frame(PCA_all$variates$X) %>% rownames_to_column("SampleID") %>% dplyr::left_join(metadata_filter)

i <- "Secretor"

PCA_plot <- PCA_all_scores %>% dplyr::filter(!(is.na(Secretor))) %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6, 
            title = paste("PCA - HMOs", i, sep = " "), label = "SampleID",
            xlab = paste("PC1 (", round(PCA_all$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_all$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_all_scores %>% dplyr::filter(!(is.na(Secretor))) %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

dist_hmos <- vegdist(HMOs_log, method = "euclidean")
permanova_hmos <- adonis2(dist_hmos ~ PCA_all_scores$Secretor, PCA_all_scores, na.action = na.omit, method = "euclidean")

mixOmics::plotVar(PCA_all, cex = 3)
