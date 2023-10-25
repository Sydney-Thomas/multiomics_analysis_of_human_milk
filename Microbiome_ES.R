library(tidyverse)
library(phyloseq)
library(vegan)
library(Hmisc)
library(ggpubr)
library(CoDaSeq)
library(ggplot2)
library(patchwork)
library(mixOmics)

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


# Collado group provided ASV table processed with decontam v 1.14.0 (remove contamination using blanks) 
# and filtered to remove "Chloroplastos, Mitochondrias and Eukariotas/Unassigned"
# They have also provided a table with HMOs quantification + secretor status

asv_table <- read_csv("./Data/ES/otu_filtered.csv")
taxa_table <- read_csv("./Data/ES/tax_filtered.csv") %>% column_to_rownames("ASV")
taxa_table$ASV <- rownames(taxa_table)


hmos_table <- read_csv("./Data/ES/HMO_Collado.csv")
metadata <- read.csv("./Data/ES/Metadata_Collado.csv") %>% 
  dplyr::select(-Secretor) %>% 
  inner_join(hmos_table, by = "Sample")

# Check taxonomy 
any(is.na(taxa_table$Genus))

# NAs are only present at genus and species level. Fill genus level with info from family
taxa_table_fill <- taxa_table
taxa_table_fill$Genus[is.na(taxa_table_fill$Genus)] <- as.character(taxa_table_fill$Family[is.na(taxa_table_fill$Genus)])

any(is.na(taxa_table_fill$Genus))

# Prepare for phyloseq
taxa_table_final <- as.matrix(taxa_table_fill)
asv_table_final <- asv_table %>% column_to_rownames("ASV") %>% as.matrix()
metadata_final <- metadata %>% dplyr::select(-Sample) %>% filter(SampleID != "") %>% column_to_rownames("SampleID")

############################
# Generate phyloseq object #
############################
ps <- phyloseq(otu_table(asv_table_final, taxa_are_rows = TRUE), tax_table(taxa_table_final), sample_data(metadata_final))

# Remove archaea
table(tax_table(ps)[, "Kingdom"], exclude = NULL)
psF <- subset_taxa(ps, Kingdom != "Archaea")

# Remove samples with less than 10000 reads
psF_filter <- prune_samples(sample_sums(psF) > 10000, psF)

########################
# CHECK AT GENUS LEVEL #
########################
psT <- tax_glom(psF_filter, taxrank = "Genus")

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

# Visualize rarefaction curves
rarcurve <- psT %>%
  otu_table() %>%     # extract abundance table
  as("matrix") %>%    # convert phyloseq object to matrix
  t()                 # transpose - samples have to be in rows

rarecurve(rarcurve,
          step = 1000, ylab = "Unique ASVs", 
          xlab = "Sequencing Depth", label = FALSE)

# Export table of total phyla frequency for entire dataset
phylum_table <- table(tax_table(psT)[, "Phylum"], exclude = NULL) %>% as.data.frame() %>% arrange(desc(Freq))
write.csv(phylum_table, "phylum_table_Collado.csv", row.names = F)

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
write_csv(phylum_2, "Phylum_abundance_Collado.csv")

phylum_table$Var1 <- factor(phylum_table$Var1, levels = phylum_table$Var1)
phylum_pie <- phylum_table %>% ggpie("Freq", label = "Var1", fill = "Var1", legend = "none")
phylum_pie

rm(ps, psF)

################################
# ALPHA DIVERSITY at ASV level #
################################
ps_asv_rar <- rarefy_even_depth(psF_filter, 
                                sample.size = min(sample_sums(psF_filter)), 
                                rngseed = 1234, # set this so that subsampling is reproducible
                                replace = FALSE, 
                                trimOTUs = TRUE)


# Calculate Shannon diversity at asv level 
shannon_asv <- ps_asv_rar %>%
  estimate_richness(measures = "Shannon") %>% rownames_to_column("SampleID")

write.csv(shannon_asv, "Collado_microbiome_shannon_asv.csv", row.names = F)
i <- "EBF"

shannon_asv %>% left_join(metadata, by = "SampleID") %>%    
  dplyr::filter(Secretor == 1) %>%
  dplyr::filter(!is.na(!!sym(i))) %>% 
  ggboxplot(x = i, y = "Shannon", add = "jitter",
            add.params = list(alpha = 0.3, color = i), legend = "right",
            xlab = FALSE, ylab = FALSE, nrow = 1, title = "Shannon Diversity - ASV") +
  stat_compare_means() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

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
  estimate_richness(measures = "Shannon") %>% rownames_to_column("SampleID")

write.csv(shannon_asv, "Collado_microbiome_shannon_genus.csv", row.names = F)

i <- "C_section"

shannon_genus %>% left_join(metadata, by = "SampleID") %>%    
  dplyr::filter(Secretor == 1) %>%
  dplyr::filter(!is.na(!!sym(i))) %>% 
  ggboxplot(x = i, y = "Shannon", add = "jitter",
            add.params = list(alpha = 0.3, color = i), legend = "right",
            xlab = FALSE, ylab = FALSE, nrow = 1, title = "Shannon Diversity - Genus") +
  stat_compare_means() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

rm(ps_asv_rar)

# Filter out genera not detected in at least 5 samples
prev <- apply(X = otu_table(psT),
              MARGIN = ifelse(taxa_are_rows(psT), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf <- data.frame(Prevalence = prev,
                     TotalAbundance = taxa_sums(psT),
                     rownames(otu_table(psT)))
psA <- prune_taxa((prev > 5), psT)

# Extract count table and metadata
OTUG <- as.data.frame(psA@otu_table) %>% t() %>% as.data.frame()
metaG <- psA@sam_data %>% as.matrix() %>% as.data.frame() %>% rownames_to_column("SampleID")

# CLR Transformation
OTUzeros <- cmultRepl(OTUG, label = 0, method = "CZM", output = "p-counts")
OTUclr <- codaSeq.clr(OTUzeros, samples.by.row = TRUE) %>% as.data.frame()

# Robust CLR transformation
OTU_rclr <- OTUG %>% decostand(method = "rclr")

# Add taxonomy information
taxa_clr_filter <- taxa_table_fill %>% as.data.frame() %>% dplyr::filter(ASV %in% colnames(OTUclr)) %>% dplyr::select(Genus)

all(taxa_clr_filter$rowname == colnames(OTUclr))
all(taxa_clr_filter$rowname == colnames(OTUG))
all(taxa_clr_filter$rowname == colnames(OTU_rclr))

colnames(OTUG) <- taxa_clr_filter$Genus
colnames(OTUclr) <- taxa_clr_filter$Genus
colnames(OTU_rclr) <- taxa_clr_filter$Genus

OTUG <- rownames_to_column(OTUG, var = "SampleID")
OTUclr <- rownames_to_column(OTUclr, var = "SampleID")
OTU_rclr <- rownames_to_column(OTU_rclr, var = "SampleID")

write_csv(OTUG, "Microbiome_table_Collado_genus.csv")
write_csv(OTUclr, "Microbiome_table_Collado_genus_clr.csv")
write_csv(OTU_rclr, "Microbiome_table_Collado_genus_rclr.csv")


##########################################################
# Remove samples for which there is not HMOs information #
##########################################################

######## Choose whether to look at secretors or non-secretors! #######################################
metadata_filter <- metaG %>% dplyr::filter(!is.na(Secretor)) %>% as.data.frame() %>% mutate(across(c("Secretor":"Fuc"), as.numeric)) %>% dplyr::filter(Secretor == 1)

OTU_filter <- OTU_rclr %>% dplyr::filter(SampleID %in% metadata_filter$SampleID) %>% column_to_rownames("SampleID") %>% as.data.frame()

# PCA
PCA_all <- mixOmics::pca(OTU_filter, ncomp = 2, center = TRUE, scale = FALSE)
PCA_all_scores <- data.frame(PCA_all$variates$X) %>% rownames_to_column("SampleID") %>% dplyr::left_join(metadata_filter)

i <- "C_section"

PCA_plot <- PCA_all_scores %>% dplyr::filter(Secretor == 1) %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6, 
            title = paste("PCA - Microbiome", i, sep = " "),
            xlab = paste("PC1 (", round(PCA_all$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_all$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_all_scores %>% dplyr::filter(Secretor == 1) %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()
PCA_plot

dist_genus <- vegdist(OTU_filter, method = "euclidean")
permanova_genus <- adonis2(dist_genus ~ PCA_all_scores[[i]], PCA_all_scores, na.action = na.omit, method = "euclidean")


# Plot HMOs abundance
pca_hmos <- list()  # new empty list

for (n in c("X6_SL", "LSTc", "DSLNH", "DFLNH")) {
  
  PCA_plot <- PCA_all_scores %>% dplyr::filter(Secretor == 1) %>%
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

microbe_hmos <- patchwork::wrap_plots(pca_hmos, ncol = 2)
microbe_hmos

############
# PCA hmos #
############

HMOs_log <- metadata_filter %>% dplyr::select(SampleID, X2_FL:DSLNH) %>% mutate(across(X2_FL:DSLNH, log2)) %>% column_to_rownames(var = "SampleID")

# PCA
PCA_all <- mixOmics::pca(HMOs_log, ncomp = 2, center = TRUE, scale = FALSE)
PCA_all_scores <- data.frame(PCA_all$variates$X) %>% rownames_to_column("SampleID") %>% dplyr::left_join(metadata_filter)

PCA_plot <- PCA_all_scores %>% dplyr::filter(Secretor == 1) %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6, 
            title = paste("PCA - HMOs", i, sep = " "),
            xlab = paste("PC1 (", round(PCA_all$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_all$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_all_scores %>% dplyr::filter(Secretor == 1) %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()
PCA_plot

dist_hmos <- vegdist(HMOs_log, method = "euclidean")
permanova_hmos <- adonis2(dist_hmos ~ PCA_all_scores$[[i]], PCA_all_scores, na.action = na.omit, method = "euclidean")

mixOmics::plotVar(PCA_all, cex = 3)
