#####################################################################
########## Bacterial isolate analysis for Thomas et al. 2024
#####################################################################

## Load libraries
library(ggplot2)
library(tidyr)
library(growthcurver)
library(ggpubr)
library(rstatix)
library(eoffice)
library(plyr)

## Set working director
setwd("D:/")

## Import data
df <- read.csv("Growth curve data.csv")
hmo <- read.csv("HMO degradation data.csv")

###################################
####### Generate growth curves
###################################

df$Treatment <- factor(df$Treatment, levels = c("No treatment", "pHMOs", "6'SL", "LSTc"))

ggplot() + 
  geom_smooth(df,
              mapping = aes(Time, OD, colour = Treatment, group = Treatment), 
              size =1) +
  facet_wrap(~Strain, ncol = 4) +
  theme_bw() +
  xlab("Time (h)") + ylab("Optical Density") +
  scale_x_continuous(breaks = c(0,6,12,18,24), labels = c("0","6","12","18","24"))

topptx(file = "Rothia.pptx", append = TRUE, width = 8, height = 2, units = "cm")

###################################
####### Calculate growth curve metrics
###################################

df$assay2 <- paste(df$Assay, df$Well, sep = "_")
df2 <- subset(df, select = c(assay2, Time, OD))
df2 <- spread(df2, assay2, OD)
df2 <- SummarizeGrowthByPlate(df2)

df3 <- subset(df2, select = c(sample, auc_e, k, r))

df3$Strain <- df$Strain[match(df3$sample, df$assay2)]
df3$Treatment <- df$Treatment[match(df3$sample, df$assay2)]

df4 <- gather(df3, Metric, Value, -sample, -Strain, -Treatment)

## Perform statistical comparisons
tukey <- data.frame(matrix(ncol = 7, nrow = 0, dimnames = list(NULL, c("Strain", "Metric", "Comparison", "Treatment.diff","Treatment.lwr","Treatment.upr","Treatment.p.adj"))))

strains <- unique(df4$Strain)
metrics <- unique(df4$Metric)

for (i in 1:length(strains)) {
  df5 <- subset(df4, df4$Strain == strains[i])
  for (j in 1:length(metrics)) {
    df6 <- subset(df5, df5$Metric == metrics[j])
    res.aov <- aov(Value ~ Treatment, data = df6)
    tukey1 <- (TukeyHSD(res.aov, which = "Treatment"))
    tukey2 <- as.data.frame(tukey1[1:1])
    tukey2$Metric <- metrics[j]
    tukey2$Strain <- strains[i]
    tukey2 <- tibble::rownames_to_column(tukey2, "Comparison")
    tukey <- rbind(tukey, tukey2)
    }
}
  
tukey$BH <- p.adjust(tukey$Treatment.p.adj, method = "BH")

tukey <- add_significance(tukey,
  p.col = "BH", output.col = NULL,
  cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
  symbols = c("****", "***", "**", "*", "ns")
)

###################################
####### Generate boxplots
###################################

ggplot(df4, aes(Treatment, Value, fill = Treatment)) +
  geom_boxplot() +
  geom_jitter(size = 1) +
  facet_grid(Metric ~ Strain, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_blank()) #+
  #scale_y_continuous(limits = c(0, 1.9))

topptx(file = "Rothia.pptx", append = TRUE, width = 7, height = 7, units = "cm")


###################################
####### HMO degradation 
###################################

## Statistically compare treatments
t <- gather(hmo, HMO, Concentration, -Sample, -Assay, -Strain, -Treatment, -Replicate)
t <- t[!(is.na(t$Concentration)), ]

assays <- c(240228, 240305)
strains <- c("Rothia 1", "Rothia 2")

t$TreatHMO <- paste(t$Treatment, t$HMO, sep = "_")
treatments <- unique(t$TreatHMO)

ttests <-  data.frame(matrix(ncol = 11, nrow = 0, dimnames = list(NULL, c(".y.", "group1", "group2", "n1", "n2", "statistic", "df", "p", "Strain", "Treatment", "HMO"))))

## Stats for Rothia 1 and Rothia 2
t2 <- subset(t, t$Assay == 240219)
  
  for (j in 1:length(strains)) {
    t3 <- subset(t2, t2$Strain %in% c(strains[j], "Negative control"))
    
    for (k in 1: length(treatments)) {
      t4 <- subset(t3, t3$TreatHMO == treatments[k])
      t5 <- t4 %>% t_test(Concentration ~ Strain)
      
      t5$Treatment <- unique(t4$Treatment)
      t5$HMO <- unique(t4$HMO)
      
      ttests <- rbind(t5, ttests)
    }
  }

## Stats for Rothia 3 and Veillonella 1
for (j in 1:length(assays)) {
  t2 <- subset(t, t$Assay == assays[j])
    
    for (k in 1: length(treatments)) {
      t4 <- subset(t2, t2$TreatHMO == treatments[k])
      t5 <- t4 %>% t_test(Concentration ~ Strain)
      
      t5$Treatment <- unique(t4$Treatment)
      t5$HMO <- unique(t4$HMO)
      
      ttests <- rbind(t5, ttests)
      
    }
  }

ttests <- adjust_pvalue(ttests, p.col = "p", output.col = "p.adj", method = "BH")
ttests <- add_significance(ttests, p.col = "p.adj", output.col = "p.adj.signif",
                           cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
                           symbols = c("****", "***", "**", "*", "ns"))
