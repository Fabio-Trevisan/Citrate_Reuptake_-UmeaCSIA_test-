library(readr)
library(ggpubr)
library(tidyverse)
library(plyr)
library(dplyr)
library(rstatix)
library(purrr)
library(agricolae)
library(reshape2)


#CHECK CORRELATION BETWEEN [C] & ENRICHMENT% AND RELATIVE ABUNDANCE & ENRICHMENT%


#OA
#Read CSV ####
table <- read.csv("20230711 3-NPH acids from in-house script_IsoCor_res.csv", sep=";", header=T)
Class <- "OA"

#AA
#Read CSV ####
table <- read.csv("20230714_AccQ-Tag_AA_peak areas_inhouse_script 2.0_IsoCor_res.csv", sep=";", header=T)
Class <- "AA"

#Dataset preparation ####
#cleaning
df <- table[,-c(3,5,7:9)] #remove un-useful columns
df <- df[df$isotopologue<=0,] #filter only isotopologue M+0
df <- df[!str_detect(df$sample, "blank"),] #remove blank samples
df <- df[!str_detect(df$sample, "std"),] #remove standards
df <- df[!str_detect(df$sample, "QC"),] #remove QCs
row.names(df) <- NULL
df <- df[,-3] #remove un-useful columns

#Time-Treatment-Labeling 
Treatment_factors <- read.csv(paste("Treatment_factors_", Class, ".csv", sep=""), sep=";", header=T) #load file with treatment_factors
df$Treatment <- Treatment_factors$Treatment
df$Time <- Treatment_factors$Time
df$Labeling <- Treatment_factors$Labeling #add them to df


#Dataset Uniformation ####
df$metabolite <- factor(df$metabolite)
df$Time <- factor(df$Time)
df$Labeling <- factor(df$Labeling)
df$Treatment <- factor(df$Treatment)

vector_metabolite <- levels(factor(df$metabolite))
vector_Treatment <- levels(factor(df$Treatment))
vector_Time <- levels(factor(df$Time))

#create Subset according to metabolite - Treatment - Time ####
#metabolite
Subset <- lapply(vector_metabolite, function(i){ 
  i <- subset(df, metabolite == i)
})
names(Subset) <- vector_metabolite #add names

#Treatment
Subset_2 <- lapply(vector_metabolite, function(m){
  lapply(vector_Treatment, function(i){
    subset(Subset[[m]], Treatment == i)
  })
})
names(Subset_2) <- vector_metabolite #add names
for(i in vector_metabolite) {
  names(Subset_2[[i]]) <- vector_Treatment
} #add names

#Time
Subset_3 <- lapply(vector_metabolite, function(m){
  lapply(vector_Treatment, function(i){
    lapply(vector_Time, function(n){
      subset(Subset_2[[m]][[i]], Time == n)
    })
  })
})
names(Subset_3) <- vector_metabolite #add names
for(i in vector_metabolite) {
  names(Subset_3[[i]]) <- vector_Treatment
} #add names
for(i in vector_metabolite) {
  for(n in vector_Treatment) {
    names(Subset_3[[i]][[n]]) <- vector_Time
  }
} #add names



# Correlation plots####
#Time as colour + Treatment as FacetWrap 
p1 <- ggscatter(df, x = "mean_enrichment", y = "area", 
          color = "Time",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Enrichment %", ylab = "Relative_abundance") +
  scale_color_manual(values=c("grey77","skyblue3", "blue", "darkblue")) +
  scale_fill_manual(values=c("grey77","skyblue3", "blue", "darkblue")) +
  facet_wrap(~metabolite + Treatment, scales="free")
p1

ggsave(filename = paste("InHouse_", Class, "_Corr-matrix_(Time).pdf", sep=""), plot = last_plot(), dpi = 600, units = "cm", width = 150, height = 110, scale = 0.5)

#Treatment as colour + Time as FacetWrap 
p2 <- ggscatter(df, x = "mean_enrichment", y = "area", 
               color = "Treatment",
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "spearman",
               xlab = "Enrichment %", ylab = "Relative_abundance") +
  scale_color_manual(values=c("grey77","orange", "blue")) +
  scale_fill_manual(values=c("grey77","orange", "blue")) +
  facet_wrap(~metabolite + Time, scales="free")
p2

ggsave(filename = paste("InHouse_", Class, "_Corr-matrix_(Treatment).pdf", sep=""), plot = last_plot(), dpi = 600, units = "cm", width = 150, height = 110, scale = 0.5)

#Metabolites as FacetWrap 
p3 <- ggscatter(df, x = "mean_enrichment", y = "area", 
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "spearman",
               xlab = "Enrichment %", ylab = "Relative_abundance") +
  facet_wrap(~metabolite, scales="free")
p3

ggsave(filename = paste("InHouse_", Class, "_Corr-matrix.pdf", sep=""), plot = last_plot(), dpi = 600, units = "cm", width = 150, height = 110, scale = 0.5)

#All together
p4 <- ggscatter(df, x = "mean_enrichment", y = "area", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "spearman",
                xlab = "Enrichment %", ylab = "Relative_abundance")
p4

ggsave(filename = paste("InHouse_", Class, "_Correlation.pdf", sep=""), plot = last_plot(), dpi = 600, units = "cm", width = 75, height = 55, scale = 0.3)



#Correlation statistic test ####
#NOT WORKING YET (use different subsets)
res <- lapply(vector_metabolite, function(m){
  lapply(names(Subset_3[[m]]), function(i){ 
    lapply(names(Subset_3[[m]][[i]]), function(n){ 
      cor.test(Subset_3[[m]][[i]][[n]], method = "spearman")
    })
  })
})
res <- cor.test(df$Relative_abundance, df$Labeled, method = "pearson")
