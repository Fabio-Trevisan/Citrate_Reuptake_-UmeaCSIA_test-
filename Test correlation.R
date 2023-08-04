library(readr)
library(ggpubr)
library(tidyverse)
library(plyr)
library(dplyr)
library(rstatix)
library(purrr)
library(agricolae)
library(reshape2)



#REPLICATE RESULTS USING UNLABELED AS DATA VARIABLE (95% INSTEAD OF 5%) AND CHECK RELIABILITY OF RESULTS
#CHECK CORRELATION BETWEEN [C] & ENRICHMENT% AND RELATIVE ABUNDANCE & ENRICHMENT%
#RUN 2-WAY, 1-WAY ANOVA AND POST HOC TEST ON TREATMENT AND TIME OF SIGNIFICANT MOLECULES

#Read CSV ####
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
df <- table[,-c(3,5:9)] #remove un-useful columns
df <- df[df$isotopologue<=0,] #filter only isotopologue M+0
df <- df[!str_detect(df$sample, "blank"),] #remove blank samples
df <- df[!str_detect(df$sample, "std"),] #remove standards
df <- df[!str_detect(df$sample, "QC"),] #remove QCs
row.names(df) <- NULL
df <- df[,-3] #remove un-useful columns


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



# Correlation ####
ggscatter(df, x = "Enrichment", y = "Relative_abundance", 
          color = "Time",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Labeling %", ylab = "Relative_abundance") +
  scale_color_manual(values=c("darkblue","grey77","skyblue3", "blue")) +
  scale_fill_manual(values=c("darkblue","grey77","skyblue3", "blue")) +
  facet_wrap(~metabolite + Treatment, scales="free")

ggsave(filename = "Corr-matrix_RelAbund_vs_Labeled_(Time).pdf", plot = last_plot(), dpi = 600, units = "cm", width = 150, height = 110, scale = 0.5)



#NOT WORKING YET (use different subsets)
res <- lapply(vector_metabolite, function(m){
  lapply(names(Subset_3[[m]]), function(i){ 
    lapply(names(Subset_3[[m]][[i]]), function(n){ 
      cor.test(Subset_3[[m]][[i]][[n]], method = "spearman")
    })
  })
})
res <- cor.test(df$Relative_abundance, df$Labeled, method = "pearson")
