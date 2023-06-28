library(readr)
library(ggpubr)
library(tidyverse)
library(plyr)
library(dplyr)
library(rstatix)
library(purrr)
library(agricolae)
library(reshape2)

#Exudates times statistics ####

#REPLICATE RESULTS USING UNLABELED AS DATA VARIABLE (95% INSTEAD OF 5%) AND CHECK RELIABILITY OF RESULTS
#CHECK CORRELATION BETWEEN [C] & ENRICHMENT% AND RELATIVE ABUNDANCE & ENRICHMENT%
#RUN 2-WAY, 1-WAY ANOVA AND POST HOC TEST ON TREATMENT AND TIME OF SIGNIFICANT MOLECULES

#Read CSV ####
df <- read.csv("DATA enrichment.csv", sep=";",
               header=T)
df <- df[,-c(2,3)]


#Dataset Uniformation ####
df$Compound <- factor(df$Compound)
df$Time <- factor(df$Time)
df$Labeling <- factor(df$Labeling)
df$Treatment <- factor(df$Treatment)

vector_Compound <- levels(factor(df$Compound))
vector_Treatment <- levels(factor(df$Treatment))
vector_Time <- levels(factor(df$Time))

#USED LABELED VALUES FOR TEST!!! (% OF 13C labelled)
#create Subset according to Compound - Treatment - Time ####
#Compound
Subset <- lapply(vector_Compound, function(i){ 
  i <- subset(df, Compound == i)
})
names(Subset) <- vector_Compound #add names

#Treatment
Subset_2 <- lapply(vector_Compound, function(m){
  lapply(vector_Treatment, function(i){
    subset(Subset[[m]], Treatment == i)
  })
})
names(Subset_2) <- vector_Compound #add names
for(i in vector_Compound) {
  names(Subset_2[[i]]) <- vector_Treatment
} #add names

#Time
Subset_3 <- lapply(vector_Compound, function(m){
  lapply(vector_Treatment, function(i){
    lapply(vector_Time, function(n){
      subset(Subset_2[[m]][[i]], Time == n)
    })
  })
})
names(Subset_3) <- vector_Compound #add names
for(i in vector_Compound) {
  names(Subset_3[[i]]) <- vector_Treatment
} #add names
for(i in vector_Compound) {
  for(n in vector_Treatment) {
    names(Subset_3[[i]][[n]]) <- vector_Time
  }
} #add names



# Correlation ####
ggscatter(df, x = "Labeled", y = "Relative_abundance", 
          color = "Time",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Labeling %", ylab = "Relative_abundance") +
  scale_color_manual(values=c("darkblue","grey77","skyblue3", "blue")) +
  scale_fill_manual(values=c("darkblue","grey77","skyblue3", "blue")) +
  facet_wrap(~Compound + Treatment, scales="free")

ggsave(filename = "Corr-matrix_RelAbund_vs_Labeled_(Time).pdf", plot = last_plot(), dpi = 600, units = "cm", width = 150, height = 110, scale = 0.5)


#NOT WORKING YET (use different subsets)
res <- lapply(vector_Compound, function(m){
  lapply(names(Subset_3[[m]]), function(i){ 
    lapply(names(Subset_3[[m]][[i]]), function(n){ 
      cor.test(Subset_3[[m]][[i]][[n]], method = "spearman")
    })
  })
})
res <- cor.test(df$Relative_abundance, df$Labeled, method = "pearson")
