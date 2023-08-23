library(readr)
library(ggpubr)
library(ggplot)
library(tidyverse)
library(plyr)
library(dplyr)
library(rstatix)
library(purrr)
library(agricolae)
library(reshape2)


#Histogram delta enrichment (1 for each molecule, time as x, treatment as colour)
#THINK ON HOW TO OBTAIN ARRORBARS ON MEAN ENRICHMENT


#Histogram for each molecule and time to show isotopologues enrichment 

#OA
# Read CSV ####
table <- read.csv("20230711 3-NPH acids from in-house script_IsoCor_res.tsv", sep="\t", header=T)
Class <- "OA"

#AA
# Read CSV ####
table <- read.csv("20230714_AccQ-Tag_AA_peak areas_inhouse_script 2.0_IsoCor_res.tsv", sep="\t", header=T)
Class <- "AA"



# Dataset preparation ####
#cleaning
df <- table[,-c(3,5:9)] #remove un-useful columns
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

#Uniformation and vectors creation
df$metabolite <- factor(df$metabolite)
df$Treatment <- factor(df$Treatment)
df$Time <- factor(df$Time)
df$Labeling <- factor(df$Labeling) #convert character and numerbers to factors

vector_metabolite <- levels(factor(df$metabolite))
vector_Treatment <- levels(factor(df$Treatment))
vector_Time <- levels(factor(df$Time))
vector_Labeling <- levels(factor(df$Labeling)) #write vectors for subseting

#Replace 0 with random values of 1/10 +/- 10% of min value
df["mean_enrichment"][df["mean_enrichment"]==0] <- NA #raplace meanenrichment 0 values with 1/10 ot 1/100 of min values
x <- min(df$mean_enrichment, na.rm = T)/10 #calculate min values/10
df[is.na(df)] <-  runif(sum(is.na(df)), min = x-(x/10), max = x+(x/10))



# Summary table ####
Summary_table <- ddply(df, c("metabolite", "Time", "Labeling", "Treatment"), summarise,
                       N    = sum(!is.na(mean_enrichment)),
                       mean = mean(mean_enrichment, na.rm=TRUE),
                       sd   = sd(mean_enrichment, na.rm=TRUE),
                       se   = sd / sqrt(N))



# Scatterplot ####
##    ALL 
Summary_table_L <- Summary_table %>% filter(str_detect(Labeling, "L"))
Summary_table_L$Time <- as.character(Summary_table_L$Time)
Summary_table_L$Time <- as.numeric(Summary_table_L$Time)

f1 <- ggplot(Summary_table_L, aes(x = Time, y = mean, group = Treatment, colour = Treatment)) + 
  geom_line(aes(group = Treatment)) + 
  geom_point(aes(shape = Treatment)) + 
  scale_shape_manual(values = c(15:18)) +
  scale_color_manual(values=c("grey77", "darkorange2", "skyblue3"))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, group = Treatment), width = 5) +
  theme_bw() + 
  scale_y_continuous(labels = percent) +
  scale_x_continuous(breaks=seq(0,120,15))
f2 <- f1 + facet_wrap(~metabolite, scales="free") + 
  ylab("% Enrichment") + 
  xlab("Time (min)") 

ggsave(filename = paste("InHouse_", Class, "_Scatterplot.pdf", sep=""), plot = f2, dpi = 600, units = "cm", width = 80, height = 80, scale = 0.5)



##   Significant only 
OA_Significant_Molecules <- "Aconitate|Citrate|Hydroxyglutarate|Isocitrate|Oxoglutarate"
#Pyruvate, Succinare, Oxalate and Oxaloacetate P<0.05 but removed because non-sense 

AA_Significant_Molecules <- "Aspartate|Citrulline|GABA|Glutamate|Glutamine"
#Lysine and Alanine 0.05<P<0.10 

if (Class == "AA") {
  Summary_table_L_Significant <- Summary_table_L %>% 
    filter(str_detect(metabolite, gsub('["]', '', AA_Significant_Molecules)))
} else { 
  Summary_table_L_Significant <- Summary_table_L %>% 
    filter(str_detect(metabolite, gsub('["]', '', OA_Significant_Molecules)))
}

f3 <- ggplot(Summary_table_L_Significant, aes(x = Time, y = mean, group = Treatment, colour = Treatment)) + 
  geom_line(aes(group = Treatment)) + 
  geom_point(aes(shape = Treatment)) + 
  scale_shape_manual(values = c(15:18)) +
  scale_color_manual(values=c("grey77", "darkorange2", "skyblue3"))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, group = Treatment), width = 5) +
  theme_bw() + 
  scale_y_continuous(labels = percent) +
  scale_x_continuous(breaks=seq(0,120,15))
f4 <- f3 + facet_wrap(~metabolite, scales="fixed") + 
  ylab("% Enrichment") + 
  xlab("Time (min)") 
f4

ggsave(filename = paste("InHouse_", Class, "_Scatterplot_Significant_2.pdf", sep=""), plot = f4, dpi = 600, units = "cm", width = 80, height = 50, scale = 0.3)



# Barplot ####
##    ALL ####
Summary_table_L <- Summary_table %>% filter(str_detect(Labeling, "L"))


f1 <- ggplot(Summary_table_L, aes(x = Time, y = mean, fill = Treatment)) + 
  geom_bar(stat = "identity", position = 'dodge', ) + 
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, group = Treatment), width = 0.5 , position = position_dodge(.9)) +
  scale_fill_manual(values=c("grey77", "darkorange2", "skyblue3")) +
  theme_bw() + 
  scale_y_continuous(labels = percent)
f2 <- f1 + facet_wrap(~metabolite, scales="free") + 
  ylab("% Enrichment") + 
  xlab("Time (min)") 
f2

ggsave(filename = paste("InHouse_", Class, "_Barplot.pdf", sep=""), plot = f2, dpi = 600, units = "cm", width = 80, height = 80, scale = 0.5)



##   Significant only ####
OA_Significant_Molecules <- "Aconitate|Citrate|Hydroxyglutarate|Isocitrate|Oxoglutarate"
#Pyruvate, Succinare, Oxalate and Oxaloacetate P<0.05 but removed because non-sense 

AA_Significant_Molecules <- "Aspartate|Citrulline|GABA|Glutamate|Glutamine"
#Lysine and Alanine 0.05<P<0.10 

if (Class == "AA") {
  Summary_table_L_Significant <- Summary_table_L %>% 
    filter(str_detect(metabolite, gsub('["]', '', AA_Significant_Molecules)))
} else { 
  Summary_table_L_Significant <- Summary_table_L %>% 
    filter(str_detect(metabolite, gsub('["]', '', OA_Significant_Molecules)))
}

f3 <- ggplot(Summary_table_L_Significant, aes(x = Time, y = mean, fill = Treatment)) + 
  geom_bar(stat = "identity", position = 'dodge', ) + 
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, group = Treatment), width = 0.5 , position = position_dodge(.9)) +
  scale_fill_manual(values=c("grey77", "darkorange2", "skyblue3")) +
  theme_bw() + 
  scale_y_continuous(labels = percent)
f4 <- f3 + facet_wrap(~metabolite, scales="free") + 
  ylab("% Enrichment") + 
  xlab("Time (min)") 
f4

ggsave(filename = paste("InHouse_", Class, "_Barplot_Significant_2.pdf", sep=""), plot = f4, dpi = 600, units = "cm", width = 80, height = 50, scale = 0.3)
