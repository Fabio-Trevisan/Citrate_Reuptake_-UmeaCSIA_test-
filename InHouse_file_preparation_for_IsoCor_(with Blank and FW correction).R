library(readr)
library(ggpubr)
library(tidyverse)
library(plyr)
library(dplyr)
library(rstatix)
library(purrr)
library(agricolae)
library(reshape2)

#OA
df <- read.csv("20230711 3-NPH acids from in-house script.csv", sep=";",
               header=T)
df <- df[-1,]
df <- data.frame(t(df)) #reverse dataframe 
names(df) <- t(df[1,]) #rename columns
df <- df[-1,]

df$metabolite <- row.names(df)
row.names(df) <- NULL
df$derivative <- 0
df$isotopologue <- 0
df$resolution <- 0 #add rown needed by UsoCor

write.table(df, file = "20230711 3-NPH acids from in-house script_Clean.csv", sep =";")


df1 <- melt(df, value.name = "area", variable.name = "sample", 
            id = c("metabolite", "derivative", "isotopologue", "resolution")) #melt columns

write.table(df1, file = "20230711 3-NPH acids from in-house script_IsoCor.csv", sep =";")



#AA
df <- read.csv("20230714_AccQ-Tag_AA_peak areas_inhouse_script 2.0.csv", sep=";",
               header=T)
df <- df[-1,]
df <- data.frame(t(df)) #reverse dataframe 
names(df) <- t(df[1,]) #rename columns
df <- df[-1,]

df$metabolite <- row.names(df)
row.names(df) <- NULL
df$derivative <- 0
df$isotopologue <- 0
df$resolution <- 0 #add rown needed by UsoCor

write.table(df, file = "20230714_AccQ-Tag_AA_peak areas_inhouse_script 2.0_Clean.csv", sep =";")


df1 <- melt(df, value.name = "area", variable.name = "sample", 
            id = c("metabolite", "derivative", "isotopologue", "resolution")) #melt columns

write.table(df1, file = "20230714_AccQ-Tag_AA_peak areas_inhouse_script 2.0_IsoCor.csv", sep =";")
