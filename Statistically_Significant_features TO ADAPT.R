library(readr)
library(ggpubr)
library(tidyverse)
library(plyr)
library(dplyr)
library(rstatix)
library(purrr)
library(agricolae)
library(reshape2)

#Selecting significant features 

#Read CSV ####
df <- read.csv("Exudates_ISTD_Area_.csv", sep=",",
               header=T)
df <- df[!str_detect(names(df), "Sample")]
nr. <- 58
nc. <- ncol(df)-1

#Dataset Uniformation ####
names(df) <- substring(names(df), 3)
colnames(df)[colnames(df) == ""] <- "Tr"
colnames(df)[colnames(df) == "me"] <- "Time"
c <- c("Tr","Time")
table <- melt(df, na.rm = FALSE, value.name = "Value", id = c, variable = "RI_mz")

vector_RImz <- dput(colnames(df[,-(1:2)]))

table$RI_mz <- factor(table$RI_mz)
table$Time <- factor(table$Time)
table$Tr <- factor(table$Tr)

#create Subsets according to RI_mz ####
Subsets <- lapply(vector_RImz, function(i){ 
  i <- subset(table, RI_mz == i)
})

names(Subsets) <- vector_RImz


#Assumptions ####
## 1. Homogeneity of variances
##Tr*Time
Levene_test <- lapply(split(table, table$RI_mz), function(i){
  levene_test(Value ~ Tr * Time, data = i)
})

##2. Normality
##Shapiro-Wilk test for all single Trs
SW_test <- table %>%
  group_by(Time, Tr, RI_mz) %>%
  shapiro_test(Value)
View(SW_test)
write.table(SW_test, file = "13C_ShapiroWilk_test_results_2.0.csv", quote = FALSE, sep = ";")

##3. Indipendency
#Data are indepent by experimental design!



#2way ANOVA ####
##Multiple RI_mz
TwoWay_Anova <- lapply(split(table, table$RI_mz), function(i){
  anova(lm(Value ~ Tr * Time, data = i))
})
write.table(TwoWay_Anova, file = "TwoWay_Anova_results.csv", quote = FALSE, sep = ";")



#P-Value print in original table ####
#P-Value df preparation
TwoWayAnova_PValue <- lapply(vector_RImz, function(i){
  as.data.frame(TwoWay_Anova[[i]][["Pr(>F)"]])
})
names(TwoWayAnova_PValue) <- vector_RImz #extract P_Value as list

for(i in vector_RImz) { 
  names(TwoWayAnova_PValue[[i]]) <- "P_Value"
} #rename sub-title (columns) to P-Value

vector_rows <- c("Tr", "Time", "Tr:Time", "Residuals")
for(i in vector_RImz) {   
  row.names(TwoWayAnova_PValue[[i]]) <- vector_rows
} #rename sub-title (rows) to Tr, Time, Tr:Time, Residuals

P_Value <- data.frame(TwoWayAnova_PValue) #transform list in dataframe
names(P_Value) <- vector_RImz #rename Columns
P_Value <- P_Value[-4,] #remove Row

#Combining df
df1 <- df
df1$Tr_Time <- paste(df1$Tr, df1$Time, sep="_") #Combine names Treatment and Time 
df2 <- df1[,-(1:2)] #remove Tr & Time
df3 <- as.data.frame(t(df2)) #reverse original dataframe
names(df3) <-  unlist(df3[nc.,]) #rename Columns
df3 <- df3[-nc.,] #remove Row (metabolites)
df4 <- cbind(df3, (t(P_Value))) #add column with P_Values to original table


#P Value filtering 
z <- 0.05 #filtering ratio
Filtered <- df4[df4$`Time`<z,] #P_Value filtering >z
Filter <- "Time"
#OR
Filtered <- df4[df4$`Tr`<z,] #P_Value filtering >z
Filter <- "Tr"
#OR
Filtered <- df4[df4$`Tr:Time`<z,] #P_Value filtering >z
Filter <- "TrTime"

Filtered <- Filtered[,1:nr.] #remove columns (samples)

#Cleaning
Clean_df <- data.frame(t(Filtered)) #invert row/columns
Clean_df <- data.frame(Names = row.names(Clean_df), Clean_df) #extract rownames into column
rownames(Clean_df) <- NULL #clean row names
Clean_df <- extract(Clean_df, Names, into = c("Tr", "Blank","Time"), "(.{1})(.{1})(.{1})", remove=T) #separate Treatment from Time
Clean_df <- Clean_df[!str_detect(names(Clean_df), "Blank")]

#save
write.csv(Clean_df, file = paste("Exudates_Significant", Filter, ".csv", sep="_"), row.names=FALSE)
