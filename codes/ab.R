library(readxl)
library(dplyr)
library (phyloseq)

df = read_xlsx("input_files/AB.xlsx", sheet=2)
colnames(df) [colnames(df)=="Weight change"] = "Weight_change"
df

#new coloumn that is the percentage weight change
df <- mutate(df, weight_change_percent = (Initial_weight-Final_weight)/Final_weight *100)

# Proportion of each treatment that survived 
#control 
control_group <- subset(df, Treatment == "control")
prop.table(table(control_group$Survived))["Yes"]
# 86% survived

#antibiotic group
AB_group <- subset(df, Treatment == "antibiotic")
prop.table(table(AB_group$Survived))["Yes"]
# 64% survived

# t-test for treatment vs weight change
control <- subset(df, Treatment == "control")
antibiotic <- subset(df, Treatment == "antibiotic")
t.test(control$Final_weight, antibiotic$Final_weight)
# t = 0.69098, df = 37.988, p-value = 0.4938

# Chi-square test for survival vs treatment.
table <- table(df$Survived, df$Treatment)
chi <- chisq.test(table)
chi
# X-squared = 2.5092, df = 1, p-value = 0.1132


############################
#############egg############
############################

justeggs <- df[df$Brood_stage == "Egg", ]

# t-test for treatment vs weight change
control <- subset(justeggs, Treatment == "control")
antibiotic <- subset(justeggs, Treatment == "antibiotic")
t.test(control$Weight_change, antibiotic$Weight_change)

# Chi-square test for survival vs treatment.
table <- table(justeggs$Survived, justeggs$Treatment)
chi <- chisq.test(table)
chi

