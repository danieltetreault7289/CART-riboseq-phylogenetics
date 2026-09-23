##TRANS EFF, but comparing long reads only TPMs in riboseq to total RNAseq TPMs

#install.packages(c("tidyverse", "janitor", "ggplot2", "ggprism"))
library(tidyverse)
library(dplyr)
library(janitor)
library(ggplot2)
library(ggprism)

##adjust to your working directory
path <- ##adjust to your working directory

setwd(path)

riboseq <- read.csv("10.28.25 all six sample tpms_longreadsonly.csv")
riboseq <- riboseq[,-1]
total <- read.csv("10.20.25 total RNA all six sample tpms.csv")
total <- total[,-1]


unique_to_riboseq <- riboseq %>%
  filter(!name %in% total$name)

riboseq <- riboseq %>%
  filter(name %in% total$name)

total <- total %>% 
  filter(name %in% riboseq$name)


trans_eff <- riboseq$name
trans_eff <- as.data.frame(trans_eff)
trans_eff$wt1 <- riboseq$tpm_wt1 / total$tpm_wt1
trans_eff$wt2 <- riboseq$tpm_wt2 / total$tpm_wt2
trans_eff$wt3 <- riboseq$tpm_wt3 / total$tpm_wt3
trans_eff$ydif1 <- riboseq$tpm_ydif1 / total$tpm_ydif1
trans_eff$ydif2 <- riboseq$tpm_ydif2 / total$tpm_ydif2
trans_eff$ydif3 <- riboseq$tpm_ydif3 / total$tpm_ydif3

avg_wt <- apply(trans_eff[,2:4], 1, mean)
avg_wt <- data.frame(avg_wt)
avg_wt <- mutate(avg_wt, "name" = riboseq$name)
head(avg_wt)
avg_ydif <- apply(trans_eff[,5:7], 1, mean)
avg_ydif <- data.frame(avg_ydif)
avg_ydif <- mutate(avg_ydif, "name" = riboseq$name)

combined_avg <- avg_wt %>% cbind(avg_ydif)
combined_avg <- combined_avg[,-4]

combined_avg_all <- combined_avg %>%  cbind(trans_eff)

combined_avg_all <- combined_avg_all %>% 
  filter(!grepl("^rrn", name)) %>% 
  filter(!grepl("^trn", name)) %>% 
  filter(!grepl("scr", name))

riboseq <- riboseq %>% 
  filter(!grepl("^rrn", name)) %>% 
  filter(!grepl("^trn", name)) %>% 
  filter(!grepl("scr", name))

total <- total %>% 
  filter(!grepl("^rrn", name)) %>% 
  filter(!grepl("^trn", name)) %>% 
  filter(!grepl("scr", name))

combined_avg_all <- combined_avg_all %>%
  mutate("tpm_wt1ribo" = riboseq$tpm_wt1) %>% 
  mutate("tpm_wt2ribo" = riboseq$tpm_wt2) %>% 
  mutate("tpm_wt3ribo" = riboseq$tpm_wt3) %>% 
  mutate("tpm_ydif1ribo" = riboseq$tpm_ydif1) %>% 
  mutate("tpm_ydif2ribo" = riboseq$tpm_ydif2) %>% 
  mutate("tpm_ydif3ribo" = riboseq$tpm_ydif3) %>% 
  mutate("tpm_wt1rnaseq" = total$tpm_wt1) %>% 
  mutate("tpm_wt2rnaseq" = total$tpm_wt2) %>% 
  mutate("tpm_wt3rnaseq" = total$tpm_wt3) %>% 
  mutate("tpm_ydif1rnaseq" = total$tpm_ydif1) %>% 
  mutate("tpm_ydif2rnaseq" = total$tpm_ydif2) %>% 
  mutate("tpm_ydif3rnaseq" = total$tpm_ydif3) 

combined_avg_all_filtered <- combined_avg_all %>% 
  filter(tpm_wt1ribo > 25 & tpm_wt2ribo > 25 & tpm_wt3ribo > 25 & tpm_ydif1ribo > 25 & tpm_ydif2ribo > 25 & tpm_ydif3ribo > 25) %>% 
  filter(tpm_wt1rnaseq > 25 & tpm_wt2rnaseq > 25 & tpm_wt3rnaseq > 25 & tpm_ydif1rnaseq > 25 & tpm_ydif2rnaseq > 25 & tpm_ydif3rnaseq > 25)

##VOLCANO STYLE
#run a T-test for each gene, keep only the p.value
#get a dataframe with all  gene annotations and the p.value from comparing the average of 3 WT biorep TPMs and 3 +YdiF(EQ2) biorep TPMs 
trans_eff_pvalues <- apply(trans_eff, 1, function(x) t.test((as.numeric(x[2:4])), (as.numeric(x[5:7])))$p.value)
trans_eff_pvalues_frame <- data.frame(trans_eff_pvalues) 
trans_eff_pvalues_frame <- mutate(trans_eff_pvalues_frame, "name" = riboseq$name)
head(trans_eff_pvalues_frame)

foldchange <- avg_ydif$avg_ydif / avg_wt$avg_wt
foldchange <- data.frame(foldchange)
foldchange <- log2(foldchange)
foldchange <- foldchange %>% 
  filter(!grepl("^rrn", name)) %>% 
  filter(!grepl("^trn", name)) %>% 
  filter(!grepl("scr", name))
foldchange <- mutate(foldchange, "name" = riboseq$name)

volcano <- mutate(foldchange, "p_value" = trans_eff_pvalues_frame$trans_eff_pvalues)
volcano <- mutate(volcano, trans_eff_pvalues_neglog10 = -1*log10(trans_eff_pvalues))
#volcano$name <- gsub(pattern = " gene", replacement = "", volcano$name)
volcano <- clean_names(volcano)

volcano <- volcano %>% cbind(trans_eff)
volcano <- volcano[,-1]
volcano <- volcano[,-4]
volcano <- volcano %>% cbind(combined_avg)
volcano <- volcano[,-11]

write.csv(volcano, file = "10.28.25 translation efficiency from long reads only.csv")


