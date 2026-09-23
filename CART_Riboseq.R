#install.packages(c("tidyverse", "janitor", "ggplot2", "ggprism"))
library(tidyverse)
library(dplyr)
library(janitor)
library(ggplot2)
library(ggprism)

##adjust to your working directory
path <- ##adjust to your working directory

setwd(path)

#-------------------------------------------------------------------------------
#analyzing all ribosome footprints larger than 23 nt
#analyze by TPM first

#have to read in .tsv files to avoid #NAME error in variant change identity for deletions and insertions, which use "-" and "+" characters in the raw Geneious export
wt1 <- read.delim('HF1 assembled to 168 Annotations.tsv', sep = "\t") %>% 
  clean_names() %>% select(name, tpm) %>% rename(tpm_wt1 = tpm)
wt2 <- read.delim('HF2 assembled to 168 Annotations.tsv', sep = "\t") %>% 
  clean_names() %>% select(name, tpm) %>% rename(tpm_wt2 = tpm)
wt3 <-read.delim('HF3 assembled to 168 Annotations.tsv', sep = "\t") %>% 
  clean_names() %>% select(name, tpm) %>% rename(tpm_wt3 = tpm)
ydif1 <-read.delim('HF4 assembled to 168 Annotations.tsv', sep = "\t") %>% 
  clean_names() %>% select(name, tpm) %>% rename(tpm_ydif1 = tpm)
ydif2 <-read.delim('HF5 assembled to 168 Annotations.tsv', sep = "\t") %>% 
  clean_names() %>% select(name, tpm) %>% rename(tpm_ydif2 = tpm)
ydif3 <-read.delim('HF6 assembled to 168 Annotations.tsv', sep = "\t") %>% 
  clean_names() %>% select(name, tpm) %>% rename(tpm_ydif3 = tpm)

combined <- wt1 %>%  bind_cols(tpm_wt2 = wt2$tpm_wt2) %>% 
  bind_cols(tpm_wt3 = wt3$tpm_wt3) %>% 
  bind_cols(tpm_ydif1 = ydif1$tpm_ydif1) %>% 
  bind_cols(tpm_ydif2 = ydif2$tpm_ydif2) %>% 
  bind_cols(tpm_ydif3 = ydif3$tpm_ydif3)
combined$name <- gsub(pattern = " gene", replacement = "", combined$name)

write.csv(combined, file = "10.28.25 all six sample tpms.csv")



avg_wt <- apply(combined[,2:4], 1, mean)
avg_wt <- data.frame(avg_wt)
avg_wt <- mutate(avg_wt, "name" = wt1$name)
avg_wt$name <- gsub(pattern = " gene", replacement = "", avg_wt$name)
head(avg_wt)
avg_ydif <- apply(combined[,5:7], 1, mean)
avg_ydif <- data.frame(avg_ydif)
avg_ydif <- mutate(avg_ydif, "name" = wt1$name)
avg_ydif$name <- gsub(pattern = " gene", replacement = "", avg_ydif$name)

combined_avg <- avg_wt %>% cbind(avg_ydif)
combined_avg <- combined_avg[,-4]

ggplot(data = combined_avg, aes(x = avg_wt, y = avg_ydif)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "10.28.25 riboseq, average TPMs",
       x = "WT (TPM)",
       y = "+YdiF(EQ2) (TPM)") +
  theme_prism()

combined_avg_filtered <- combined_avg %>%
  filter(!grepl("^rrn", name))

combined <- combined %>%
  filter(!grepl("^rrn", name))

combined_avg_final <- combined_avg_filtered %>% 
  bind_cols(combined)

combined_avg_final <- combined_avg_final[,-4]

write.csv(combined_avg_final, file = "10.28.25 combined_avg_filtered_tpm.csv")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#>69 nt footprints only

#install.packages(c("tidyverse", "janitor", "ggplot2", "ggprism"))
library(tidyverse)
library(dplyr)
library(janitor)
library(ggplot2)
library(ggprism)


##adjust to your working directory
path <- ##adjust to your working directory

setwd(path)


wt1 <- read.delim('HF1_largeonly assembled to 168 Annotations.tsv', sep = "\t") %>% 
  clean_names() %>% select(name, tpm) %>% rename(tpm_wt1 = tpm)
wt2 <- read.delim('HF2_largeonly assembled to 168 Annotations.tsv', sep = "\t") %>% 
  clean_names() %>% select(name, tpm) %>% rename(tpm_wt2 = tpm)
wt3 <-read.delim('HF3_largeonly assembled to 168 Annotations.tsv', sep = "\t") %>% 
  clean_names() %>% select(name, tpm) %>% rename(tpm_wt3 = tpm)
ydif1 <-read.delim('HF4_largeonly assembled to 168 Annotations.tsv', sep = "\t") %>% 
  clean_names() %>% select(name, tpm) %>% rename(tpm_ydif1 = tpm)
ydif2 <-read.delim('HF5_largeonly assembled to 168 Annotations.tsv', sep = "\t") %>% 
  clean_names() %>% select(name, tpm) %>% rename(tpm_ydif2 = tpm)
ydif3 <-read.delim('HF6_largeonly assembled to 168 Annotations.tsv', sep = "\t") %>% 
  clean_names() %>% select(name, tpm) %>% rename(tpm_ydif3 = tpm)

combined <- wt1 %>%  bind_cols(tpm_wt2 = wt2$tpm_wt2) %>% 
  bind_cols(tpm_wt3 = wt3$tpm_wt3) %>% 
  bind_cols(tpm_ydif1 = ydif1$tpm_ydif1) %>% 
  bind_cols(tpm_ydif2 = ydif2$tpm_ydif2) %>% 
  bind_cols(tpm_ydif3 = ydif3$tpm_ydif3)
combined$name <- gsub(pattern = " gene", replacement = "", combined$name)

write.csv(combined, file = "10.28.25 all six sample tpms_longreadsonly.csv")


avg_wt <- apply(combined[,2:4], 1, mean)
avg_wt <- data.frame(avg_wt)
avg_wt <- mutate(avg_wt, "name" = wt1$name)
avg_wt$name <- gsub(pattern = " gene", replacement = "", avg_wt$name)
head(avg_wt)
avg_ydif <- apply(combined[,5:7], 1, mean)
avg_ydif <- data.frame(avg_ydif)
avg_ydif <- mutate(avg_ydif, "name" = wt1$name)
avg_ydif$name <- gsub(pattern = " gene", replacement = "", avg_ydif$name)
head(avg_ydif)

combined_avg <- avg_wt %>% cbind(avg_ydif)
combined_avg <- combined_avg[,-4]

combined_avg <- combined_avg %>%
  filter(!grepl("^rrn", name))

combined <- combined %>%
  filter(!grepl("^rrn", name))

combined_avg_final <- combined_avg %>% 
  bind_cols(combined)

combined_avg_final <- combined_avg_final[,-4]

write.csv(combined_avg_final, file = "10.28.25 combined_avg_long_footprints_only.csv")

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#footprint length analysis
##footprint length analysis BUT only with "used reads" after mapping!!
path <- ##adjust to your working directory

setwd(path)

wt1_lengths <- read.csv("HF1 assembled to 168 Used Reads lengths graph.csv") %>% 
  clean_names() %>% rename(length = sequence_length) %>% rename(sequences = number_of_sequences)
wt2_lengths <- read.csv("HF2 assembled to 168 Used Reads lengths graph.csv") %>% 
  clean_names() %>% rename(length = sequence_length) %>% rename(sequences = number_of_sequences)
wt3_lengths <- read.csv("HF3 assembled to 168 Used Reads lengths graph.csv") %>% 
  clean_names() %>% rename(length = sequence_length) %>% rename(sequences = number_of_sequences)
ydif1_lengths <- read.csv("HF4 assembled to 168 Used Reads lengths graph.csv") %>% 
  clean_names() %>% rename(length = sequence_length) %>% rename(sequences = number_of_sequences)
ydif2_lengths <- read.csv("HF5 assembled to 168 Used Reads lengths graph.csv") %>% 
  clean_names() %>% rename(length = sequence_length) %>% rename(sequences = number_of_sequences)
ydif3_lengths <- read.csv("HF6 assembled to 168 Used Reads lengths graph.csv") %>% 
  clean_names() %>% rename(length = sequence_length) %>% rename(sequences = number_of_sequences)

all_lengths <- wt1_lengths$length
all_lengths <- as.data.frame(all_lengths)
all_lengths$wt1 <- wt1_lengths$sequences
all_lengths$wt2 <- wt2_lengths$sequences
all_lengths$wt3 <- wt3_lengths$sequences
all_lengths$ydif1 <- ydif1_lengths$sequences
all_lengths$ydif2 <- ydif2_lengths$sequences
all_lengths$ydif3 <- ydif3_lengths$sequences


#reshape to long format for geom_bar
all_lengths_long <- all_lengths %>%
  pivot_longer(
    cols = c(wt1, wt2, wt3, ydif1, ydif2, ydif3),
    names_to = "sample",
    values_to = "count"
  ) %>%
  mutate(group = ifelse(grepl("^wt", sample), "WT", "YdiF"))

ggplot(all_lengths_long, aes(x = all_lengths, y = count, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5) +
  scale_fill_manual(values = c("WT" = "darkcyan", "YdiF" = "maroon4")) +
  theme_classic() +
  labs(
    x = "footprint length (nt)",
    y = "# of footprint sequences",
    fill = "sample"
  ) +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(filename = "10.28.25 Footprint seq. lengths after trimming and mapping WT vs YdiF.tif",
       width = 6, height = 3,
       bg = 'white')

ggsave(filename = "10.28.25 Footprint seq. lengths after trimming and mapping WT vs YdiF.jpg",
       width = 6, height = 3,
       bg = 'white')


############################################################################
############################################################################

############################################################################

#making an ssrA coverage figure, like a clean IGV

#install.packages(c("tidyverse", "janitor", "ggplot2", "ggprism"))
library(tidyverse)
library(dplyr)
library(janitor)
library(ggplot2)
library(ggprism)

##adjust to your working directory
path <- ##adjust to your working directory

setwd(path)

wt1coverage <- read.csv("HF1_largeonly_ssrAcoverage.csv")
wt2coverage <- read.csv("HF2_largeonly_ssrAcoverage.csv")
wt3coverage <- read.csv("HF3_largeonly_ssrAcoverage.csv")
ydif1coverage <- read.csv("HF4_largeonly_ssrAcoverage.csv")
ydif2coverage <- read.csv("HF5_largeonly_ssrAcoverage.csv")
ydif3coverage <- read.csv("HF6_largeonly_ssrAcoverage.csv")


#####capturing +/- 100 bp around the ssrA CDS
wt1ssrA <- wt1coverage[3450612:3451171,]
wt2ssrA <- wt2coverage[3450612:3451171,]
wt3ssrA <- wt3coverage[3450612:3451171,]
ydif1ssrA <- ydif1coverage[3450612:3451171,]
ydif2ssrA <- ydif2coverage[3450612:3451171,]
ydif3ssrA <- ydif3coverage[3450612:3451171,]

largeonlyssracov <- wt1ssrA %>% 
  left_join(wt2ssrA, by = "Position") %>% 
  left_join(wt3ssrA, by = "Position") %>% 
  left_join(ydif1ssrA, by = "Position") %>% 
  left_join(ydif2ssrA, by = "Position") %>% 
  left_join(ydif3ssrA, by = "Position") 

colnames(largeonlyssracov) <- c("position", "wt1", "wt2","wt3","ydif1","ydif2","ydif3")

write.csv(largeonlyssracov, file = "ssrA coverage of all samples large reads only.csv")




wt1allcov <- read.csv("HF1_all_coverage.csv")
wt2allcov <- read.csv("HF2_all_coverage.csv")
wt3allcov <- read.csv("HF3_all_coverage.csv")
ydif1allcov <- read.csv("HF4_all_coverage.csv")
ydif2allcov <- read.csv("HF5_all_coverage.csv")
ydif3allcov <- read.csv("HF6_all_coverage.csv")

wt1allssra <- wt1allcov[3450612:3451171,]
wt2allssra <- wt2allcov[3450612:3451171,]
wt3allssra <- wt3allcov[3450612:3451171,]
ydif1allssra <- ydif1allcov[3450612:3451171,]
ydif2allssra <- ydif2allcov[3450612:3451171,]
ydif3allssra <- ydif3allcov[3450612:3451171,]

allssracov <- wt1allssra %>% 
  left_join(wt2allssra, by = "Position") %>% 
  left_join(wt3allssra, by = "Position") %>% 
  left_join(ydif1allssra, by = "Position") %>% 
  left_join(ydif2allssra, by = "Position") %>% 
  left_join(ydif3allssra, by = "Position") 

colnames(allssracov) <- c("position", "wt1", "wt2","wt3","ydif1","ydif2","ydif3")

write.csv(allssracov, file = "ssrA coverage of all samples all read lengths.csv")
