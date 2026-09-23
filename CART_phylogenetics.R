# Load packages and set working directory.
library(tidyverse)
library(treeio)
library(ggtree)
library(phangorn)
library(ggnewscale)
library(castor)
library(ggbreak)
library(ggprism)

setwd("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Daniel")

# Define function to format tables for the tree.
table_for_tree = function(phylum, gene, accessions){
  table = as.data.frame(prop.table(table(phylum, gene), margin = 1)*100)
  table_new = data.frame(cbind(gene = table$Freq[table$gene == "TRUE"]))
  rownames(table_new) = accessions
  return(table_new)
}

# Upload dataframe containing NCBI Assembly and Nucleotide accession numbers, assembly statistics, taxonomy, and gene presence.(Table S1).
df_GCF_nc_compressed = read.csv("df_daniel_rrfs_052826.csv") 

# Separate into one Nucleotide accession number per row.
df_GCF_nc = separate_rows(df_GCF_nc_compressed, nuccore, sep = ";")

# Make a smaller dataframe with only data from phyla with >10 genomes.
dfShort = df_GCF_nc_compressed %>%
  group_by(phylum) %>%
  filter(n() > 10)

n_vals = dfShort %>%
  group_by(phylum) %>%
  summarize(n = n())

# Randomly sample one assembly per phylum for the collapsed tree. This was performed once and the assembly accessions were saved to the file "random_genomes_for_tree_new.csv" for reuse. 
random = dfShort %>%
  group_by(phylum) %>%
  filter(n() > 10) %>%
  sample_n(1) %>%
  select(-X) %>%
  ungroup()

#write.csv(random, "random_genomes_dt_tree_051826.csv")

random = read.csv("random_genomes_dt_tree_051826.csv") %>%
  mutate(X = assembly) %>%
  filter(phylum != "Cyanobacteriota" & phylum != "Bacteroidota" & phylum != "Planctomycetota" & phylum != "Verrucomicrobiota")

random = column_to_rownames(random, 'X')

####### ---- Collapsed tree ---- #######

# Upload tree.
tree = read.newick("phyla_051826.tre")

# Fix GCF accession labels and drop two leaves that are duplicated phyla.
tree$tip.label = sub("_([A-Z]|_).*", "", tree$tip.label)
tree_trim = drop.tip(tree, tip = c("GCF_041929205.1", "GCF_046000445.1", "GCF_038400315.1", "GCF_033878955.1"))

# Midpoint root.
tree_mid = midpoint(tree_trim)

# Format tables for gene presence to append to the tree.

table_smpB = table_for_tree(dfShort$phylum, dfShort$smpB, random$assembly)
table_ssrA = table_for_tree(dfShort$phylum, dfShort$ssrA, random$assembly)
table_rqcH = table_for_tree(dfShort$phylum, dfShort$rqcH, random$assembly)
table_mutS2 = table_for_tree(dfShort$phylum, dfShort$mutS2, random$assembly)
table_hrpA = table_for_tree(dfShort$phylum, dfShort$hrpA, random$assembly)
table_smrB = table_for_tree(dfShort$phylum, dfShort$smrB, random$assembly)
table_rae1 = table_for_tree(dfShort$phylum, dfShort$rae1, random$assembly)


genes = data.frame(cbind(smpB = table_smpB$gene, ssrA = table_ssrA$gene, rqcH = table_rqcH$gene, mutS2 = table_mutS2$gene,  hrpA = table_hrpA$gene, smrB = table_smrB$gene, rae1 = table_rae1$gene))
rownames(genes) = rownames(random)


# Plot the tree and heatmap.
p = ggtree(tree_mid) %<+% random + 
  xlim(0, 10) + 
  geom_tiplab(aes(label=phylum), size = 4) + 
  geom_nodepoint(aes(fill = as.numeric(label)*100), size = 2, shape = 21) + 
  scale_fill_gradient(low = "white", high = "black", name = "Bootstrap\npercentage") + 
  new_scale_fill()

p1 = gheatmap(p, genes, offset = 6.7, width=2, font.size=3.5, colnames = TRUE, color=NA, colnames_angle = 45, colnames_offset_y = -0.2, colnames_offset_x = -0.15) + 
  scale_fill_viridis_c(option="G", direction = -1, name="Percent\nwith gene")


ggsave("dt_rrf_tree_052826.png", p1, width = 8.5, height = 6.4, dpi = 600, units = "in")

svg("dt_rrf_tree_052826.svg", width = 8.5, height = 6.4)
p1
dev.off()

####### ---- % phyla with each gene ---- #######

phy_table = df_GCF_nc_compressed %>%
  group_by(phylum) %>%
  summarize(total = n(), smpB = sum(smpB), ssrA = sum(ssrA), rqcH = sum(rqcH), mutS2 = sum(mutS2), hrpA = sum(hrpA), smrB = sum(smrB), rae1 = sum(rae1))

phy_table_long = rbind(phy_table, c("All phyla", sum(phy_table$total), sum(phy_table$smpB), sum(phy_table$ssrA), sum(phy_table$rqcH), sum(phy_table$mutS2), sum(phy_table$hrpA), sum(phy_table$smrB), sum(phy_table$rae1)))

write.csv(phy_table_long, "phylo_table_052826.csv", row.names = FALSE)

