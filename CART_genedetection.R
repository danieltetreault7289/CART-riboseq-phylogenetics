library(tidyverse)
library(janitor)
library(jsonlite)

setwd("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Daniel")

# Define function to read in HMMER files.
read_nhmmer = function(file){
  df = read.table(file, skip = 2, sep = "", colClasses = c("character", rep("NULL", 3), rep("character", 11), rep("NULL", 8)), fill = TRUE, row.names = NULL)
  df = df[,1:12]
  df = replace(df, df=='', NA)
  colnames(df) = c("accessions", "hmmfrom", "hmmto", "alifrom", "alito", "envfrom", "envto", "sqlen", "strand", "eval", "score", "bias")
  df = drop_na(df)
  df = df %>% mutate_at(c("hmmfrom", "hmmto", "alifrom", "alito", "envfrom", "envto", "sqlen", "eval", "score", "bias"), as.numeric)
  
  return(df)
}

read_phmmer = function(file) {
  df = data.frame(readLines(file)) %>%
  slice(4:n()) %>%
  separate(readLines.file., into = c(as.character(1:34)), sep = " +") %>%
  select(-2,-3,-4) %>%
  replace(is.na(.), " ") %>% #Replace NAs with spaces.
  unite("description", 16:31, sep = " ") %>%
  mutate_if(is.character, str_trim) %>% #Trim spaces off description.
  filter(`1` != "#") %>%
  mutate(species = str_extract(description, "(?<=\\[).*(?=\\])")) %>% #Extract species name from between [] in description. 
  `colnames<-`(c("accession", "full_eval", "full_score", "full_bias", "dom_eval", "dom_score", "dom_bias", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "description", "species")) %>% #Rename columns.
  mutate_at(vars(full_eval, full_score, full_bias, dom_eval, dom_score, dom_bias), as.numeric)
  return(df)
  }

# Upload assembly metadata for NCBI assemblies and select relevant data columns.
df_metadata = read_tsv("ref_prok_genome_metadata_050526.tsv") %>%
  mutate(taxID  = as.character(`Organism Taxonomic ID`)) %>%
  select(-`Annotation Name`, -`Assembly Release Date`, -`WGS project accession`, -`Assembly Paired Assembly Accession`, -`Organism Taxonomic ID`) %>%
  clean_names()

# Upload lineages acquired from NCBI taxdump and taxonkit based on each assembly's taxid.
lineages = read.csv("taxonomy_full_050826.tsv", col.names = "taxID", header = FALSE) %>% 
  separate(taxID, into = c("organism", "domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "merge") %>% 
  separate(organism, into = c("taxID", "organism"), sep = "\\s", extra = "merge") %>%
  drop_na()


df_GCF_tax_select = left_join(df_metadata, lineages, by = join_by(tax_id == taxID)) 

# Upload NCBI RefSeq Assembly database accessions (beginning with "GCF_") and corresponding NCBI Nucleotide (nuccore) database accessions (beginning with "NC_" or "NZ_".

df_GCF_nc = read.csv("GCF_accessions_bact_clean_051226.csv")

# Join the NCBI accessions with their metadata.
df_GCF_nc_tax = left_join(df_GCF_nc, df_GCF_tax_select, by = join_by(assembly == assembly_accession)) %>%
  mutate(check_m_completeness = as.numeric(check_m_completeness)) %>%
  mutate(check_m_contamination = as.numeric(check_m_contamination)) 

# Some were missing taxonomy data. Got this from NCBI datasets. Join with other tax data.
lines = readLines("missing_taxonomies.jsonl")
lines = lapply(lines, fromJSON)
lines = lapply(lines, unlist)
df_GCF_tax_missing = bind_rows(lines) %>%
  drop_na(taxonomy.classification.phylum.name) %>%
  select(tax_id = query, domain = taxonomy.classification.domain.name, phylum = taxonomy.classification.phylum.name, class = taxonomy.classification.class.name, order = taxonomy.classification.order.name, family = taxonomy.classification.family.name, genus = taxonomy.classification.genus.name, species = taxonomy.classification.species.name)


df_GCF_nc_tax = left_join(df_GCF_nc_tax, df_GCF_tax_missing, by = join_by(tax_id)) %>%
  mutate(phylum = coalesce(phylum.x, phylum.y),
         class = coalesce(class.x, class.y),
         order = coalesce(order.x, order.y),
         family = coalesce(family.x, family.y),
         genus = coalesce(genus.x, genus.y),
         species = coalesce(species.x, species.y),
         .keep = "unused") 

# ssrA

df_ssrA_some = read_nhmmer("ssrA_HMMER_5e2.csv")
df_ssrA_alpha = read_nhmmer("ssrA_alpha_HMMER_5e2.csv")

df_ssrA = rbind(df_ssrA_some, df_ssrA_alpha)

# Import the GCF+prot accessions for hrpA and rae1, and merge with length data.
df_hrpA_gcf = read_tsv("hrpA_GCF_WP_accs.tsv", col_names = c("assembly", "prot_acc")) # Pre-filtered for lengths >1200 aa!
df_rae1_gcf = read_tsv("rae1_GCF_WP_accs.tsv", col_names = c("assembly", "prot_acc")) 
df_smpB_gcf = read_tsv("smpB_GCF_WP_accs.tsv", col_names = c("assembly", "prot_acc")) 
df_mutS2_gcf = read_tsv("mutS2_GCF_WP_accs.tsv", col_names = c("assembly", "prot_acc")) # Pre-filtered for full_score >400!
df_rqcH_gcf = read_tsv("rqcH_GCF_WP_accs.tsv", col_names = c("assembly", "prot_acc")) 
df_smrB_gcf = read_tsv("smrB_GCF_WP_accs.tsv", col_names = c("assembly", "prot_acc")) # Pre-filtered for full_score >100!


df = data.frame(cbind(df_GCF_nc_tax, smpB = (df_GCF_nc_tax$assembly %in% df_smpB_gcf$assembly), ssrA = (df_GCF_nc_tax$nuccore %in% df_ssrA$accessions), rqcH = (df_GCF_nc_tax$assembly %in% df_rqcH_gcf$assembly), mutS2 = (df_GCF_nc_tax$assembly %in% df_mutS2_gcf$assembly), hrpA = (df_GCF_nc_tax$assembly %in% df_hrpA_gcf$assembly), smrB = (df_GCF_nc_tax$assembly %in% df_smrB_gcf$assembly), rae1 = (df_GCF_nc_tax$assembly %in% df_rae1_gcf$assembly)))


# Multiple contigs, scaffolds, or chromosomes can make up an assembly. Identify gene hits across all contigs/scaffolds/chromosomes in each assembly. 
df_presence = df %>%
  group_by(assembly) %>%
  summarise(across(smpB:rae1, ~sum(.))) %>%
  mutate_if(is.numeric, ~1 * (. != 0))

bools = ifelse(df_presence[-1] == 1,"TRUE","FALSE")
bools = data.frame(cbind(assembly = df_presence$assembly, bools))

# Save the names of all Nucleotide accession numbers that were searched in each assembly.
names =  df %>%
  distinct(nuccore, .keep_all = TRUE) %>%
  group_by(assembly) %>%
  summarize(nuccore=paste(nuccore, collapse=";"))

# Prepare the final dataframe with Assembly and Nucleotide accession numbers, assembly statistics, taxonomy, and gene presence.
# Remove archaeal genomes. Filter based on CheckM completeness and contamination. Remove any duplicate genomes.
df_distinct = inner_join(names, bools, by = "assembly") %>%
  left_join(df[,1:21], by = "assembly", multiple = "any") %>%
  select(-X, -nuccore.y) %>%
  rename(nuccore = nuccore.x) %>% 
  filter(check_m_completeness > 80) %>%
  filter(check_m_contamination < 10) %>%
  distinct(assembly, .keep_all = TRUE) %>%
  drop_na(phylum) %>%
  mutate(across(c(smpB:rae1), as.logical))
  
# Rename taxonomy to be consistent with Coleman et al 2021. Same as Katrina's paper for now.
df_distinct$phylum[df_distinct$phylum == "delta/epsilon subdivisions"] = "Pseudomonadota"
df_distinct$phylum[df_distinct$phylum == "Pseudomonadota"] = "Proteobacteria"
df_distinct$phylum[df_distinct$phylum == "Terrabacteria group"] = df_distinct$class[df_distinct$phylum == "Terrabacteria group"]
df_distinct$phylum[df_distinct$phylum == "Bacillota"] = "Firmicutes"
#df_distinct$phylum[df_distinct$phylum == "Actinomycetota"] = "Actinobacteriota"
df_distinct$phylum[df_distinct$phylum == "Abditibacteriota"] = "Armatimonadota"
df_distinct$phylum[df_distinct$phylum == "Aquificota"] = "Aquificota + Campylobacterota + Deferribacterota"
df_distinct$phylum[df_distinct$phylum == "Campylobacterota"] = "Aquificota + Campylobacterota + Deferribacterota"
df_distinct$phylum[df_distinct$phylum == "Deferribacterota"] = "Aquificota + Campylobacterota + Deferribacterota"
df_distinct$phylum[df_distinct$phylum == "Thermodesulfobacteriota"] = "Desulfuromonadota + Desulfobacterota"

df_distinct$phylum[df_distinct$phylum == "Calditrichota"] = "FCB group"
df_distinct$phylum[df_distinct$phylum == "Fibrobacterota"] = "FCB group"
df_distinct$phylum[df_distinct$phylum == "Bacteroidota"] = "FCB group"

df_distinct$phylum[df_distinct$phylum == "Verrucomicrobiota"] = "PVC group"
df_distinct$phylum[df_distinct$phylum == "Planctomycetota"] = "PVC group"
df_distinct$phylum[df_distinct$phylum == "Chlamydiota"] = "PVC group"

df_distinct$phylum[df_distinct$phylum == "Cyanobacteriota"] = "Cyanobacteriota/Melainabacteria group"

df_distinct$phylum[df_distinct$phylum == "Thermomicrobiota"] = "Chloroflexota"
df_distinct$phylum[df_distinct$phylum == "Proteobacteria"] = df_distinct$class[df_distinct$phylum == "Proteobacteria"]


df_distinct %>%
  summarize(total = n(), smpB = sum(smpB), ssrA = sum(ssrA), rqcH = sum(rqcH), mutS2 = sum(mutS2), hrpA = sum(hrpA), smrB = sum(smrB), rae1 = sum(rae1))


# Write Table S#.
write.csv(df_distinct, file = "df_daniel_rrfs_052826.csv")

