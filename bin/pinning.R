#!/usr/bin/env Rscript
package_repository <- "http://cran.us.r-project.org"
if (!require('data.table')) install.packages('data.table', repos = package_repository); library('data.table')
if (!require('dplyr')) install.packages('dplyr', repos = package_repository); library('dplyr')
if (!require('readr')) install.packages('readr', repos = package_repository); library('readr')
if (!require('tidyr')) install.packages('tidyr', repos = package_repository); library('tidyr')
if (!require('purrr')) install.packages('purrr', repos = package_repository); library('purrr')
if (!require('ggplot2')) install.packages('ggplot2', repos = package_repository); library('ggplot2')
if (!require('stringr')) install.packages('stringr', repos = package_repository); library('stringr')
if (!require('optparse')) install.packages('optparse', repos = package_repository); library('optparse')

# 2024-02-15 christian
# edit: mads

# Manual path definitions for development
vcf_folder <- "./results/variants"
pooltable_dir <- "assets/test_data/simulated_samples/sim_2x2/sampletable.tsv"
decodetable_dir <- "./assets/test_data/simulated_samples/sim_2x2/decodetable.tsv"
out_folder <- "./results/pinned"

option_list <- list(
  make_option(c("-n", "--nextflow"), help = "Script running in Nextflow [default: false]",
              type = "logical", action = "store_true", default = FALSE)
)

opt_parser <- OptionParser(option_list=option_list)
args <- parse_args(opt_parser)

# Nextflow defined paths
if (args$nextflow) {
  vcf_folder <- "variants"
  pooltable_dir <- "pooltable.tsv"
  decodetable_dir <- "decodetable.tsv"
  out_folder <- "./"
}

if (!dir.exists(out_folder)) {
  dir.create(out_folder)
}

decodetable <- read_delim(decodetable_dir,
                          delim="\t",
                          col_names = c('sample_id', 'h_pool', 'v_pool'),
                          col_types = cols(sample_id = col_character(),
                                           h_pool = col_character(),
                                           v_pool = col_character()))

pooltable <- read_delim(pooltable_dir,
                        col_select=1,
                        delim="\t",
                        col_names = c('pool_id'),
                        col_types = cols(pool_id = col_character()))

GATK.vcf.list <- 
  data.frame(GATK.vcf = list.files(vcf_folder, "*.GATK.vcf.gz$")) %>%
  mutate(pool_id = gsub(".GATK.vcf.gz", "", GATK.vcf))

Lofreq.vcf.list <- 
  data.frame(Lofreq.vcf = list.files(vcf_folder, "*.lofreq.vcf.gz$")) %>%
  mutate(pool_id = gsub(".lofreq.vcf.gz", "", Lofreq.vcf))

VCF.list <- 
  left_join(GATK.vcf.list, Lofreq.vcf.list, by = "pool_id")


if (!setequal(VCF.list$pool_id,pooltable$pool_id) || anyNA(VCF.list)) {
  stop('ERROR: Missing VCFs files')
}

#
## Get variants for each individual pool
#

collect.Pools <- list()
for(i in 1:nrow(VCF.list)){
  pool_id <- VCF.list$pool_id[i]
  gatk_filename <- file.path(vcf_folder, VCF.list$GATK.vcf[i])
  lofreq_filename <- file.path(vcf_folder, VCF.list$Lofreq.vcf[i])

  # Read and process GATK VCF file. (multi-allelic calls are on same row.)
  GATK.vcf <- 
    fread(cmd = paste0("gunzip -c ", gatk_filename, " |  grep -v '##'"))  %>%  
    rename(CHROM = 1) %>% 
    separate(10, c("GT","AD","DP","GQ","PL"), sep = ":") %>% 
    separate(AD, c("Ref.n", "Alt.n"), sep = ",") %>% 
    mutate(AF = as.numeric(Alt.n)/(as.numeric(Ref.n) + as.numeric(Alt.n))) %>% 
    mutate(var.type = ifelse(nchar(REF) > nchar(ALT), "deletion", ifelse(nchar(REF) < nchar(ALT), "insertion", "SNV"))) %>% 
    mutate(multiallelic = grepl(",",ALT)) %>% 
    mutate(var.type = ifelse(multiallelic == "TRUE", "multiallelic", var.type)) %>%  
    mutate(multiallelic = grepl(",", ALT)) %>% 
    mutate(var.type = ifelse(multiallelic == "TRUE", "multiallelic", var.type)) %>%  
    mutate(chr.pos = paste0(CHROM,":", POS)) %>% 
    mutate(var.ID = paste0(chr.pos, ":", REF, ":", ALT)) %>% 
    mutate(MQRankSum = gsub("MQRankSum=","" ,str_extract(INFO, "MQRankSum=[[:digit:]]+.[[:digit:]]+|MQRankSum=-[[:digit:]]+.[[:digit:]]+"))) %>%
    mutate(BaseQRankSum = gsub("BaseQRankSum=","" ,str_extract(INFO, "BaseQRankSum=[[:digit:]]+.[[:digit:]]+|BaseQRankSum=-[[:digit:]]+.[[:digit:]]+"))) %>% 
    mutate(SOR = gsub("SOR=","" ,str_extract(INFO, "SOR=[[:digit:]]+.[[:digit:]]+|SOR=-[[:digit:]]+.[[:digit:]]+"))) %>% 
    mutate(FS = gsub("FS=","" ,str_extract(INFO, "FS=[[:digit:]]+.[[:digit:]]+|SOR=-[[:digit:]]+.[[:digit:]]+"))) %>% 
    mutate(GT = str_count(GT, "1")) %>% 
    select(chr.pos, var.ID, GQ, GT, DP, AF,Alt.n, var.type, MQRankSum, FS, BaseQRankSum, SOR)

  # Read and process LoFreq VCF file
  VCF.names <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER" , "INFO")
  Lofreq.vcf <- fread(cmd = paste0("gunzip -c ", lofreq_filename, " |  grep -v '##'"))
  colnames(Lofreq.vcf) <- VCF.names
  
  Lofreq.vcf.annotated <- 
    Lofreq.vcf %>% 
    separate(INFO, c("DP", "AF", "SB", "DP4"), sep = ";") %>% 
    mutate(AF = gsub("AF=","", AF), DP4 = gsub("DP4=", "", DP4), DP = gsub("DP=", "", DP)) %>% 
    separate(DP4, c("FR","RR","FA","RA"), sep=",") %>%
    mutate(Alt.n = as.numeric(FA) + as.numeric(RA)) %>%
    mutate(chr.pos = paste0(CHROM,":", POS))
  
  # If multiple variants are called at the same position, merge on same row.
  # Comma separate values.
  ALT <-
    Lofreq.vcf.annotated %>%
    select(chr.pos, ALT) %>%
    group_by(chr.pos) %>%
    mutate(order = seq_along(chr.pos)) %>%
    spread(order, ALT) %>%
    unite(ALT, 2:length(.), sep = ",", na.rm = TRUE)
  
  REF <-
    Lofreq.vcf.annotated %>%
    select(chr.pos, REF) %>%
    group_by(chr.pos) %>%
    filter(!duplicated(REF)) %>%
    mutate(order = seq_along(chr.pos)) %>%
    spread(order, REF) %>%
    unite(REF, 2:length(.), sep = ",", na.rm = TRUE)
  
  ALT.n <-
    Lofreq.vcf.annotated %>%
    select(chr.pos, Alt.n) %>%
    group_by(chr.pos) %>%
    mutate(order = seq_along(chr.pos)) %>%
    spread(order, Alt.n) %>%
    unite(Alt.n, 2:length(.), sep = "/", na.rm = TRUE)
  
  DP <-
    Lofreq.vcf.annotated %>%
    select(chr.pos, DP) %>%
    group_by(chr.pos) %>%
    mutate(order = seq_along(chr.pos)) %>%
    spread(order, DP) %>%
    unite(DP, 2:length(.), sep = "/", na.rm = TRUE)
  
  AF <-
    Lofreq.vcf.annotated %>%
    select(chr.pos, AF) %>%
    group_by(chr.pos) %>%
    mutate(order = seq_along(chr.pos)) %>%
    spread(order, AF) %>%
    unite(AF, 2:length(.), sep = "/", na.rm = TRUE)
  
  Lofreq.combine <- 
    list(REF,ALT,ALT.n,DP, AF)
  
  Lofreq.vcf.annotated.combined <- 
    reduce(Lofreq.combine, left_join, by = "chr.pos")
  
  colnames(Lofreq.vcf.annotated.combined)[2:length(Lofreq.vcf.annotated.combined)] <- 
    paste0(colnames(Lofreq.vcf.annotated.combined)[2:length(Lofreq.vcf.annotated.combined)],".lofreq")

  # Combine calls by chromosome position (ONLY variants in common)
  Pool_Calls <-
    inner_join(GATK.vcf, Lofreq.vcf.annotated.combined, by = "chr.pos") %>%
    mutate(Pool_ID = pool_id)

  collect.Pools[[i]] <- Pool_Calls

  print(i)
  rm(ALT, REF, ALT.n, DP, AF, Lofreq.vcf, Lofreq.combine, Lofreq.vcf.annotated, Lofreq.vcf.annotated.combined, GATK.vcf)
}

invisible(gc())

collect.Pools_df <- do.call("rbind",collect.Pools) 

# Pool decode contains (individual) sampleID, horizontal_pool_ID, vertical_pool_ID

#
## Chr.pos match:
#

collect.chrPos <- list()

# Go through each (individual) cell position.
for(i in 1:nrow(decodetable)){
  
  # Extract all variants from horizontal pool
  Target.1 <-
    collect.Pools_df %>% 
    filter(Pool_ID == decodetable$h_pool[i]) %>% 
    mutate(FAMID = decodetable$sample_id[i]) %>% 
    select(chr.pos, everything())
  
  # Append ".A" to columns
  colnames(Target.1)[2:length(Target.1)] <- 
    paste0(colnames(Target.1)[2:length(Target.1)], ".A")
  
  # Extract all variants from vertical pool
  Target.2 <- collect.Pools_df %>% 
    filter(Pool_ID == decodetable$v_pool[i]) %>% 
    mutate(FAMID = decodetable$sample_id[i]) %>% 
    select(chr.pos, everything())
  
  # Append ".B" to columns
  colnames(Target.2)[2:length(Target.2)] <- 
    paste0(colnames(Target.2)[2:length(Target.2)], ".B")
  
  # Join on chromosome position
  collect.chrPos[[i]] <- 
    Target.1 %>% 
    inner_join(., Target.2, by = "chr.pos") 
}

# Stack
collect.chrPos_df <- 
  do.call("rbind", collect.chrPos)

# For each chromosome position:
#   Count the number of pools with variant calls (in each dimension)
#   Add genotype concordance. GT = number of alleles with alternative.
#                             var.id = chromosome, position, ref and alt.
#
#   Add allele concordance. var.id = chromosome, position, ref and alt

collect.chrPos_df.UN <- 
  collect.chrPos_df %>% 
  group_by(chr.pos) %>% 
  mutate(P1.unique = length(unique(Pool_ID.A)),
         P2.unique = length(unique(Pool_ID.B))) %>% 
  mutate(GT.concordance = ifelse(GT.A == GT.B & var.ID.A == var.ID.B,
                                 "match",
                                 "mismatch"), 
         Allele.concordance = ifelse(var.ID.A == var.ID.B,
                                     "match",
                                     "mismatch"))

# Pinnables
# Variant is pinnable if variant is unique in one dimension of pools.
collect.chrPos_df.UN.pinnables <- 
  collect.chrPos_df.UN %>%
  filter(P1.unique == 1 | P2.unique == 1) %>%
  ungroup()

# Private pinnables
# Variant is pinnable and private if variant is unique in both dimensions of pools.
Private.chrPos <- 
  collect.chrPos_df.UN.pinnables %>%
  filter(P1.unique == 1 & P2.unique == 1)

# Private sample count - count variants per sample (individual).
Private.chrPos.count <-
  Private.chrPos %>%
  count(FAMID.A)

# Find outliers based on number of private variants.
Q1 <- quantile(Private.chrPos.count$n, 0.25)
Q3 <- quantile(Private.chrPos.count$n, 0.75)
IQR <- Q3 - Q1

Private.chrPos.count.outliers <-
  Private.chrPos.count %>%
  mutate(outliers = ifelse(n < (Q1 - 1.5 * IQR) | n > (Q3 + 1.5 * IQR),
                           FAMID.A,
                           ""))

chr_pos_outlier_plot <-
  Private.chrPos.count.outliers %>%
    ggplot(., aes("", n)) +
    geom_boxplot() +
    xlab("Mapped reads") + 
    geom_text(aes(label=outliers),
              position = position_jitter(width = 0.3, height = 0),
              color = "red", size = 3) +
    ggtitle("chr:pos # private", subtitle = "IQR outliers marked")


#
## Var.ID match:
#

# Same as above but based on var.id, ie. both chromosome, position, ref and alt (only GATK output).
collect.varID <- list()
for(i in 1:nrow(decodetable)){
  Target.1 <-
    collect.Pools_df %>%
    filter(Pool_ID == decodetable$h_pool[i]) %>%
    mutate(FAMID = decodetable$sample_id[i]) %>%
    select(var.ID, everything())
  
  colnames(Target.1)[2:length(Target.1)] <-
    paste0(colnames(Target.1)[2:length(Target.1)], ".A")
  
  Target.2 <-
    collect.Pools_df %>%
    filter(Pool_ID == decodetable$v_pool[i]) %>%
    mutate(FAMID = decodetable$sample_id[i]) %>%
    select(var.ID, everything())
  
  colnames(Target.2)[2:length(Target.2)] <-
    paste0(colnames(Target.2)[2:length(Target.2)], ".B")
  
  collect.varID[[i]] <-
    Target.1 %>%
    inner_join(., Target.2, by = "var.ID") 
}

collect.varID_df <-
  do.call("rbind", collect.varID)

collect.varID_df.UN <-
  collect.varID_df %>%
  group_by(var.ID) %>%
  mutate(P1.unique = length(unique(Pool_ID.A)),
         P2.unique =  length(unique(Pool_ID.B)))

#Pinnables
collect.varID_df.UN.pinnables <-
  collect.varID_df.UN %>%
  filter(P1.unique == 1 | P2.unique == 1) %>%
  ungroup()

#Private pinnables
Private.varID <-
  collect.varID_df.UN.pinnables %>%
  filter(P1.unique == 1 & P2.unique == 1)

#Private sample count:
Private.varID.count <-
  Private.varID %>%
  count(FAMID.A)

Q1 <- quantile(Private.varID.count$n, 0.25)
Q3 <- quantile(Private.varID.count$n, 0.75)
IQR <- Q3 - Q1

Private.varID.count.outliers <-
  Private.varID.count %>%
  mutate(outliers = ifelse(n < (Q1 - 1.5 * IQR) | n > (Q3 + 1.5 * IQR), FAMID.A, ""))

var_id_outlier_plot <-
  Private.varID.count.outliers %>%
    ggplot(., aes("", n)) +
    geom_boxplot() +
    xlab("Mapped reads") + 
    geom_text(aes(label=outliers),
              position = position_jitter(width = 0.3, height = 0),
              color = "red", size = 3) +
    ggtitle("var:id # private", subtitle = "IQR outliers marked")

write_delim(collect.chrPos_df, file = file.path(out_folder, "chr_pos.tsv"), delim = "\t")
write_delim(collect.chrPos_df.UN, file = file.path(out_folder,"chr_pos_unique.tsv"), delim = "\t")
write_delim(collect.chrPos_df.UN.pinnables, file = file.path(out_folder,"chr_pos_unique_pin.tsv"), delim = "\t")

write_delim(collect.varID_df, file = file.path(out_folder ,"var_id.tsv"), delim = "\t")
write_delim(collect.varID_df.UN, file = file.path(out_folder, "var_id_unique.tsv"), delim = "\t")
write_delim(collect.varID_df.UN.pinnables, file = file.path(out_folder, "var_id_unique_pin.tsv"), delim = "\t")

ggsave(file.path(out_folder, "var_id_outliers.pdf"), plot = var_id_outlier_plot)
ggsave(file.path(out_folder, "chr_pos_outlier.pdf"), plot = chr_pos_outlier_plot)
