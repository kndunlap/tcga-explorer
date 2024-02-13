# Load Packages test -----------------------------------------------------------

library(UCSCXenaTools)
library(tidyverse)

all <- read_csv("TCGA_all_tidy_with_metadata.csv")
all <- all |> select(!1)

# Import Files ------------------------------------------------------------
# SKIP DOWN TO #1 - scatter IF YOU ALREADY HAVE THE TIDY FILE.

downloadTCGA(project = "KIRC", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "LGG", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "KIRP", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "PANCAN", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "CHOL", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "COADREAD", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "ACC", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "CESC", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "READ", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "SARC", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "DLBC", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "PRAD", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "LUNG", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "LIHC", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "KICH", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "HNSC", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "PCPG", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "LUSC", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "TGCT", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "ESCA", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "THCA", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "LUAD", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "LAML", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "BLCA", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "SKCM", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "BRCA", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "PAAD", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "GBM", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "STAD", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "MESO", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "UVM", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "GBMLGG", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "THYM", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "UCEC", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "BRCA", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "UCS", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")
downloadTCGA(project = "COAD", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = "xena2")


# Function to clean files (2/8/24 changed to not shorten the patient code.) ------------------------------------------------


TCGAFun <- function(inputFile, codeType, outputFile) {
  library(tidyverse)
  code <- read_tsv(inputFile)
  
  code <- code |>
    pivot_longer(
      cols = starts_with("TCGA-"),
      names_to = "patient", 
      values_to = "log2exp"
    ) |>
    pivot_wider(
      names_from = sample,
      values_from = log2exp
    ) |>
    mutate(
      Type = !!codeType
    ) |>
    relocate(Type) 
  
  write_csv(code, outputFile)
}

# Clean and Preserve Files ------------------------------------------------

ACC <- TCGAFun("HiSeqV2ACC.tsv", "ACC", "ACC_tidy.csv")
BLCA <- TCGAFun("HiSeqV2BLCA.tsv", "BLCA", "BLCA_tidy.csv")
BRCA <- TCGAFun("HiSeqV2BRCA.tsv", "BRCA", "BRCA_tidy.csv")
CESC <- TCGAFun("HiSeqV2CESC.tsv", "CESC", "CESC_tidy.csv")
CHOL <- TCGAFun("HiSeqV2CHOL.tsv", "CHOL", "CHOL_tidy.csv")
COAD <- TCGAFun("HiSeqV2COAD.tsv", "COAD", "COAD_tidy.csv")
COADREAD <- TCGAFun("HiSeqV2COADREAD.tsv", "COADREAD", "COADREAD_tidy.csv")
DLBC <- TCGAFun("HiSeqV2DLBC.tsv", "DLBC", "DLBC_tidy.csv")
ESCA <- TCGAFun("HiSeqV2ESCA.tsv", "ESCA", "ESCA_tidy.csv")
GBM <- TCGAFun("HiSeqV2GBM.tsv", "GBM", "GBM_tidy.csv")
GBMLGG <- TCGAFun("HiSeqV2GBMLGG.tsv", "GBMLGG", "GBMLGG_tidy.csv")
HNSC <- TCGAFun("HiSeqV2HNSC.tsv", "HNSC", "HNSC_tidy.csv")
KICH <- TCGAFun("HiSeqV2KICH.tsv", "KICH", "KICH_tidy.csv")
KIRC <- TCGAFun("HiSeqV2KIRC.tsv", "KIRC", "KIRC_tidy.csv")
KIRP <- TCGAFun("HiSeqV2KIRP.tsv", "KIRP", "KIRP_tidy.csv")
LAML <- TCGAFun("HiSeqV2LAML.tsv", "LAML", "LAML_tidy.csv")
LGG <- TCGAFun("HiSeqV2LGG.tsv", "LGG", "LGG_tidy.csv")
LIHC <- TCGAFun("HiSeqV2LIHC.tsv", "LIHC", "LIHC_tidy.csv")
LUAD <- TCGAFun("HiSeqV2LUAD.tsv", "LUAD", "LUAD_tidy.csv")
LUNG <- TCGAFun("HiSeqV2LUNG.tsv", "LUNG", "LUNG_tidy.csv")
LUSC <- TCGAFun("HiSeqV2LUSC.tsv", "LUSC", "LUSC_tidy.csv")
MESO <- TCGAFun("HiSeqV2MESO.tsv", "MESO", "MESO_tidy.csv")
OV <- TCGAFun("HiSeqV2OV.tsv", "OV", "OV_tidy.csv")
PAAD <- TCGAFun("HiSeqV2PAAD.tsv", "PAAD", "PAAD_tidy.csv")
PCPG <- TCGAFun("HiSeqV2PCPG.tsv", "PCPG", "PCPG_tidy.csv")
PRAD <- TCGAFun("HiSeqV2PRAD.tsv", "PRAD", "PRAD_tidy.csv")
READ <- TCGAFun("HiSeqV2READ.tsv", "READ", "READ_tidy.csv")
SARC <- TCGAFun("HiSeqV2SARC.tsv", "SARC", "SARC_tidy.csv")
SKCM <- TCGAFun("HiSeqV2SKCM.tsv", "SKCM", "SKCM_tidy.csv")
STAD <- TCGAFun("HiSeqV2STAD.tsv", "STAD", "STAD_tidy.csv")
TGCT <- TCGAFun("HiSeqV2TGCT.tsv", "TGCT", "TGCT_tidy.csv")
THCA <- TCGAFun("HiSeqV2THCA.tsv", "THCA", "THCA_tidy.csv")
THYM <- TCGAFun("HiSeqV2THYM.tsv", "THYM", "THYM_tidy.csv")
UCEC <- TCGAFun("HiSeqV2UCEC.tsv", "UCEC", "UCEC_tidy.csv")
UCS <- TCGAFun("HiSeqV2UCS.tsv", "UCS", "UCS_tidy.csv")
UVM <- TCGAFun("HiSeqV2UVM.tsv", "UVM", "UVM_tidy.csv")


# Merge -------------------------------------------------------------------
library(purrr)

paths <- list.files("TSV Files from TCGA Site", pattern = "[.]csv$")
paths

files <- map(paths, read_csv)
files

all <- list_rbind(files)
write.csv(all, "TCGA_all_tidy.csv")


# Import and merge with metadata ------------------------------------------

metadata <- read_tsv("tcga _metadata.csv")
View(metadata)

metadata <- metadata |>
  select(!samples) |>
  rename(patient = sample)

all <- all |>
  inner_join(metadata)

write.csv(all, "TCGA_all_tidy_with_metadata.csv")

all |>
  count(sample_type)

# 1. scatter - Make a scatter plot of two genes with a specific cancer type --------

scatter <- function(gene1, gene2, code) {
  
  gene1 <- ensym(gene1)
  gene2 <- ensym(gene2)
  code <- ensym(code)

  filtered_data <- all |>
    select({{gene1}}, {{gene2}}, Type, sample_type) |>
    filter(if (code != "TCGA") Type == code else TRUE) |>
    filter(sample_type != "Solid Tissue Normal") |>
    filter(sample_type != "Recurrent Tumor") |>
    filter(sample_type != "Additional - New Primary") |>
    filter(sample_type != "Additional Metastatic")
  
  cor_value <- cor(filtered_data[[gene1]], filtered_data[[gene2]])
  
  ggplot(filtered_data, aes(x = {{gene1}}, y = {{gene2}})) + 
    geom_point(alpha = 0.2) +
    annotate("text", x = min(filtered_data[[gene1]]), y = max(filtered_data[[gene2]]),
             label = paste(code, "Correlation Coefficient: r =", round(cor_value, 3)),
             hjust = 0, vjust = 0, size = 5.5) +
    geom_smooth(method = "lm") +
    theme_classic() +
    theme(axis.text.x = element_text(size = 15)) +
    theme(axis.text.y = element_text(size = 15)) +
    theme(axis.title.x = element_text(size = 15)) +
    theme(axis.title.y = element_text(size = 15)) 
}

scatter(ASS1, ASL, BRCA)


# 2. scatter_facet - makes scatter plots of 2 genes for all 36 cancer types --------

scatter_facet <- function(gene1, gene2) {
  gene1 <- ensym(gene1)
  gene2 <- ensym(gene2)
  
  
  library(ggpubr)
  filtered_data <- all |>
    select({{gene1}}, {{gene2}}, Type, sample_type) |>
    filter(sample_type != "Solid Tissue Normal") |>
    filter(sample_type != "Recurrent Tumor") |>
    filter(sample_type != "Additional - New Primary") |>
    filter(sample_type != "Additional Metastatic")
  
  cor_value <- filtered_data |>
    group_by(Type) |>
    summarize(
      cor = cor({{gene1}}, {{gene2}})
    )
  
  plot <- ggplot(filtered_data, aes(x = {{gene1}}, y = {{gene2}})) + 
    geom_point(alpha = 0.3) +
    facet_wrap(~Type) +
    stat_cor(label.x.npc = .06,
             label.y.npc = 1.0,
             vjust = 1,
             size = 3) +
    geom_smooth(method = "lm")
    
    cor_value |>
      arrange(desc(cor)) |>
      print(n = Inf)
    
    print(plot)
}
scatter_facet(ASS1, SLC7A5)



# 3. Gives you the correlation value of two genes within your desired cancer type (or all) -----

single_cor <- function(Gene1, Gene2, code) {
  code <- ensym(code)
  all1 <- all |>
    select({{Gene1}}, {{Gene2}}, Type, sample_type) |>
    filter(sample_type != "Solid Tissue Normal") |>
    filter(sample_type != "Recurrent Tumor") |>
    filter(sample_type != "Additional - New Primary") |>
    filter(sample_type != "Additional Metastatic") |>
    filter(if (code != "TCGA") Type == code else TRUE) |>
    summarize(
      corvalue = cor({{Gene1}}, {{Gene2}})
    )
  return(all1)
}

single_cor(ASS1, ASL, TCGA)

# 4. boxplot_many - Give any number of genes and the TCGA type and get a single boxplot -----------

boxplot_many <- function(code, ...) {
  genes <- ensyms(...)
  code <- ensym(code)
  genes <- as.character(genes)
  all |>
    select(genes, Type, sample_type) |>
    filter(sample_type != "Solid Tissue Normal") |>
    filter(sample_type != "Recurrent Tumor") |>
    filter(sample_type != "Additional - New Primary") |>
    filter(sample_type != "Additional Metastatic") |>
    filter(if (code != "TCGA") Type == code else TRUE) |>
    pivot_longer(
      cols = (genes),
      names_to = "Gene", 
      values_to = "log2exp"
    ) |>
    mutate(Gene = factor(Gene, levels = genes)) |>
    ggplot(aes(x = Gene, y = log2exp, fill = Gene)) + 
    geom_boxplot() +
    labs(title = code) +
    theme_minimal() +
    theme(axis.line = element_line(linewidth = .5)) +
    theme(axis.text = element_text(size = 12)) +
    theme(legend.position = "none")
}

boxplot_many(TCGA, OTC, ASS1, SLC7A5, ARG1, ARG2)

# 5. boxplot_onegene -  One gene boxplot - all 36 cancer types ----------------------------------

boxplot_onegene <- function(Gene) {
  all |>
    select({{Gene}}, Type, sample_type) |>
    filter(sample_type != "Solid Tissue Normal") |>
    filter(sample_type != "Recurrent Tumor") |>
    filter(sample_type != "Additional - New Primary") |>
    filter(sample_type != "Additional Metastatic") |>
    ggplot(aes(x = reorder(Type, {{Gene}}, FUN=median), y = {{Gene}}, fill = Type)) +
    geom_boxplot() +
    theme(axis.text = element_text(size = 8.25)) +
    labs(x = "TCGA Code") +
    theme_minimal() +
    theme(axis.line = element_line(linewidth = .5)) +
    theme(legend.position = "none") 
}

boxplot_onegene(ASS1)


# 6. onegene_corloop - Correlate over all genes giving one gene. (Fixed) ---------

onegene_corloop <- function(gene, code) {
  code <- ensym(code)
  alltest <- all |>
    filter(sample_type != "Solid Tissue Normal") |>
    filter(sample_type != "Recurrent Tumor") |>
    filter(sample_type != "Additional - New Primary") |>
    filter(sample_type != "Additional Metastatic") |>
    select(!20533:20571) |>
    select(!patient) |>
    filter(if (code != "TCGA") Type == code else TRUE) |>
    relocate({{gene}}, .after = Type)
  
  alllist <- lapply(c(2:ncol(alltest)), function(x) cor(alltest[,2], alltest[,x], method = "pearson", use = "complete.obs"))
  listframe <- data.frame(alllist)
  listframe |>
    pivot_longer(
      cols = 1:ncol(listframe),
      names_to = "gene",
      values_to = "cor"
    ) |>
    arrange(desc(cor)) |>
    print(n = 25)
  
}

onegene_corloop(SLC7A5, TCGA) 


# 7. normal_cancer - Create a Boxplot of Solid Tissue vs. Cancer in Your Type --------

normal_cancer <- function(Gene, code){
  code <- ensym(code)
  title1 <- paste0("Comparison of ", code, " to nearby Non-Cancerous Tissue")
  all |>
    filter(sample_type != "Additional - New Primary") |>
    filter(sample_type != "Additional Metastatic") |>
    filter(sample_type != "Recurrent Tumor") |>
    select({{Gene}}, Type, sample_type) |>
    filter(if (code != "TCGA") Type == code else TRUE) |>
    ggplot(aes(x = sample_type, y = {{Gene}}, fill = sample_type)) + 
    geom_boxplot() +
    labs(title = title1) +
    theme_minimal() +
    theme(axis.line = element_line(linewidth = .5)) +
    theme(axis.text = element_text(size = 12)) +
    theme(axis.text.x = element_text(size = 8.5)) +
    theme(legend.position = "none")
  
}
normal_cancer(BRCA1, BRCA)


# 8. rank - Ranks cancer types by expression of given gene --------------

rank <- function(gene) {
  gene <- ensym(gene)
  all |>
  select({{gene}}, Type) |>
  group_by(Type) |>
  summarize(
    mean = mean({{gene}}, na.rm = TRUE) 
  ) |>
  arrange(desc(mean)) |>
  print(n = Inf)
}
rank(OTC)


# Messing Around ----------------------------------------------------------


all |>
  select(BRCA1, sample_type, Type) |>
  filter(Type == "BRCA") |>
  group_by(sample_type) |>
  summarize(
    mean = mean(BRCA1)
  ) |>
  ggplot(aes(x = sample_type, y = mean, fill = sample_type)) +
  geom_col()


### genelist ###
genelist <- read.csv("genelist.csv", skip = 2)
genelist <- genelist |>
  select(AC008770.3) |>
  t() 
genevec <- as.vector(genelist)

### function 6 with specific list ###
SLC7A5cor |>
  filter(gene %in% genevec) |>
  View()


