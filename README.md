TCGA Explorer
================
Kyle Dunlap
2024-02-16

## Introduction

This repository is used to store information about TCGA datasets as well
as code on how to mine through them. In this README, I will share with
you several functions that will hopefully help you mine through TCGA
data to get some useful results. These functions use a dataset compiled
from the PANCANCER subset in Xenabrowser. If you want to see how I
cleaned and put this dataset together, please check out the code in this
repository. This pre-cleaned dataset is stored on the University of Utah
Biochemistry Fileserver. It is comprised of 12,724 patients, 20,532
genes, and 2 columns at the beginning indicating the TCGA cancer code
(<https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations>)
and the abbreviated patient ID. Below is the code used to load in data,
and what the first six rows of the table looks like.

To access the data, please download from this link.

https://drive.google.com/file/d/1LzjOF6AJxkDOs8mn0xr_jeHAYgo8fpoC/view?usp=drive_link

Then, run the below code to initialize the dataset. Substitute "filename" with whatever you named the file when you downloaded it.

``` r
library(tidyverse)
all <- read_csv("filename")
all <- all |> select(!1)
```

    ## # A tibble: 6 × 20,571
    ##   Type  patient      ARHGEF10L HIF3A RNF17 RNF10 RNF11 RNF13 GTF2IP1  REM1 MTVR2
    ##   <chr> <chr>            <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl> <dbl> <dbl>
    ## 1 ACC   TCGA-OR-A5L…      6.25  4.74     0  12.3  10.2  9.91    12.6  4.81 0.747
    ## 2 ACC   TCGA-OR-A5J…      8.65  4.46     0  12.2  10.3  9.94    12.7  4.27 0.532
    ## 3 ACC   TCGA-OR-A5K…      8.08  5.22     0  13.4  10.4 10.4     12.7  1.84 0    
    ## 4 ACC   TCGA-PK-A5H…      8.53  2.44     0  12.3  10.3 10.4     12.6  5.04 0    
    ## 5 ACC   TCGA-OR-A5L…      8.92  5.85     0  12.4  10.5 11.0     13.1  2.49 0    
    ## 6 ACC   TCGA-OR-A5J…      7.43  5.55     0  12.5  10.9  9.92    12.6  4.50 0    
    ## # ℹ 20,560 more variables: RTN4RL2 <dbl>, C16orf13 <dbl>, C16orf11 <dbl>,
    ## #   FGFR1OP2 <dbl>, TSKS <dbl>, ATRX <dbl>, PMM2 <dbl>, LOC100272146 <dbl>,
    ## #   ASS1 <dbl>, NCBP1 <dbl>, ZNF709 <dbl>, ZNF708 <dbl>, RBM14 <dbl>,
    ## #   NCBP2 <dbl>, DISC1 <dbl>, CAMK1 <dbl>, RPL37 <dbl>, SPR <dbl>,
    ## #   ZNF700 <dbl>, ZNF707 <dbl>, CAMK4 <dbl>, ZNF704 <dbl>, LOC339240 <dbl>,
    ## #   GOLGA6B <dbl>, RNF115 <dbl>, RNF112 <dbl>, ZC3H14 <dbl>, SPN <dbl>,
    ## #   HMGCLL1 <dbl>, NACAP1 <dbl>, LRRTM1 <dbl>, GRIN1 <dbl>, RBMY1A3P <dbl>, …

This dataset is named “all” and will be the main input for all of the
following functions.

# Function 1 - scatter() - Makes a scatterplot of two genes.

To use this function, you input two genes and then the cancer type (ACC,
BLCA, LIHC) for example. If you want to make a large plot utilizing all
36 cancer types, you will plug in “TCGA” as the third argument.

``` r
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
```

Running this function gives the output below.

``` r
scatter(PHGDH, PSAT1, TCGA)
```

<img src="README_files/figure-gfm/function 1 result-1.png" width="75%" height="75%" />

# Function 2 - scatter_facet() - Makes scatter plots of 2 genes for all 36 cancer types.

This function takes advantage of the facet_wrap() function in ggplot2.
You will input two genes, and get 36 “mini-scatter plots”, one for each
cancer type. It will return both a table and a plot. To save space, I
only show a tiny snippet of the plot, but you can see the whole thing by
changing print(n = 3) to print(n = Inf).

``` r
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
      print(n = 3)
    
    print(plot)
}
```

Running this function gives the output below.

``` r
scatter_facet(PHGDH, PSAT1)
```

    ## # A tibble: 36 × 2
    ##   Type     cor
    ##   <chr>  <dbl>
    ## 1 KICH   0.779
    ## 2 GBMLGG 0.680
    ## 3 READ   0.663
    ## # ℹ 33 more rows

![](README_files/figure-gfm/function%202%20result-1.png)<!-- -->

# Function 3 - single_cor() Gives you the correlation value of two genes.

This function outputs a single number: The correlation value between two
genes within one cancer type or across all of them.

``` r
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
```

Running this code will give you your single value, in this case, the
correlation between PSAT1 and PHGDH over all cancer types is 0.476

``` r
single_cor(PSAT1, PHGDH, TCGA)
```

    ## # A tibble: 1 × 1
    ##   corvalue
    ##      <dbl>
    ## 1    0.476

# Function 4 - boxplot_many() - Give any number of genes and the TCGA type and get a single boxplot

This can be done pretty easily with XenaBrowser, but this method saves
clicking time and improves the graphics. You can input any number of
genes (greater than 1) and get a boxplot. As we’ve been doing, you can
input either one cancer type or all of them with “TCGA”.

``` r
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
```

Running this code with any number of genes will give you the boxplot
below. The cancer type (or TCGA) must be in the first argument.

``` r
boxplot_many(TCGA, OTC, ASS1, SLC7A5, ARG1, ARG2)
```

    ## Warning: Using an external vector in selections was deprecated in tidyselect 1.1.0.
    ## ℹ Please use `all_of()` or `any_of()` instead.
    ##   # Was:
    ##   data %>% select(genes)
    ## 
    ##   # Now:
    ##   data %>% select(all_of(genes))
    ## 
    ## See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](README_files/figure-gfm/function%204%20result-1.png)<!-- -->

# Function 5 - boxplot_onegene() - Give one gene and get a boxplot of all 36 cancer types.

This plot takes one gene as an input and gives you a boxplot of 36
boxes, each corresponding to a cancer type. This is a good way to see
what cancers your gene is highly or lowly expressed in.

``` r
boxplot_onegene <- function(Gene) {
  all |>
    select({{Gene}}, Type, sample_type) |>
    filter(sample_type != "Solid Tissue Normal") |>
    filter(sample_type != "Recurrent Tumor") |>
    filter(sample_type != "Additional - New Primary") |>
    filter(sample_type != "Additional Metastatic") |>
    ggplot(aes(x = reorder(Type, {{Gene}}, FUN=median), y = {{Gene}}, fill = Type)) +
    geom_boxplot() +
    labs(x = "TCGA Code") +
    theme_minimal() +
    theme(axis.line = element_line(linewidth = .5), 
          axis.text.x = element_text(size = 9, angle = 45, hjust = 1, face = "bold"),
          legend.position = "none")
}
```

Here, you input one gene to get the boxplot.

``` r
boxplot_onegene(PHGDH)
```

![](README_files/figure-gfm/function%205%20result-1.png)<!-- -->

# Function 6 - onegene_corloop() - Correlate over all genes after giving one gene.

This function takes one gene as an input and returns a table of the most
correlated gene to the input gene. Like previous functions, you can
choose what cancer type to look at. This is a good way to see what genes
may have some sort of co-regulation with the input gene.

``` r
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
    print(n = 5)
  
}
```

Give the input gene and cancer code. You can print as many genes as you
want, but I will only show 5.

``` r
onegene_corloop(PHGDH, TCGA) 
```

    ## # A tibble: 20,530 × 2
    ##   gene     cor
    ##   <chr>  <dbl>
    ## 1 PHGDH  1    
    ## 2 CBS    0.485
    ## 3 PSAT1  0.476
    ## 4 ERI3   0.451
    ## 5 MAP6D1 0.441
    ## # ℹ 20,525 more rows

# 7. normal_cancer - Create a Boxplot of Solid Tissue vs. Cancer in Your Type

What’s cool about this dataset is that there is also tissue next to the
diseased area that we have access to. Here, you can plot the values of a
gene in a cancer tissue vs. the healthy tissue nearby.

``` r
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
```

Give a gene and a cancer type to get a boxplot.

``` r
normal_cancer(PHGDH, BRCA)
```

![](README_files/figure-gfm/function%207%20result-1.png)<!-- -->

# 8. rank - Ranks cancer types by expression of given gene

This function will take one gene and rank its expression within the 36
available cancer types in TCGA.

``` r
rank <- function(gene) {
  gene <- ensym(gene)
  all |>
  select({{gene}}, Type) |>
  group_by(Type) |>
  summarize(
    mean = mean({{gene}}, na.rm = TRUE) 
  ) |>
  arrange(desc(mean)) |>
  print(n = 5)
}
```

Input one gene. You will only see the top 5 results.

``` r
rank(PHGDH)
```

    ## # A tibble: 36 × 2
    ##   Type    mean
    ##   <chr>  <dbl>
    ## 1 LGG     12.3
    ## 2 UCS     12.1
    ## 3 GBMLGG  12.1
    ## 4 MESO    11.7
    ## 5 PRAD    11.5
    ## # ℹ 31 more rows
