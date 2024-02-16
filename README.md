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
and the abbreviated patient ID. Below is the code used to very quickly
do a final clean on this sheet, and what the output looks like.

``` r
library(tidyverse)
load("all.RData")
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

![](README_files/figure-gfm/function%201%20result-1.png)<!-- -->

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
