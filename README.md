# Kyle note - need to unclean the "all" dataset to have their full code!! Or else can't join with metadata.
# Introduction
This repository is used to store information about TCGA datasets as well as code on how to mine through them.
In this README, I will share with you six functions that will hopefully help you mine through TCGA data to get some useful results.
These functions use a dataset compiled from the PANCANCER subset in Xenabrowser. If you want to see how I cleaned and put this dataset together, please check out the code in this repository. This pre-cleaned dataset is stored on the University of Utah Biochemistry Fileserver. It is comprised of 12,724 patients, 20,532 genes, and 2 columns at the beginning indicating the TCGA cancer code (https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations) and the abbreviated patient ID. Below is the code used to very quickly do a final clean on this sheet, and what the output looks like.

```
all <- read.csv("TCGA_all_tidy.csv")
all <- as.tibble(all)
all <- all |> select(!1)
```

<img width="897" alt="image" src="https://github.com/kndunlap/tcga-explorer/assets/61035909/526fca6a-9f26-416a-8795-32402f90eb22">

This dataset is termed "all" and will be the main input for all six of the following functions.

# Function 1 - scatter() - Makes a scatterplot of two genes.
To use this function, you input two genes and then the cancer type (ACC, BLCA, LIHC) for example. If you want to make a large plot utilizing all 36 cancer types, you will plug in "TCGA" as the third argument.

```
scatter <- function(gene1, gene2, code) {
  filtered_data <- all |>
  select({{gene1}}, {{gene2}}, Type) |>
  filter(if (code != "TCGA") Type == code else TRUE)
  
  cor_value <- cor(filtered_data[[gene1]], filtered_data[[gene2]])
  
  ggplot(filtered_data, aes(x = !!sym(gene1), y = !!sym(gene2))) + 
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
Using the code below gives the output below.
```
scatter("SLC7A5", "FOXM1", "TCGA")
```

<img width="1325" alt="image" src="https://github.com/kndunlap/tcga-explorer/assets/61035909/c6872988-0dd6-4416-b81b-e6256f1b19d5">

Changing "TCGA" to "BRCA" then shows you just the breast cancer samples, instead of the entire dataset.

```
scatter("SLC7A5", "FOXM1", "BRCA")
```
<img width="1324" alt="image" src="https://github.com/kndunlap/tcga-explorer/assets/61035909/4803d98f-4800-491b-81d1-d71e34883bee">

# Function 2 - scatter_facet() - Makes scatter plots of 2 genes for all 36 cancer types.
This function takes advantage of the facet_wrap() function in ggplot2. You will input two genes, and get 36 "mini-scatter plots", one for each cancer type.

```
scatter_facet <- function(gene1, gene2) {
  library(ggpubr)
  filtered_data <- all |>
    select({{gene1}}, {{gene2}}, Type) 
  
  cor_value <- cor(filtered_data[[gene1]], filtered_data[[gene2]])
  
  ggplot(filtered_data, aes(x = !!sym(gene1), y = !!sym(gene2))) + 
    geom_point(alpha = 0.5) +
    facet_wrap(~Type) +
    stat_cor(label.x.npc = .06,
             label.y.npc = 1.0,
             vjust = 1,
             size = 3) +
    geom_smooth(method = "lm")
}
```
```
scatter_facet("ASL", "ASS1")
```
<img width="1325" alt="image" src="https://github.com/kndunlap/tcga-explorer/assets/61035909/45616744-5ce1-4281-a0e6-122ad6229071">

# Function 3 - single_cor() Gives you the correlation value of two genes.
This function outputs a single number: The correlation value between two genes within one cancer type or across all of them.
```
single_cor <- function(Gene1, Gene2, code) {
  all1 <- all |>
  select({{Gene1}}, {{Gene2}}, Type) |>
  filter(if (code != "TCGA") Type == code else TRUE) |>
  summarize(
    corvalue = cor({{Gene1}}, {{Gene2}})
  )
  return(all1)
}
```
The below function call returns the correlation between ASL and ASS1 in SKCM (skin cancer melanoma)

```
single_cor(ASL, ASS1, "SKCM")
```
<img width="134" alt="image" src="https://github.com/kndunlap/tcga-explorer/assets/61035909/4f0a3419-88cd-47a7-b72f-56727685de26">

And this call returns the correlation, taking all patients into account. This is done by making "TCGA" the third argument.
```
single_cor(ASL, ASS1, "TCGA")
```
<img width="135" alt="image" src="https://github.com/kndunlap/tcga-explorer/assets/61035909/257199de-4ce6-4585-a61e-3b036a78988c">

# Function 4 - boxplot_many() - Give any number of genes and the TCGA type and get a single boxplot 
This can be done pretty easily with XenaBrowser, but this method saves clicking time and improves the graphics. You can input any number of genes (greater than 1) and get a boxplot. As we've been doing, you can input either one cancer type or all of them with "TCGA".

```
boxplot_many <- function(code, ...) {
  genes <- c(...)
  all |>
  select(genes, Type) |>
  filter(if (code != "TCGA") Type == code else TRUE) |>
  pivot_longer(
    cols = (genes),
    names_to = "Gene", 
    values_to = "log2exp"
  ) |>
  ggplot(aes(x = reorder(Gene, log2exp), y = log2exp, fill = Gene)) + 
  geom_boxplot(varwidth = TRUE) +
  labs(title = code, x = "Gene List", y = "Log2 Expression") +
  theme_minimal() +
  theme(axis.line = element_line(linewidth = .5)) +
  theme(axis.text = element_text(size = 12)) +
  theme(legend.position = "none")
}
```
```
boxplot_many("TCGA", "OTC", "ASS1", "SLC7A5", "ARG1", "ARG2")
```
<img width="697" alt="image" src="https://github.com/kndunlap/tcga-explorer/assets/61035909/2353138c-eff0-4c2d-a899-ac5cc96b12fd">

# Function 5 - boxplot_onegene() -  Give one gene and get a boxplot of all 36 cancer types.
This plot takes one gene as an input and gives you a boxplot of 36 boxes, each corresponding to a cancer type. This is a good way to see what cancers your gene is highly or lowly expressed in.
```
boxplot_onegene <- function(Gene) {
  all |>
  select({{Gene}}, Type) |>
  ggplot(aes(x = reorder(Type, {{Gene}}, FUN=median), y = {{Gene}}, fill = Type)) +
           geom_boxplot() +
  labs(x = "TCGA Code") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8.2),
    axis.text.y = element_text(size = 13),
    axis.line = element_line(linewidth = .5),
    legend.position = "none"
    ) 
}
```
```
boxplot_onegene(ASS1)
```
<img width="1869" alt="image" src="https://github.com/kndunlap/tcga-explorer/assets/61035909/4c0efcc0-6b2c-4a71-871e-c2b59524f5f4">





