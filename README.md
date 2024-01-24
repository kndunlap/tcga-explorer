# Introduction
This repository is used to store information about TCGA datasets as well as code on how to mine through them.
In this README, I will share with you six functions that will hopefully help you mine through TCGA data to get some useful results.
These functions use a dataset compiled from the PANCANCER subset in Xenabrowser. If you want to see how I put this dataset together, please check out the code in this repository. This pre-cleaned dataset is stored on the University of Utah Biochemistry Fileserver. It is comprised of 12,724 patients, 20,532 genes, and 2 columns at the beginning indicating the TCGA cancer code (https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations) and the abbreviated patient ID. Below is the code used to very quickly do a final clean on this sheet, and what the output looks like.

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

