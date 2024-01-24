# tcga-explorer
This repository is used to store information about TCGA datasets as well as code on how to mine through them.
In this README, I will share with you six functions that will hopefully help you mine through TCGA data to get some useful results.
These functions use a dataset compiled from the PANCANCER subset in Xenabrowser. This pre-cleaned dataset is stored on the University of Utah Biochemistry Fileserver. It is comprised of 12,724 patients, 20,532 genes, and 2 columns at the beginning indicating the TCGA cancer code and the abbreviated patient ID. Below is the code used to very quickly do a final clean on this sheet, and what the output looks like.

```
all <- read.csv("TCGA_all_tidy.csv")
all <- as.tibble(all)
all <- all |> select(!1)
```

<img width="897" alt="image" src="https://github.com/kndunlap/tcga-explorer/assets/61035909/526fca6a-9f26-416a-8795-32402f90eb22">

