# 6. gene_correlation_mapper - Correlate over all genes giving one gene.  ---------

gene_correlation_mapper <- function(gene, code) {
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
  
  output <- list()
  for (i in c(2:ncol(alltest))) {            
    output[[i]] <- cor(alltest[,2], alltest[,i])
  }
  output <- output[-1]
  
  listframe <- data.frame(output)
  listframe |>
    pivot_longer(
      cols = 1:ncol(listframe),
      names_to = "gene",
      values_to = "cor"
    ) |>
    arrange(desc(cor)) |>
    print(n = 25)
  
}

# Run

gene_correlation_mapper(SLC7A5, TCGA)