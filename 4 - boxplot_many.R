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

# Run

boxplot_many(TCGA, OTC, ASS1, SLC7A5, ARG1, ARG2)