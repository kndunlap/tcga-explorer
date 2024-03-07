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
    labs(x = "TCGA Code") +
    theme_minimal() +
    theme(axis.line = element_line(linewidth = .5), 
          axis.text.x = element_text(size = 9, angle = 45, hjust = 1, face = "bold"),
          legend.position = "none")
}

# Run

boxplot_onegene(ASS1)