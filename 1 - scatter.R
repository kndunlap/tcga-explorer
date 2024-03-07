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

# Run

scatter(ASS1, ASL, BRCA)