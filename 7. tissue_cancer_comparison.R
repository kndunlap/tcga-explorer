# 7. tissue_cancer_comparison - Create a Boxplot of Solid Tissue vs. Cancer in Your Type --------

tissue_cancer_comparison <- function(Gene, code){
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

# Run

tissue_cancer_comparison(BRCA1, BRCA)