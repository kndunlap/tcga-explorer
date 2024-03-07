# 2. scatter_facet - makes scatter plots of 2 genes for all 36 cancer types --------

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
    print(n = Inf)
  
  print(plot)
}

# Run 

scatter_facet(ASS1, SLC7A5)