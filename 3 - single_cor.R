# 3. Gives you the correlation value of two genes within your desired cancer type (or all) -----

single_cor <- function(Gene1, Gene2, code) {
  code <- ensym(code)
  all1 <- all |>
    select({{Gene1}}, {{Gene2}}, Type, sample_type) |>
    filter(sample_type != "Solid Tissue Normal") |>
    filter(sample_type != "Recurrent Tumor") |>
    filter(sample_type != "Additional - New Primary") |>
    filter(sample_type != "Additional Metastatic") |>
    filter(if (code != "TCGA") Type == code else TRUE) |>
    summarize(
      corvalue = cor({{Gene1}}, {{Gene2}})
    )
  return(all1)
}

# Run

single_cor(ASS1, ASL, TCGA)