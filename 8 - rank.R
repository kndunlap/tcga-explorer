# 8. rank - Ranks cancer types by expression of given gene --------------

rank <- function(gene) {
  gene <- ensym(gene)
  all |>
    select({{gene}}, Type) |>
    group_by(Type) |>
    summarize(
      mean = mean({{gene}}, na.rm = TRUE) 
    ) |>
    arrange(desc(mean)) |>
    print(n = Inf)
}

# Rank

rank(OTC)