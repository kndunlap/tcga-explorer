### Import packages and files and kidney patients

library(tidyverse)

kidney <- read_tsv("gtex_RSEM_Hugo_norm_count.tsv")

patients <- c("sample",
  "GTEX-11OF3-1326-SM-5N9FJ",
  "GTEX-11PRG-2226-SM-5GU5R",
  "GTEX-12696-0926-SM-5FQTV",
  "GTEX-12WSG-0826-SM-5EQ5A",
  "GTEX-13112-2126-SM-5GCO4",
  "GTEX-1399S-0526-SM-5IJG8",
  "GTEX-13NYB-1726-SM-5N9G2",
  "GTEX-13O1R-2526-SM-5N9FW",
  "GTEX-13OW6-1826-SM-5N9F9",
  "GTEX-13RTJ-2226-SM-5S2Q1",
  "GTEX-145MN-0326-SM-5QGQI",
  "GTEX-147F4-2626-SM-5Q5CS",
  "GTEX-1497J-0826-SM-5NQAJ",
  "GTEX-N7MS-1626-SM-3LK5F",
  "GTEX-NPJ8-2226-SM-3TW8D",
  "GTEX-QLQW-1626-SM-4R1K1",
  "GTEX-T5JC-1526-SM-4DM68",
  "GTEX-T6MN-1826-SM-5CHQB",
  "GTEX-WI4N-2026-SM-4OOS7",
  "GTEX-XPVG-0526-SM-4B65N",
  "GTEX-Y5V6-2026-SM-5IFHO",
  "GTEX-ZC5H-1726-SM-5HL7X",
  "GTEX-ZDXO-0226-SM-4WKH7",
  "GTEX-ZE9C-1426-SM-4WKGM",
  "GTEX-ZLFU-0926-SM-5P9F8",
  "GTEX-ZYFD-1526-SM-5NQ7T",
  "GTEX-ZYFG-1626-SM-5GZYY",
  "GTEX-ZYT6-2226-SM-5GIC9"
)

### Collect only SLC genes and make tidy table
kidneyOnly_SLC <- kidney |>
  select(patients) |>
  filter(startsWith(sample, "SLC"))

cleanKidney <- kidneyOnly_SLC |>
  pivot_longer(
    cols = starts_with("GTEX-"),
    names_to = "patient", 
    values_to = "log2exp"
  ) |>
  pivot_wider(
    names_from = sample,
    values_from = log2exp
  ) |>
  relocate(SLC5A2) |>
  select(!patient)


### Filter out non-expressed SLCs - cor is not necessary but it can be done at this point

cleanKidney_expr <- kidneyOnly_SLC |>
  pivot_longer(
    cols = starts_with("GTEX-"),
    names_to = "patient", 
    values_to = "log2exp"
  ) |>
  pivot_wider(
    names_from = sample,
    values_from = log2exp
  ) |>
  relocate(SLC6A19) |>
  select(!patient) |>
  select(where(~ mean(.x, na.rm = TRUE) >= 1))


output <- list()
for (i in c(2:ncol(cleanKidney_expr))) {            
  output[[i]] <- cor(cleanKidney_expr[,1], cleanKidney_expr[,i])
}
output <- output[-1]

listframe <- data.frame(output)
names <- listframe |>
  pivot_longer(
    cols = 1:ncol(listframe),
    names_to = "gene",
    values_to = "cor"
  ) |>
  arrange(desc(cor)) |>
  print(n = 25) 


ALLSLCgenes <- names$gene


### Took above list and ran through metascape. Starting here is a list of potentially neutral transporters

slc_vector <- c("SLC6A19", "SLC5A2", "SLC7A8", "SLC7A9", "SLC16A9", "SLC7A7", "SLC47A1", "SLC6A17", "SLC6A13", "SLC16A10", "SLC4A9",
                "SLC3A1", "SLC38A3", "SLC43A2", "SLC3A2", "SLC44A3", "SLC38A7", "SLC46A3", "SLC25A15", "SLC38A4", "SLC1A5", "SLC38A6",
                "SLC7A14", "SLC38A9", "SLC22A18", "SLC7A2", "SLC41A3", "SLC38A2", "SLC15A4", "SLC36A4", "SLC7A10", "SLC7A6", "SLC43A1",
                "SLC7A5", "SLC7A13", "SLC38A11", "SLC1A4", "SLC7A3", "SLC38A5", "SLC14A1", "SLC14A2", "SLC6A7", "SLC44A5", "SLC22A31",
                "SLC7A4", "SLC7A1", "SLC6A15", "SLC25A2")

cleanKidney_neutral <- cleanKidney |>
  select(slc_vector)

output <- list()
for (i in c(1:ncol(cleanKidney_neutral))) {            
  output[[i]] <- cor(cleanKidney_neutral[,1], cleanKidney_neutral[,i])
}

listframe <- data.frame(output)
neutralNames <- listframe |>
  pivot_longer(
    cols = 1:ncol(listframe),
    names_to = "gene",
    values_to = "cor_with_SLC5A2"
  ) |>
  arrange(desc(cor_with_SLC5A2)) |>
  print(n = 50) 

## Join output of mean expression with neutralNames

neutralMeans <- cleanKidney_neutral |>
  summarise(across(everything(), mean)) |>
  pivot_longer(
    cols = 1:ncol(cleanKidney_neutral),
    names_to = "gene",
    values_to = "expr"
  )

final <- neutralNames |>
  left_join(neutralMeans)

annotation <- read.csv("SLCannotation.csv")

final_with_annotation <- final |>
  left_join(annotation)


### Smaller list of potentially neutral

annotation_neutral <- read.csv("SLCannotation_neutral.csv")

final_with_annotation_neutral <- final |>
  inner_join(annotation_neutral)

write_csv(final_with_annotation_neutral, "neutral_AAs.csv")


### Looking at CCLE now.
