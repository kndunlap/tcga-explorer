install.packages("survival")
library(survival)

all |>
  select(!3:20531) |>
  ggplot(aes(x = OS.time, y = SELS)) +
  geom_point(aes(color = vital_status))

install.packages(c("lubridate", "ggsurvfit", "gtsummary", "tidycmprsk"))
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(condsurv)

lung1 <- lung |>
  mutate(
    status = recode(status, `1` = 0, `2` = 1)
  )

lung <- lung1

survfit2(Surv(time, status) ~ 1, data = lung) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  )



survfit2(Surv(OS.time, vital_status_num) ~ 1, data = all_vital) |>
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  )


### Function code would start here

kaplan_meier <- function(gene, code) {

library(survival)
library(ggsurvfit)
  
gene <- ensym(gene)
code <- ensym(code)

all_vital <- all |>
  mutate(
    vital_status_num = ifelse(
      vital_status == "Dead", 0, 1)
  ) 

med <- all_vital |>
  filter(if (code != "TCGA") Type == code else TRUE) |>
  summarize(
    median = median({{gene}}, na.rm = TRUE)
  ) |>
  pull(median)

split <- all_vital |>
  filter(if (code != "TCGA") Type == code else TRUE) |>
  mutate(
    class = case_when(
      {{gene}} > med ~ "High",
      {{gene}} < med ~ "Low"
    )
  ) 

survfit2(Surv(OS.time, vital_status_num) ~ class, data = split) |>
  ggsurvfit() +
  labs(
    x = "Days",
    title = gene
  )

}

kaplan_meier(RNF10, BRCA)

survfit2(Surv(OS.time, vital_status_num) ~ class, data = split)

