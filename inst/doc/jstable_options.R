## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE, message = F, warning = F
)
library(jstable)
library(survival)
library(dplyr)

## -----------------------------------------------------------------------------
lung %>%
  mutate(
    status = as.integer(status == 1),
    sex = factor(sex),
    kk = factor(as.integer(pat.karno >= 70)),
    kk1 = factor(as.integer(pat.karno >= 60)),
    ph.ecog = factor(ph.ecog)
  ) -> lung
lung.label <- mk.lev(lung)
lung.label <- lung.label %>%
  mutate(
    val_label = case_when(
      variable == "sex" & level == "1" ~ "Male",
      variable == "sex" & level == "2" ~ "Female",
      variable == "kk" & level == "0" ~ "No",
      variable == "kk" & level == "1" ~ "Yes",
      variable == "kk1" & level == "0" ~ "No",
      variable == "kk1" & level == "1" ~ "Yes",
      TRUE ~ val_label
    )
  )

## -----------------------------------------------------------------------------
TableSubgroupMultiCox(Surv(time, status) ~ sex, var_subgroups = c("kk", "kk1"), data = lung, time_eventrate = 100, line = TRUE, cluster = "inst", strata = "inst", weights = "age", event = FALSE, count_by = "sex", labeldata = lung.label)

## -----------------------------------------------------------------------------
TableSubgroupMultiCox(Surv(time, status) ~ sex, var_subgroups = c("kk", "kk1"), data = lung, time_eventrate = 100, line = TRUE, cluster = "inst", strata = "inst", weights = "age", event = TRUE, count_by = "sex", labeldata = lung.label)

## -----------------------------------------------------------------------------
TableSubgroupMultiCox(Surv(time, status) ~ sex, var_subgroups = c("kk", "kk1"), data = lung, time_eventrate = 100, line = TRUE, cluster = "inst", strata = "inst", weights = "age", event = TRUE, count_by = NULL, labeldata = lung.label)

## -----------------------------------------------------------------------------
CreateTableOneJS(vars = names(lung), strata = "ph.ecog", data = lung, showAllLevels = F, labeldata = lung.label, Labels = T, pairwise = T)

## -----------------------------------------------------------------------------
CreateTableOneJS(vars = names(lung), strata = "ph.ecog", data = lung, showAllLevels = F, labeldata = lung.label, Labels = T, pairwise = T, pairwise.showtest = T)

## -----------------------------------------------------------------------------
library(survey)
data(nhanes)
nhanes$SDMVPSU <- as.factor(nhanes$SDMVPSU)
nhanes$race <- as.factor(nhanes$race)
nhanes$RIAGENDR <- as.factor(nhanes$RIAGENDR)
a.label <- mk.lev(nhanes)
a.label <- a.label %>%
  dplyr::mutate(val_label = case_when(
    variable == "race" & level == "1" ~ "White",
    variable == "race" & level == "2" ~ "Black",
    variable == "race" & level == "3" ~ "Hispanic",
    variable == "race" & level == "4" ~ "Asian",
    TRUE ~ val_label
  ))
nhanesSvy <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC2YR, nest = TRUE, data = nhanes)

svyCreateTableOneJS(
  vars = c("HI_CHOL", "race", "agecat", "RIAGENDR"),
  strata = "race", data = nhanesSvy, factorVars = c("HI_CHOL", "race", "RIAGENDR"), labeldata = a.label, Labels = T, pairwise = T
)

## -----------------------------------------------------------------------------
svyCreateTableOneJS(
  vars = c("HI_CHOL", "race", "agecat", "RIAGENDR"),
  strata = "race", data = nhanesSvy, factorVars = c("HI_CHOL", "race", "RIAGENDR"), labeldata = a.label, Labels = T, pairwise = T, pairwise.showtest = T
)

