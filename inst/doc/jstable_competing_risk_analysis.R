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
data <- mgus2
data$etime <- with(data, ifelse(pstat == 0, futime, ptime))
data$event <- with(data, ifelse(pstat == 0, 2 * death, 1))
data$event <- factor(data$event, 0:2, labels = c("censor", "pcm", "death"))
data$age65 <- with(data, ifelse(age > 65, 1, 0))
data$age65 <- factor(data$age65)
pdata <- survival::finegray(survival::Surv(etime, event) ~ ., data = data)
TableSubgroupMultiCox(formula = Surv(fgstart, fgstop, fgstatus) ~ sex, data = pdata, var_cov = "age", weights = "fgwt", var_subgroups = c("age65"))

## -----------------------------------------------------------------------------
fgfit <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age + sex,
  weight = fgwt, data = pdata, model = T
)
cox2.display(fgfit)

