## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE, message = F, warning = F
)
library(jstable)

## ----eval = F-----------------------------------------------------------------
#  remotes::install_github("jinseob2kim/jstable")
#  library(jstable)

## -----------------------------------------------------------------------------
## Gaussian
glm_gaussian <- glm(mpg ~ cyl + disp, data = mtcars)
glmshow.display(glm_gaussian, decimal = 2)

## Binomial
glm_binomial <- glm(vs ~ cyl + disp, data = mtcars, family = binomial)
glmshow.display(glm_binomial, decimal = 2)

## -----------------------------------------------------------------------------
library(geepack) ## for dietox data
data(dietox)
dietox$Cu <- as.factor(dietox$Cu)
dietox$ddn <- as.numeric(rnorm(nrow(dietox)) > 0)
gee01 <- geeglm(Weight ~ Time + Cu, id = Pig, data = dietox, family = gaussian, corstr = "ex")
geeglm.display(gee01)

gee02 <- geeglm(ddn ~ Time + Cu, id = Pig, data = dietox, family = binomial, corstr = "ex")
geeglm.display(gee02)

## -----------------------------------------------------------------------------
library(lme4)
l1 <- lmer(Weight ~ Time + Cu + (1 | Pig), data = dietox)
lmer.display(l1, ci.ranef = T)

l2 <- glmer(ddn ~ Time + Cu + (1 | Pig), data = dietox, family = "binomial")
lmer.display(l2)

## -----------------------------------------------------------------------------
library(survival)
fit1 <- coxph(Surv(time, status) ~ ph.ecog + age, cluster = inst, lung, model = T) ## model = T: to extract original data
fit2 <- coxph(Surv(time, status) ~ ph.ecog + age + frailty(inst), lung, model = T)
cox2.display(fit1)
cox2.display(fit2)

## -----------------------------------------------------------------------------
library(coxme)
fit <- coxme(Surv(time, status) ~ ph.ecog + age + (1 | inst), lung)
coxme.display(fit)

## -----------------------------------------------------------------------------
library(survey)
data(api)
apistrat$tt <- c(rep(1, 20), rep(0, nrow(apistrat) - 20))
apistrat$tt2 <- factor(c(rep(0, 40), rep(1, nrow(apistrat) - 40)))

dstrat <- svydesign(id = ~1, strata = ~stype, weights = ~pw, data = apistrat, fpc = ~fpc)
ds <- svyglm(api00 ~ ell + meals + tt2, design = dstrat)
ds2 <- svyglm(tt ~ ell + meals + tt2, design = dstrat, family = quasibinomial())
svyregress.display(ds)
svyregress.display(ds2)

## -----------------------------------------------------------------------------
data(pbc, package = "survival")
pbc$sex <- factor(pbc$sex)
pbc$stage <- factor(pbc$stage)
pbc$randomized <- with(pbc, !is.na(trt) & trt > 0)
biasmodel <- glm(randomized ~ age * edema, data = pbc, family = binomial)
pbc$randprob <- fitted(biasmodel)

if (is.null(pbc$albumin)) pbc$albumin <- pbc$alb ## pre2.9.0

dpbc <- svydesign(id = ~1, prob = ~randprob, strata = ~edema, data = subset(pbc, randomized))

model <- svycoxph(Surv(time, status > 0) ~ sex + protime + albumin + stage, design = dpbc)
svycox.display(model)

## -----------------------------------------------------------------------------
library(dplyr)
lung %>%
  mutate(
    status = as.integer(status == 1),
    sex = factor(sex),
    kk = factor(as.integer(pat.karno >= 70)),
    kk1 = factor(as.integer(pat.karno >= 60))
  ) -> lung

# TableSubgroupMultiCox(Surv(time, status) ~ sex, var_subgroups = c("kk", "kk1"), data=lung, line = T)

## Survey data
library(survey)
data.design <- svydesign(id = ~1, data = lung, weights = ~1)
# TableSubgroupMultiCox(Surv(time, status) ~ sex, var_subgroups = c("kk", "kk1"), data = data.design, line = F)

