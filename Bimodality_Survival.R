# This code can be used to plot survival curves.
# The data used was downloaded from GDC.
# Datasets include: 
# 1) Liver_survival.csv: liver cancer patients stratified using miR-105 and miR-767 expression 
# 2) Lung_survival.csv: lung cancer patients stratified using miR-105 and miR-767 expression

library(survival)
library(ggplot2)
library(survminer)

# Read datasets
load("liver_lung_survival.RData")

# Fit using Cox proportional hazards
mod_liver = coxph(Surv(Days, Delta) ~ Cluster,
                  data = liver, method = "breslow")

mod_lung = coxph(Surv(Days, Delta) ~ Cluster,
                 data = lung, method = "breslow")

summary(mod_liver)
summary(mod_lung)


liver_plot = survfit(Surv(Days, Delta) ~ Cluster,
                     data = liver)
lung_plot = survfit(Surv(Days, Delta) ~ Cluster,
                    data = lung)


## Kaplan Meier curves

# Plot survival in liver cancer patients
ggsurvplot(liver_plot, data = liver,
           risk.table = TRUE,
           break.time.by = 1000,
           pval = TRUE,
           pval.coord = c(0, 0),
           ggtheme = theme_pubr(),
           risk.table.y.text.col = TRUE,
           risk.table.y.text = FALSE,
           legend.labs = c("Low", "High"),
           xlab = "Time (Days)",
           palette = c("dodgerblue", "red"),
           xlim = c(0,max(liver$Days)+500))


# Plot survival in lung cancer patients
ggsurvplot(lung_plot, data = lung,
           risk.table = TRUE,
           break.time.by = 1000,
           pval = TRUE,
           pval.coord = c(0, 0),
           ggtheme = theme_pubr(),
           risk.table.y.text.col = TRUE,
           risk.table.y.text = FALSE,
           legend.labs = c("Low", "High"),
           xlab = "Time (Days)",
           palette = c("dodgerblue", "red"),
           xlim = c(0,max(lung$Days)+500))

