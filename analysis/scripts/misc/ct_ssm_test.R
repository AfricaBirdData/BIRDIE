library(BIRDIE)
library(tidyverse)

rm(list = ls())

counts <- barberspan

ssmcounts <- prepCtSsmData(counts, species = 4)

fit <- fitCwacCtSsm(ssmcounts, mod_file = "analysis/models/cwac_ct_ssm_dyn.jags",
                    param = c("beta", "lambda", "sig.B", "sig.zeta", "sig.eps",
                              "sig.alpha", "sig.e", "sig.o", "mu"),
                    jags_control = list(ncores = 3))

summary(fit)
fit

out <- plotCtSsm(fit, ssmcounts)

plot(out$plot)

out$data
