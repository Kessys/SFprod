rm(list = ls())
setwd("C:/Users/kessy/Dropbox/Artigo Kessys/codes")
library(Matrix); library(MASS); library(car); library(numDeriv)
load("dados.RData")
reg <- factor(dados$regiao, levels = c(5, 1, 2, 3, 4))
source("SFprod.R")

mod.h0 <- SF.prod(model = "half",
                  fr.form  = ly ~ lxterra + lxtrab + lxcap + reg,
                  s2u.form = ~ assistec + finan + gini,
                  s2w.form = ~ 1,
                  data = dados)

mod.h1 <- SF.prod(model = "half1s",
                  fr.form  = ly ~ lxterra + lxtrab + lxcap + reg,
                  s2u.form = ~ assistec + finan + gini,
                  s2w.form = ~ 1,
                  end.form = cbind(assistec, finan) ~ lxterra + lxtrab + lxcap + reg + social + demo + ambi + gini,
                  data = dados)

mod.h2 <- SF.prod(model = "half2s",
                  fr.form  = ly ~ lxterra + lxtrab + lxcap + reg,
                  s2u.form = ~ assistec + finan + gini,
                  s2w.form = ~ 1,
                  end.form = cbind(assistec, finan) ~ lxterra + lxtrab + lxcap + reg + social + demo + ambi + gini,
                  data = dados)

mod.e0 <- SF.prod(model = "exp",
                  fr.form  = ly ~ lxterra + lxtrab + lxcap + reg,
                  s2u.form = ~ assistec + finan + gini,
                  s2w.form = ~ 1,
                  data = dados)

mod.e1 <- SF.prod(model = "exp1s",
                  fr.form  = ly ~ lxterra + lxtrab + lxcap + reg,
                  s2u.form = ~ assistec + finan + gini,
                  s2w.form = ~ 1,
                  end.form = cbind(assistec, finan) ~ lxterra + lxtrab + lxcap + reg + social + demo + ambi + gini,
                  data = dados) #53280 iterações

mod.e2 <- SF.prod(model = "exp2s",
                  fr.form  = ly ~ lxterra + lxtrab + lxcap + reg,
                  s2u.form = ~ assistec + finan + gini,
                  s2w.form = ~ 1,
                  end.form = cbind(assistec, finan) ~ lxterra + lxtrab + lxcap + reg + social + demo + ambi + gini,
                  data = dados) #4514 iterações

mod.t0 <- SF.prod(model = "trunc",
                  fr.form  = ly ~ lxterra + lxtrab + lxcap + reg,
                  mu.form = ~ assistec + finan + gini,
                  data = dados)

mod.t1 <- SF.prod(model = "trunc1s",
                  fr.form  = ly ~ lxterra + lxtrab + lxcap + reg,
                  mu.form  = ~ assistec + finan + gini,
                  end.form = cbind(assistec, finan) ~ lxterra + lxtrab + lxcap + reg + social + demo + ambi + gini,
                  data = dados)

mod.t2 <- SF.prod(model = "trunc2s",
                  fr.form  = ly ~ lxterra + lxtrab + lxcap + reg,
                  mu.form  = ~ assistec + finan + gini,
                  end.form = cbind(assistec, finan) ~ lxterra + lxtrab + lxcap + reg + social + demo + ambi + gini,
                  data = dados)
