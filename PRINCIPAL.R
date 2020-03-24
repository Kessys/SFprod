rm(list = ls())
setwd("C:/Users/kessy/OneDrive/Área de Trabalho/Mestrado UnB/Dissertação/Scripts/SFscript")
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

{
        boxplot(dados$ly~reg, ylim=c(6,16))
        boxplot(mod.h1$fitted.y_without.corretion ~ reg, ylim=c(6,16))
        boxplot(mod.h1$fitted.y_with.corretion ~ reg, ylim=c(6,16))
        
        library(xtable)
        xtable(mod.h2$reg_IV[[1]], type = "latex", digits = 3)
        
        AIC_BIC <- matrix(c(mod.h1$value, mod.h1$AIC, mod.h1$BIC, 
                            mod.e1$value, mod.e1$AIC, mod.e1$BIC, 
                            mod.t1$value, mod.t1$AIC, mod.t1$BIC,
                            mod.h2$value, mod.h2$AIC, mod.h2$BIC,
                            mod.e2$value, mod.e2$AIC, mod.e2$BIC,
                            mod.t2$value, mod.t2$AIC, mod.t2$BIC),
                          nrow = 6, ncol = 3, byrow = T)
        colnames(AIC_BIC) <- c("Log-verossimilhança", "AIC", "BIC")
        rownames(AIC_BIC) <- c("Seminormal", "Exponencial", "Normal truncada",
                               "Seminormal", "Exponencial", "Normal truncada")
        round(AIC_BIC,1)
        xtable(AIC_BIC, type = "latex", digits = 1)
        
        corr <- matrix(c(mod.h1$cor.pearson, mod.e1$cor.pearson, mod.t1$cor.pearson,
                         mod.h2$cor.pearson, mod.e2$cor.pearson, mod.t2$cor.pearson, 
                         diag(mod.h2$`cor pearson IV`)[1], diag(mod.h2$`cor pearson IV`)[2],
                         mod.h1$cor.spearman, mod.e1$cor.spearman, mod.t1$cor.spearman,
                         mod.h2$cor.spearman, mod.e2$cor.spearman, mod.t2$cor.spearman, 
                         diag(mod.h2$`cor spearman IV`)[1], diag(mod.h2$`cor spearman IV`)[2]),
                       nrow = 8, ncol = 2)
        rownames(corr) <- c("Seminormal", "Exponencial", "Normal truncada",
                            "Seminormal", "Exponencial", "Normal truncada",
                            "IV assistec", "IV finan")
        colnames(corr) <- c("Pearson", "Spearman")
        round(corr,4)
        xtable(corr, type = "latex", digits = 4)
        
        vies_REQM <- matrix(c(mod.h1$bias, mod.h1$RMSE, 
                             mod.e1$bias, mod.e1$RMSE, 
                             mod.t1$bias, mod.t1$RMSE,
                             mod.h2$bias, mod.h2$RMSE,
                             mod.e2$bias, mod.e2$RMSE,
                             mod.t2$bias, mod.t2$RMSE),
                           nrow = 6, ncol = 2, byrow = T)
        colnames(vies_REQM) <- c("Viés", "REQM")
        rownames(vies_REQM) <- c("Seminormal", "Exponencial", "Normal truncada",
                                "Seminormal", "Exponencial", "Normal truncada")
        round(vies_REQM,2)
        xtable(vies_REQM, type = "latex", digits = 2)
        
        ajuste <- cbind(AIC_BIC, "correlação" = corr[-c(7,8),-2], vies_REQM)
        round(ajuste,3)
        xtable(ajuste, type = "latex", digits = 3)
        
        Wald_TRV <- matrix(c(mod.h1$Wald, mod.h1$LR,
                             mod.e1$Wald, mod.e1$LR,
                             mod.t1$Wald, mod.t1$LR,
                             mod.h2$Wald, mod.h2$LR,
                             mod.e2$Wald, mod.e2$LR,
                             mod.t2$Wald, mod.t2$LR),
                           nrow = 6, ncol = 4, byrow = T)
        colnames(Wald_TRV) <- c("Wald", "p-valor", "TRV", "p-valor")
        rownames(Wald_TRV) <- c("Seminormal", "Exponencial", "Normal truncada",
                            "Seminormal", "Exponencial", "Normal truncada")
        round(Wald_TRV,0)
        xtable(Wald_TRV, type = "latex", digits = 1)
        
        mean1 <- apply(exp(dados[,c(1:4)]),2,mean)
        quantil1 <- t(apply(exp(dados[,c(1:4)]),2, quantile))
        sd1 <- apply(exp(dados[,c(1:4)]),2,sd)
        CV1 <- (sd1/mean1)*100
        
        mean <- apply(dados[,-c(1:4)],2,mean)
        quantil <- t(apply(dados[,-c(1:4)],2, quantile))
        sd <- apply(dados[,-c(1:4)],2,sd)
        CV <- (sd/mean)*100
        
        library(moments)
        skewness(dados)
        kurtosis(dados)
        
        resumo1 <- cbind('Mínimo'=quantil1[,1],'1º Quartil'=quantil1[,2],'Mediana'=quantil1[,3],'Média'=mean1,'3º Quartil'=quantil1[,4],'Máximo'=quantil1[,5],'Coeficiente de variação'=CV1)
        resumo2 <- cbind('Mínimo'=quantil[,1],'1º Quartil'=quantil[,2],'Mediana'=quantil[,3],'Média'=mean,'3º Quartil'=quantil[,4],'Máximo'=quantil[,5],'Coeficiente de variação'=CV)
        resumo <- rbind(resumo1,resumo2)
        xtable(resumo1, type = "latex", digits = 0)
        xtable(resumo2, type = "latex", digits = 2)
        
        frontier <- ly ~ lxterra + lxtrab + lxcap + reg
        mf <- model.frame(frontier, dados)
        X_fr <- model.matrix(frontier, mf)
        Y_fr <- model.response(mf, "numeric")
        
        # Variaveis endogenas
        End <- cbind(assistec, finan) ~ lxterra + lxtrab + lxcap + reg + social + demo + ambi + gini
        m <- model.frame(End, dados)
        xend <- model.matrix(End, m)
        yend <- model.response(m, "numeric")
        
        cor_IV <- cor(xend[,-1], yend)
        cor_h1 <- cor(xend[,-1], mod.h1$error)
        cor_h2 <- cor(xend[,-1], mod.h2$error)
        corIV <- cbind(cor_IV, cor_h1, cor_h2)
        colnames(corIV) <- c("assistec", "finan", "e.MVIC", "e.MVIL")
        round(corIV,2)
        xtable(corIV, type = "latex", digits = 2)
}
