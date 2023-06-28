# Conventional stochastic production frontier

SF.half <- function(fr.form, s2u.form, s2w.form, data = sys.parent()) {
        {
                # Frontier
                frontier <- fr.form
                mf <- model.frame(frontier, data)
                X_fr <- model.matrix(frontier, mf)
                Y_fr <- model.response(mf, "numeric")
                
                # Sigma_u and Sigma_w
                Sigma_u <- s2u.form
                sig_u <- model.frame(Sigma_u, data)
                sigma_u <- model.matrix(Sigma_u, sig_u)
                
                Sigma_w <- s2w.form
                sig_w <- model.frame(Sigma_w, data)
                sigma_w <- model.matrix(Sigma_w, sig_w)
                
                # Initial values
                fr <- lm(frontier, data = data); names(fr$coefficients)[1] <- c("constant")
                u.p <- c("constant" = 1, rep(0, ncol(sigma_u)-1)); if(length(names(u.p)) > 1){names(u.p)[2:length(u.p)] <- colnames(sigma_u)[2:length(u.p)]}
                w.p <- c("constant" = 1, rep(0, ncol(sigma_w)-1)); if(length(names(w.p)) > 1){names(w.p)[2:length(w.p)] <- colnames(sigma_w)[2:length(w.p)]}
                
                z <- c(coef(fr), u.p, w.p)
                names(z) <- make.names(names(z), unique=T)
                
                kk <- length(z)
                k1 <- length(coef(fr))
                k2 <- k1 + length(u.p)
        }
        ll <- function(z) {
                x_p <- z[1:k1]
                u_p <- z[c(k1 + 1):k2]
                w_p <- z[c(k2 + 1):kk]
                
                Y <- X_fr %*% x_p
                U <- sigma_u %*% u_p
                W <- sigma_w %*% w_p
                
                Sigma_u <- exp(U)
                Sigma_w <- exp(W)
                lambda <- sqrt(Sigma_u / Sigma_w)
                sigma <- Sigma_u + Sigma_w
                e <- Y_fr - Y
                
                zz <- - e * lambda / sqrt(sigma)
                pz <- pmax(pnorm(zz), 9.88131291682493e-324)
                ly <- 0.5 * log(2 / pi) - 0.5 * log(sigma) + log(pz) - 0.5 * e^2 / sigma
                val <- sum(ly)
                
                return(-val)
        }
        G <- function(z) {
                x_p <- z[1:k1]
                u_p <- z[c(k1 + 1):k2]
                w_p <- z[c(k2 + 1):kk]
                
                Y <- X_fr %*% x_p
                U <- sigma_u %*% u_p
                W <- sigma_w %*% w_p
                
                Sigma_u <- exp(U)
                Sigma_w <- exp(W)
                lambda <- sqrt(Sigma_u / Sigma_w)
                sigma <- Sigma_u + Sigma_w
                e <- Y_fr - Y
                
                zz <- - e * lambda / sqrt(sigma)
                dz <- pmax(dnorm(zz), 9.88131291682493e-324)
                pz <- pmax(pnorm(zz), 9.88131291682493e-324)
                fdp <- dz; cdf <- pz; fdp_cdf <- fdp / cdf
                valor <- e / sigma + lambda / sqrt(sigma) * fdp_cdf
                
                g.b <- t(X_fr) %*% valor
                g.u <- t(sigma_u) %*% ((0.5 / sigma * (e^2 / sigma - e / (lambda * sqrt(sigma)) * fdp_cdf - 1)) * Sigma_u)
                g.w <- t(sigma_w) %*% ((0.5 / sigma * (e^2 / sigma + e * lambda * (2 + lambda ^ 2) * fdp_cdf / sqrt(sigma) - 1)) * Sigma_w)
                
                g <- c(g.b, g.u, g.w)
                
                return(-g)
        }
        est <- optim(par = z, fn = ll, gr = G, hessian = TRUE, method = "BFGS",
                     control = list(fnscale = 1, trace = TRUE, REPORT = 1, maxit = 150000))
        #round(data.frame("numerical" = grad(ll,est$par), "analytical" = G(est$par)), 6)
        {
                ep <- sqrt(diag(ginv(est$hessian)))
                estZ <- lapply(length(est$par), function(x) est$par/ep)
                p_valor <- lapply(length(est$par), function(x) (1 - pnorm(abs(est$par/ep)))*2)
                ics <- lapply(length(est$par), function(x) cbind(est$par - qnorm(.975) * ep, est$par + qnorm(.975) * ep))
                result <- data.frame("Coefficient" = cbind(est$par), "Std.Err" = ep, "z" = estZ[[1]], 
                                     "P-value" = p_valor[[1]], "Lower limit" = ics[[1]][,1], "Upper limit" = ics[[1]][,2])
                
                lns2u <- est$par[c(k1 + 1):k2]
                lns2w <- est$par[c(k2 + 1):kk]
                if (length(lns2u) == 1 & length(lns2w) == 1){
                        names(lns2u) <- c('lns2u')
                        cov.u <- ginv(est$hessian)[c(k1 + 1):k2,c(k1 + 1):k2]
                        S2U <- deltaMethod(lns2u, "exp(lns2u)", vcov = cov.u)
                        zS2U <- S2U$Estimate/S2U$SE
                        pS2U <- (1 - pnorm(abs(zS2U)))*2
                        icS2U <- lapply(length(S2U$Estimate), function(x) cbind(S2U$Estimate - qnorm(.975) * S2U$SE, S2U$Estimate + qnorm(.975) * S2U$SE))
                        S2U0 <- data.frame(S2U$Estimate, S2U$SE, zS2U, pS2U, icS2U[[1]][,1], icS2U[[1]][,2])
                        names(S2U0) <- names(result)
                        rownames(S2U0) <- "s2u"
                        
                        names(lns2w) <- c('lns2w')
                        cov.w <- ginv(est$hessian)[c(k2 + 1):kk,c(k2 + 1):kk]
                        S2W <- deltaMethod(lns2w, "exp(lns2w)", vcov = cov.w)
                        zS2W <- S2W$Estimate/S2W$SE
                        pS2W <- (1 - pnorm(abs(zS2W)))*2
                        icS2W <- lapply(length(S2W$Estimate), function(x) cbind(S2W$Estimate - qnorm(.975) * S2W$SE, S2W$Estimate + qnorm(.975) * S2W$SE))
                        S2W0 <- data.frame(S2W$Estimate, S2W$SE, zS2W, pS2W, icS2W[[1]][,1], icS2W[[1]][,2])
                        names(S2W0) <- names(result)
                        rownames(S2W0) <- "s2w"
                        
                        s20 <- c(S2U$Estimate, S2W$Estimate)
                        names(s20) <- c('s2u','s2w')
                        cov.s <- ginv(est$hessian)[c(k1 + 1):kk,c(k1 + 1):kk]
                        S2 <- deltaMethod(s20, "s2u + s2w", vcov = cov.s)
                        zS2 <- S2$Estimate/S2$SE
                        pS2 <- (1 - pnorm(abs(zS2)))*2
                        icS2 <- lapply(length(S2$Estimate), function(x) cbind(S2$Estimate - qnorm(.975) * S2$SE, S2$Estimate + qnorm(.975) * S2$SE))
                        S20 <- data.frame(S2$Estimate, S2$SE, zS2, pS2, icS2[[1]][,1], icS2[[1]][,2])
                        names(S20) <- names(result)
                        rownames(S20) <- "s2"
                        
                        Lambda <- deltaMethod(s20, "sqrt(s2u/s2w)", vcov = cov.s)
                        zLam <- Lambda$Estimate/Lambda$SE
                        pLam <- (1 - pnorm(abs(zLam)))*2
                        icLam <- lapply(length(Lambda$Estimate), function(x) cbind(Lambda$Estimate - qnorm(.975) * Lambda$SE, Lambda$Estimate + qnorm(.975) * Lambda$SE))
                        LAM <- data.frame(Lambda$Estimate, Lambda$SE, zLam, pLam, icLam[[1]][,1], icLam[[1]][,2])
                        names(LAM) <- names(result)
                        rownames(LAM) <- "lambda"
                        
                        RESP <- data.frame(rbind(result, S2U0, S2W0, S20, LAM))
                } else if (length(lns2u) == 1){
                        lns2u <- est$par[c(k1 + 1):k2]
                        names(lns2u) <- c('lns2u')
                        cov.u <- ginv(est$hessian)[c(k1 + 1):k2,c(k1 + 1):k2]
                        S2U <- deltaMethod(lns2u, "exp(lns2u)", vcov = cov.u)
                        zS2U <- S2U$Estimate/S2U$SE
                        pS2U <- (1 - pnorm(abs(zS2U)))*2
                        icS2U <- lapply(length(S2U$Estimate), function(x) cbind(S2U$Estimate - qnorm(.975) * S2U$SE, S2U$Estimate + qnorm(.975) * S2U$SE))
                        S2U0 <- data.frame(S2U$Estimate, S2U$SE, zS2U, pS2U, icS2U[[1]][,1], icS2U[[1]][,2])
                        names(S2U0) <- names(result)
                        rownames(S2U0) <- "s2u"
                        
                        RESP <- data.frame(rbind(result, S2U0))
                } else if (length(lns2w) == 1){
                        lns2w <- est$par[c(k2 + 1):kk]
                        names(lns2w) <- c('lns2w')
                        cov.w <- ginv(est$hessian)[c(k2 + 1):kk,c(k2 + 1):kk]
                        S2W <- deltaMethod(lns2w, "exp(lns2w)", vcov = cov.w)
                        zS2W <- S2W$Estimate/S2W$SE
                        pS2W <- (1 - pnorm(abs(zS2W)))*2
                        icS2W <- lapply(length(S2W$Estimate), function(x) cbind(S2W$Estimate - qnorm(.975) * S2W$SE, S2W$Estimate + qnorm(.975) * S2W$SE))
                        S2W0 <- data.frame(S2W$Estimate, S2W$SE, zS2W, pS2W, icS2W[[1]][,1], icS2W[[1]][,2])
                        names(S2W0) <- names(result)
                        rownames(S2W0) <- "s2w"
                        
                        RESP <- data.frame(rbind(result, S2W0))
                } else {
                        RESP <- result
                }
                
                Yest <- X_fr %*% est$par[1:k1]
                s2u <- exp(sigma_u %*% est$par[c(k1 + 1):k2])
                s2w <- exp(sigma_w %*% est$par[c(k2 + 1):kk])
                erro <- Y_fr - Yest
                s2 <- s2u + s2w
                
                # Efficiency estimate by E(exp(-u)|e)
                mu.mod <- - erro * s2u / s2 # media
                s.mod <- sqrt(s2w * s2u/s2) # desvio padrao
                uf <- mu.mod + s.mod * ((dnorm(-mu.mod / s.mod)) / (pnorm(mu.mod / s.mod)))
                ef <- ((pnorm(-s.mod + mu.mod / s.mod)) / (pnorm(mu.mod / s.mod))) * exp(- mu.mod + 0.5 * s.mod^2)
                
                bias <- mean(Yest) - mean(Y_fr)
                RMSE <- sqrt(var(Yest) + bias^2)
                pearson <- cor(Yest, Y_fr, method = "pearson")
                spearman <- cor(Yest, Y_fr, method = "spearman")
                
                n <- length(Y_fr)
                K <- length(est$par)
                lnL <- -est$value
                AIC <- - 2 * lnL + 2 * K
                BIC <- - 2 * lnL + log(n) * K
        }
        lista <- list("efficiency" = ef, "error" = erro, "fitted.y" = Yest,
                      "table" = RESP, "summary.ef" = summary(ef), "sd.ef" = sd(ef),
                      "cor.pearson" = c(pearson), "cor.spearman" = c(spearman),
                      "value" = -est$value, "AIC" = AIC, "BIC" = BIC, 
                      "bias" = c(bias), "RMSE" = c(RMSE), 
                      "sample size" = n, "estimated parameters" = K)
        return(lista)
}

SF.exp <- function(fr.form, s2u.form, s2w.form, data = sys.parent()) {
        {
                # Frontier
                frontier <- fr.form
                mf <- model.frame(frontier, data)
                X_fr <- model.matrix(frontier, mf)
                Y_fr <- model.response(mf, "numeric")
                
                # Sigma_u and Sigma_w
                Sigma_u <- s2u.form
                sig_u <- model.frame(Sigma_u, data)
                sigma_u <- model.matrix(Sigma_u, sig_u)
                
                Sigma_w <- s2w.form
                sig_w <- model.frame(Sigma_w, data)
                sigma_w <- model.matrix(Sigma_w, sig_w)
                
                # Initial values
                fr <- lm(frontier, data = data); names(fr$coefficients)[1] <- c("constant")
                u.p <- c("constant" = 1, rep(0, ncol(sigma_u)-1)); if(length(names(u.p)) > 1){names(u.p)[2:length(u.p)] <- colnames(sigma_u)[2:length(u.p)]}
                w.p <- c("constant" = 1, rep(0, ncol(sigma_w)-1)); if(length(names(w.p)) > 1){names(w.p)[2:length(w.p)] <- colnames(sigma_w)[2:length(w.p)]}
                
                z <- c(coef(fr), u.p, w.p)
                names(z) <- make.names(names(z), unique=T)
                
                kk <- length(z)
                k1 <- length(coef(fr))
                k2 <- k1 + length(u.p)
        }
        ll <- function(z) {
                x_p <- z[1:k1]
                u_p <- z[c(k1 + 1):k2]
                w_p <- z[c(k2 + 1):kk]
                
                Y <- X_fr %*% x_p
                U <- sigma_u %*% u_p
                W <- sigma_w %*% w_p
                
                Sigma_u <- exp(U)
                Sigma_w <- exp(W)
                lambda <- sqrt(Sigma_u / Sigma_w)
                sigma <- Sigma_u + Sigma_w
                e <- Y_fr - Y
                
                zz <- (-e - Sigma_w / sqrt(Sigma_u)) / sqrt(Sigma_w)
                pz <- pnorm(zz, log = TRUE)
                ly <- -0.5 * log(Sigma_u) + 0.5 * Sigma_w / Sigma_u + pz + e / sqrt(Sigma_u)
                val <- sum(ly)
                
                return(-val)
        }
        G <- function(z) {
                x_p <- z[1:k1]
                u_p <- z[c(k1 + 1):k2]
                w_p <- z[c(k2 + 1):kk]
                
                Y <- X_fr %*% x_p
                U <- sigma_u %*% u_p
                W <- sigma_w %*% w_p
                
                Sigma_u <- exp(U)
                Sigma_w <- exp(W)
                lambda <- sqrt(Sigma_u / Sigma_w)
                sigma <- Sigma_u + Sigma_w
                e <- Y_fr - Y
                
                zz <- (-e - Sigma_w / sqrt(Sigma_u)) / sqrt(Sigma_w)
                dz <- pmax(dnorm(zz), 9.88131291682493e-324)
                pz <- pmax(pnorm(zz), 9.88131291682493e-324)
                fdp <- dz; cdf <- pz; fdp_cdf <- fdp / cdf
                valor <- fdp_cdf / sqrt(Sigma_w) - 1 / sqrt(Sigma_u)
                
                g.b <- t(X_fr) %*% valor
                g.u <- t(sigma_u) %*% as.vector((0.5 / Sigma_u * (fdp_cdf * sqrt(Sigma_w / Sigma_u) - Sigma_w / Sigma_u - e / sqrt(Sigma_u) - 1)) * Sigma_u)
                g.w <- t(sigma_w) %*% as.vector((0.5 * (1 / Sigma_u + fdp_cdf * (e / Sigma_w - 1 / sqrt(Sigma_u)) / sqrt(Sigma_w))) * Sigma_w)
                
                g <- c(g.b, g.u, g.w)
                
                return(-g)
        }
        est <- optim(par = z, fn = ll, gr = G, hessian = TRUE, method = "Nelder-Mead",
                     control = list(fnscale = 1, trace = TRUE, REPORT = 1, maxit = 150000))
        #round(data.frame("numerical" = grad(ll,est$par), "analytical" = G(est$par)), 6)
        {
                ep <- sqrt(diag(ginv(est$hessian)))
                estZ <- lapply(length(est$par), function(x) est$par/ep)
                p_valor <- lapply(length(est$par), function(x) (1 - pnorm(abs(est$par/ep)))*2)
                ics <- lapply(length(est$par), function(x) cbind(est$par - qnorm(.975) * ep, est$par + qnorm(.975) * ep))
                result <- data.frame("Coefficient" = cbind(est$par), "Std.Err" = ep, "z" = estZ[[1]], 
                                     "P-value" = p_valor[[1]], "Lower limit" = ics[[1]][,1], "Upper limit" = ics[[1]][,2])
                
                lns2u <- est$par[c(k1 + 1):k2]
                lns2w <- est$par[c(k2 + 1):kk]
                if (length(lns2u) == 1 & length(lns2w) == 1){
                        names(lns2u) <- c('lns2u')
                        cov.u <- ginv(est$hessian)[c(k1 + 1):k2,c(k1 + 1):k2]
                        S2U <- deltaMethod(lns2u, "exp(lns2u)", vcov = cov.u)
                        zS2U <- S2U$Estimate/S2U$SE
                        pS2U <- (1 - pnorm(abs(zS2U)))*2
                        icS2U <- lapply(length(S2U$Estimate), function(x) cbind(S2U$Estimate - qnorm(.975) * S2U$SE, S2U$Estimate + qnorm(.975) * S2U$SE))
                        S2U0 <- data.frame(S2U$Estimate, S2U$SE, zS2U, pS2U, icS2U[[1]][,1], icS2U[[1]][,2])
                        names(S2U0) <- names(result)
                        rownames(S2U0) <- "s2u"
                        
                        names(lns2w) <- c('lns2w')
                        cov.w <- ginv(est$hessian)[c(k2 + 1):kk,c(k2 + 1):kk]
                        S2W <- deltaMethod(lns2w, "exp(lns2w)", vcov = cov.w)
                        zS2W <- S2W$Estimate/S2W$SE
                        pS2W <- (1 - pnorm(abs(zS2W)))*2
                        icS2W <- lapply(length(S2W$Estimate), function(x) cbind(S2W$Estimate - qnorm(.975) * S2W$SE, S2W$Estimate + qnorm(.975) * S2W$SE))
                        S2W0 <- data.frame(S2W$Estimate, S2W$SE, zS2W, pS2W, icS2W[[1]][,1], icS2W[[1]][,2])
                        names(S2W0) <- names(result)
                        rownames(S2W0) <- "s2w"
                        
                        s20 <- c(S2U$Estimate, S2W$Estimate)
                        names(s20) <- c('s2u','s2w')
                        cov.s <- ginv(est$hessian)[c(k1 + 1):kk,c(k1 + 1):kk]
                        S2 <- deltaMethod(s20, "s2u + s2w", vcov = cov.s)
                        zS2 <- S2$Estimate/S2$SE
                        pS2 <- (1 - pnorm(abs(zS2)))*2
                        icS2 <- lapply(length(S2$Estimate), function(x) cbind(S2$Estimate - qnorm(.975) * S2$SE, S2$Estimate + qnorm(.975) * S2$SE))
                        S20 <- data.frame(S2$Estimate, S2$SE, zS2, pS2, icS2[[1]][,1], icS2[[1]][,2])
                        names(S20) <- names(result)
                        rownames(S20) <- "s2"
                        
                        Lambda <- deltaMethod(s20, "sqrt(s2u/s2w)", vcov = cov.s)
                        zLam <- Lambda$Estimate/Lambda$SE
                        pLam <- (1 - pnorm(abs(zLam)))*2
                        icLam <- lapply(length(Lambda$Estimate), function(x) cbind(Lambda$Estimate - qnorm(.975) * Lambda$SE, Lambda$Estimate + qnorm(.975) * Lambda$SE))
                        LAM <- data.frame(Lambda$Estimate, Lambda$SE, zLam, pLam, icLam[[1]][,1], icLam[[1]][,2])
                        names(LAM) <- names(result)
                        rownames(LAM) <- "lambda"
                        
                        RESP <- data.frame(rbind(result, S2U0, S2W0, S20, LAM))
                } else if (length(lns2u) == 1){
                        lns2u <- est$par[c(k1 + 1):k2]
                        names(lns2u) <- c('lns2u')
                        cov.u <- ginv(est$hessian)[c(k1 + 1):k2,c(k1 + 1):k2]
                        S2U <- deltaMethod(lns2u, "exp(lns2u)", vcov = cov.u)
                        zS2U <- S2U$Estimate/S2U$SE
                        pS2U <- (1 - pnorm(abs(zS2U)))*2
                        icS2U <- lapply(length(S2U$Estimate), function(x) cbind(S2U$Estimate - qnorm(.975) * S2U$SE, S2U$Estimate + qnorm(.975) * S2U$SE))
                        S2U0 <- data.frame(S2U$Estimate, S2U$SE, zS2U, pS2U, icS2U[[1]][,1], icS2U[[1]][,2])
                        names(S2U0) <- names(result)
                        rownames(S2U0) <- "s2u"
                        
                        RESP <- data.frame(rbind(result, S2U0))
                } else if (length(lns2w) == 1){
                        lns2w <- est$par[c(k2 + 1):kk]
                        names(lns2w) <- c('lns2w')
                        cov.w <- ginv(est$hessian)[c(k2 + 1):kk,c(k2 + 1):kk]
                        S2W <- deltaMethod(lns2w, "exp(lns2w)", vcov = cov.w)
                        zS2W <- S2W$Estimate/S2W$SE
                        pS2W <- (1 - pnorm(abs(zS2W)))*2
                        icS2W <- lapply(length(S2W$Estimate), function(x) cbind(S2W$Estimate - qnorm(.975) * S2W$SE, S2W$Estimate + qnorm(.975) * S2W$SE))
                        S2W0 <- data.frame(S2W$Estimate, S2W$SE, zS2W, pS2W, icS2W[[1]][,1], icS2W[[1]][,2])
                        names(S2W0) <- names(result)
                        rownames(S2W0) <- "s2w"
                        
                        RESP <- data.frame(rbind(result, S2W0))
                } else {
                        RESP <- result
                }
                
                Yest <- X_fr %*% est$par[1:k1]
                s2u <- exp(sigma_u %*% est$par[c(k1 + 1):k2])
                s2w <- exp(sigma_w %*% est$par[c(k2 + 1):kk])
                erro <- Y_fr - Yest
                s2 <- s2u + s2w
                
                # Efficiency estimate by E(exp(-u)|e)
                mu.mod <- - erro - (s2w/sqrt(s2u)) # media
                s.mod <- sqrt(s2w) #  desvio padrao
                uf <- mu.mod + s.mod * ((dnorm(-mu.mod / s.mod)) / (pnorm(mu.mod / s.mod)))
                ef <- ((pnorm(-s.mod + mu.mod / s.mod)) / (pnorm(mu.mod / s.mod))) * exp(- mu.mod + 0.5 * s.mod^2)
                
                bias <- mean(Yest) - mean(Y_fr)
                RMSE <- sqrt(var(Yest) + bias^2)
                pearson <- cor(Yest, Y_fr, method = "pearson")
                spearman <- cor(Yest, Y_fr, method = "spearman")
                
                n <- length(Y_fr)
                K <- length(est$par)
                lnL <- -est$value
                AIC <- - 2 * lnL + 2 * K
                BIC <- - 2 * lnL + log(n) * K
        }
        lista <- list("efficiency" = ef, "error" = erro, "fitted.y" = Yest,
                      "table" = RESP, "summary.ef" = summary(ef), "sd.ef" = sd(ef),
                      "cor.pearson" = c(pearson), "cor.spearman" = c(spearman),
                      "value" = -est$value, "AIC" = AIC, "BIC" = BIC, 
                      "bias" = c(bias), "RMSE" = c(RMSE), 
                      "sample size" = n, "estimated parameters" = K)
        return(lista)
}

SF.trunc <- function(fr.form, mu.form, s2u.form = ~1, s2w.form = ~1, data = sys.parent()) {
        {
                # Frontier
                frontier <- fr.form
                mf <- model.frame(frontier, data)
                X_fr <- model.matrix(frontier, mf)
                Y_fr <- model.response(mf, "numeric")
                
                Mu <- mu.form
                Mi <- model.frame(Mu, data)
                mi <- model.matrix(Mu, Mi)
                
                # Sigma_u and Sigma_w
                Sigma_u <- s2u.form
                sig_u <- model.frame(Sigma_u, data)
                sigma_u <- model.matrix(Sigma_u, sig_u)
                
                Sigma_w <- s2w.form
                sig_w <- model.frame(Sigma_w, data)
                sigma_w <- model.matrix(Sigma_w, sig_w)
                
                # Initial values
                fr <- lm(frontier, data = data); names(fr$coefficients)[1] <- c("constant")
                mu.p <- c("constant" = 1, rep(0, ncol(mi) - 1)); if(length(names(mu.p)) > 1){names(mu.p)[2:length(mu.p)] <- colnames(mi)[2:length(mu.p)]}
                u.p <- c("constant" = 1)
                w.p <- c("constant" = 1)
                
                z <- c(coef(fr), mu.p, u.p, w.p)
                names(z) <- make.names(names(z), unique=T)
                
                kk <- length(z)
                k1 <- length(coef(fr))
                k2 <- k1 + length(mu.p)
                k3 <- k2 + length(u.p)
        }
        ll <- function(z) {
                x_p <- z[1:k1]
                mu_p <- z[c(k1 + 1):k2]
                u_p <- z[c(k2 + 1):k3]
                w_p <- z[c(k3 + 1):kk]
                
                Y <- X_fr %*% x_p
                mu <- mi %*% mu_p
                U <- sigma_u %*% u_p
                W <- sigma_w %*% w_p
                
                Sigma_u <- exp(U)
                Sigma_w <- exp(W)
                lambda <- sqrt(Sigma_u / Sigma_w)
                sigma <- Sigma_u + Sigma_w
                e <- Y_fr - Y
                
                aa <- mu / sqrt(sigma) * sqrt(1 + lambda^-2)
                pa <- pmax(pnorm(aa), 9.88131291682493e-324)
                bb <- (mu / lambda - e * lambda) / sqrt(sigma)
                pb <- pmax(pnorm(bb), 9.88131291682493e-324)
                cc <- (e + mu)^2 / sigma
                ly <- -0.5 * log(2 * pi) - 0.5 * log(sigma) - log(pa) + log(pb) - 0.5 * cc
                val <- sum(ly)
                
                return(-val)
        }
        G <- function(z) {
                x_p <- z[1:k1]
                mu_p <- z[c(k1 + 1):k2]
                u_p <- z[c(k2 + 1):k3]
                w_p <- z[c(k3 + 1):kk]
                
                Y <- X_fr %*% x_p
                mu <- mi %*% mu_p
                U <- sigma_u %*% u_p
                W <- sigma_w %*% w_p
                
                Sigma_u <- exp(U)
                Sigma_w <- exp(W)
                lambda <- sqrt(Sigma_u / Sigma_w)
                sigma <- Sigma_u + Sigma_w
                e <- Y_fr - Y
                
                aa <- mu / sqrt(sigma) * sqrt(1 + lambda^-2)
                da <- pmax(dnorm(aa), 9.88131291682493e-324)
                pa <- pmax(pnorm(aa), 9.88131291682493e-324)
                bb <- (mu / lambda - e * lambda) / sqrt(sigma)
                db <- pmax(dnorm(bb), 9.88131291682493e-324)
                pb <- pmax(pnorm(bb), 9.88131291682493e-324)
                cc <- (e + mu)^2 / sigma
                da_pa <- da/pa; db_pb <- db/pb
                valor <- (e + mu) / sigma + db_pb * lambda / sqrt(sigma)
                
                g.b <- t(X_fr) %*% valor
                
                g.mu <- t(mi) %*% as.vector(1 / sqrt(sigma) * (db_pb / lambda - da_pa * sqrt(1 + lambda^(-2)) - (e + mu) / sqrt(sigma)))
                
                Partu <- (2 * mu + e + mu / lambda^2) * db_pb / (lambda * sqrt(sigma))
                g.u <- t(sigma_u) %*% as.vector((0.5 * mu / Sigma_u^1.5 * da_pa + 0.5 / sigma * (cc - Partu - 1)) * Sigma_u)
                
                Partw <- lambda / sqrt(sigma) * db_pb * (mu + 2 * e + e * lambda^2)
                g.w <- t(sigma_w) %*% as.vector((0.5 / sigma * (cc + Partw - 1)) * Sigma_w)
                
                g <- c(g.b, g.mu, g.u, g.w)
                
                return(-g)
        }
        est <- optim(par = z, fn = ll, gr = G, hessian = TRUE, method = "BFGS",
                     control = list(fnscale = 1, trace = TRUE, REPORT = 1, maxit = 150000))
        #round(data.frame("numerical" = grad(ll,est$par), "analytical" = G(est$par)), 6)
        {
                ep <- sqrt(diag(ginv(est$hessian)))
                estZ <- lapply(length(est$par), function(x) est$par/ep)
                p_valor <- lapply(length(est$par), function(x) (1 - pnorm(abs(est$par/ep)))*2)
                ics <- lapply(length(est$par), function(x) cbind(est$par - qnorm(.975) * ep, est$par + qnorm(.975) * ep))
                result <- data.frame("Coefficient" = cbind(est$par), "Std.Err" = ep, "z" = estZ[[1]], 
                                     "P-value" = p_valor[[1]], "Lower limit" = ics[[1]][,1], "Upper limit" = ics[[1]][,2])
                
                lns2u <- est$par[c(k2 + 1):k3]
                names(lns2u) <- c('lns2u')
                cov.u <- ginv(est$hessian)[c(k2 + 1):k3,c(k2 + 1):k3]
                S2U <- deltaMethod(lns2u, "exp(lns2u)", vcov = cov.u)
                zS2U <- S2U$Estimate/S2U$SE
                pS2U <- (1 - pnorm(abs(zS2U)))*2
                icS2U <- lapply(length(S2U$Estimate), function(x) cbind(S2U$Estimate - qnorm(.975) * S2U$SE, S2U$Estimate + qnorm(.975) * S2U$SE))
                S2U0 <- data.frame(S2U$Estimate, S2U$SE, zS2U, pS2U, icS2U[[1]][,1], icS2U[[1]][,2])
                names(S2U0) <- names(result)
                rownames(S2U0) <- "s2u"
                
                lns2w <- est$par[c(k3 + 1):kk]
                names(lns2w) <- c('lns2w')
                cov.w <- ginv(est$hessian)[c(k3 + 1):kk,c(k3 + 1):kk]
                S2W <- deltaMethod(lns2w, "exp(lns2w)", vcov = cov.w)
                zS2W <- S2W$Estimate/S2W$SE
                pS2W <- (1 - pnorm(abs(zS2W)))*2
                icS2W <- lapply(length(S2W$Estimate), function(x) cbind(S2W$Estimate - qnorm(.975) * S2W$SE, S2W$Estimate + qnorm(.975) * S2W$SE))
                S2W0 <- data.frame(S2W$Estimate, S2W$SE, zS2W, pS2W, icS2W[[1]][,1], icS2W[[1]][,2])
                names(S2W0) <- names(result)
                rownames(S2W0) <- "s2w"
                
                s20 <- c(S2U$Estimate, S2W$Estimate)
                names(s20) <- c('s2u','s2w')
                cov.s <- ginv(est$hessian)[c(k2 + 1):kk,c(k2 + 1):kk]
                S2 <- deltaMethod(s20, "s2u + s2w", vcov = cov.s)
                zS2 <- S2$Estimate/S2$SE
                pS2 <- (1 - pnorm(abs(zS2)))*2
                icS2 <- lapply(length(S2$Estimate), function(x) cbind(S2$Estimate - qnorm(.975) * S2$SE, S2$Estimate + qnorm(.975) * S2$SE))
                S20 <- data.frame(S2$Estimate, S2$SE, zS2, pS2, icS2[[1]][,1], icS2[[1]][,2])
                names(S20) <- names(result)
                rownames(S20) <- "s2"
                
                Lambda <- deltaMethod(s20, "sqrt(s2u/s2w)", vcov = cov.s)
                zLam <- Lambda$Estimate/Lambda$SE
                pLam <- (1 - pnorm(abs(zLam)))*2
                icLam <- lapply(length(Lambda$Estimate), function(x) cbind(Lambda$Estimate - qnorm(.975) * Lambda$SE, Lambda$Estimate + qnorm(.975) * Lambda$SE))
                LAM <- data.frame(Lambda$Estimate, Lambda$SE, zLam, pLam, icLam[[1]][,1], icLam[[1]][,2])
                names(LAM) <- names(result)
                rownames(LAM) <- "lambda"
                
                RESP <- data.frame(rbind(result, S2U0, S2W0, S20, LAM))
                
                Yest <- X_fr %*% est$par[1:k1]
                mui <- mi %*% est$par[c(k1 + 1):k2]
                s2u <- exp(sigma_u %*% est$par[c(k2 + 1):k3])
                s2w <- exp(sigma_w %*% est$par[c(k3 + 1):kk])
                erro <- Y_fr - Yest
                s2 <- s2u + s2w
                
                # Efficiency estimate by E(exp(-u)|e)
                mu.mod <- (- erro * s2u + mui * s2w)/ s2 # media
                s.mod <- sqrt((s2u * s2w)/s2) # desvio padrao
                uf <- mu.mod + s.mod * ((dnorm(-mu.mod / s.mod)) / (pnorm(mu.mod / s.mod)))
                ef <- ((pnorm(-s.mod + mu.mod / s.mod)) / (pnorm(mu.mod / s.mod))) * exp(- mu.mod + 0.5 * s.mod^2)
                
                bias <- mean(Yest) - mean(Y_fr)
                RMSE <- sqrt(var(Yest) + bias^2)
                pearson <- cor(Yest, Y_fr, method = "pearson")
                spearman <- cor(Yest, Y_fr, method = "spearman")
                
                n <- length(Y_fr)
                K <- length(est$par)
                lnL <- -est$value
                AIC <- - 2 * lnL + 2 * K
                BIC <- - 2 * lnL + log(n) * K
        }
        lista <- list("efficiency" = ef, "error" = erro, "fitted.y" = Yest,
                      "table" = RESP, "summary.ef" = summary(ef), "sd.ef" = sd(ef),
                      "cor.pearson" = c(pearson), "cor.spearman" = c(spearman),
                      "value" = -est$value, "AIC" = AIC, "BIC" = BIC, 
                      "bias" = c(bias), "RMSE" = c(RMSE), 
                      "sample size" = n, "estimated parameters" = K)
        return(lista)
}


# Stochastic production frontier with endogenous variables in one-step

SF.half1S <- function(fr.form, end.form, s2u.form, s2w.form, data = sys.parent()) {
        {# Frontier
                frontier <- fr.form
                mf <- model.frame(frontier, data)
                X_fr <- model.matrix(frontier, mf)
                Y_fr <- model.response(mf, "numeric")
                
                # Endogenous variables
                End <- end.form
                m <- model.frame(End, data)
                xend <- model.matrix(End, m)
                yend <- model.response(m, "numeric")
                p <- ncol(as.matrix(yend))
                
                # Sigma_u and Sigma_w
                Sigma_u <- s2u.form
                sig_u <- model.frame(Sigma_u, data)
                sigma_u <- model.matrix(Sigma_u, sig_u)
                
                Sigma_w <- s2w.form
                sig_w <- model.frame(Sigma_w, data)
                sigma_w <- model.matrix(Sigma_w, sig_w)
                
                # Initial values
                ols <- lm(End, data); if(length(names(ols$coefficients)) > 1){names(ols$coefficients)[1] <- c("constant")} else{rownames(ols$coefficients)[1] <- c("constant")}
                fr <- lm(frontier, data = data); names(fr$coefficients)[1] <- c("constant")
                u.p <- c("constant" = 1, rep(0, ncol(sigma_u)-1)); if(length(names(u.p)) > 1){names(u.p)[2:length(u.p)] <- colnames(sigma_u)[2:length(u.p)]}
                w.p <- c("constant" = 1, rep(0, ncol(sigma_w)-1)); if(length(names(w.p)) > 1){names(w.p)[2:length(w.p)] <- colnames(sigma_w)[2:length(w.p)]}
                eta.p <- rep(0, p); names(eta.p) <- paste("eta", seq(1:p), sep = ".")
                
                ols.p <- c(ols$coefficients); names(ols.p) <- rep(rownames(ols$coefficients), p); if(length(names(ols$coefficients)) > 1){names(ols.p) <- names(ols$coefficients)}
                z <- c(ols.p, fr$coefficients, u.p, w.p, eta.p)
                names(z) <- make.names(names(z), unique=T)
                
                kk <- length(z)
                k1 <- length(coef(ols))
                k2 <- k1 + length(coef(fr))
                k3 <- k2 + length(u.p)
                k4 <- k3 + length(w.p)
        }
        ll <- function(z){
                z_p <- z[1:k1]
                y_p <- z[c(k1 + 1):k2]
                u_p <- z[c(k2 + 1):k3]
                w_p <- z[c(k3 + 1):k4]
                eta_p <- z[c(k4 + 1):kk]
                
                # Part x
                Z <- xend
                delta <- matrix(z_p, ncol = p)
                residual <- yend - Z %*% delta
                n <- nrow(residual)
                V <- crossprod(residual) / n
                
                ## Likelihood of x
                lx <- n * 0.5 * (- p * log(2 * pi) - log(det(V))) - 0.5 * sum(mahalanobis(residual, 0, V))
                
                ## Part y|x
                Y <- X_fr %*% y_p
                U <- sigma_u %*% u_p
                W <- sigma_w %*% w_p
                W0 <- sigma_w[,1] %*% as.matrix(w_p[1])
                Eta <- eta_p
                
                Sigma_u <- exp(U)
                Sigma_w <- exp(W)
                lambda <- sqrt(Sigma_u/Sigma_w)
                sigma <- Sigma_u + Sigma_w
                sigma_cw <- exp(W0)
                cw <- sqrt(Sigma_w/sigma_cw)
                omega <- residual %*% Eta
                
                e <- Y_fr - Y - cw * omega
                
                ## Likelihood of y|x, where x is the endogenous variable
                zz <- - e * lambda / sqrt(sigma)
                pz <- pmax(pnorm(zz), 9.88131291682493e-324)
                ly.x.half <- 0.5 * log(2 / pi) - 0.5 * log(sigma) + log(pz) - 0.5 * e^2 / sigma
                # Sum of the two likelihoods
                val <- sum(ly.x.half) + lx
                
                return(-val)
        }
        G <- function(z){
                z_p <- z[1:k1]
                y_p <- z[c(k1 + 1):k2]
                u_p <- z[c(k2 + 1):k3]
                w_p <- z[c(k3 + 1):k4]
                eta_p <- z[c(k4 + 1):kk]
                
                # Part x
                Z <- xend
                delta <- matrix(z_p, ncol = p)
                residual <- yend - Z %*% delta
                n <- nrow(residual)
                V <- crossprod(residual) / n
                
                ## Part y|x
                Y <- X_fr %*% y_p
                U <- sigma_u %*% u_p
                W <- sigma_w %*% w_p
                W0 <- sigma_w[,1] %*% as.matrix(w_p[1])
                Eta <- eta_p
                
                Sigma_u <- as.vector(exp(U))
                Sigma_w <- as.vector(exp(W))
                lambda <- sqrt(Sigma_u/Sigma_w)
                sigma <- Sigma_u + Sigma_w
                sigma_cw <- as.vector(exp(W0))
                cw <- sqrt(Sigma_w/sigma_cw)
                omega <- as.vector(residual %*% Eta)
                e <- as.vector(Y_fr - Y - cw * omega)
                
                zz <- - e * lambda / sqrt(sigma)
                dz <- pmax(dnorm(zz), 9.88131291682493e-324)
                pz <- pmax(pnorm(zz), 9.88131291682493e-324)
                fdp <- dz; cdf <- pz; fdp_cdf <- fdp / cdf
                valor <- e / sigma + lambda / sqrt(sigma) * fdp_cdf
                
                Part1 <- t(Z) %*% residual %*% ginv(V)
                Part2 <- - t(Eta %*% t(cw * valor) %*% Z)
                g.end <- Part1 + Part2
                
                g.b <- t(X_fr) %*% valor
                
                g.u <- t(sigma_u) %*% ((-0.5 / sigma * (1 - e^2 / sigma + e / sqrt(sigma) * 1/lambda * fdp/cdf)) * Sigma_u)
                
                sigma_w0 <- sigma_w
                sigma_w0[,1] <- sigma_w[,1] - c(1)
                s2w_x3 <- sigma_w * Sigma_w
                g.w <- colSums(((0.5 * s2w_x3 / sigma) * (e^2/sigma - 1 + 
                                                                  (fdp_cdf * e * lambda / sqrt(sigma)) * (2 + lambda^2))) +
                                       (sigma_w0 * 0.5 * omega * cw * valor))
                
                g.eta <- t(residual) %*% (cw * valor)
                
                g <- c(g.end, g.b, g.u, g.w, g.eta)
                
                return(-g)
        }
        est <- optim(par = z, fn = ll, gr = G, hessian = TRUE, method = "BFGS", 
                     control = list(fnscale = 1, trace = TRUE, REPORT = 1, maxit = 150000))
        #round(data.frame("numerical" = grad(ll,est$par), "analytical" = G(est$par)), 6)
        {
                ep <- sqrt(diag(ginv(est$hessian)))
                estZ <- lapply(length(est$par), function(x) est$par/ep)
                p_valor <- lapply(length(est$par), function(x) (1 - pnorm(abs(est$par/ep)))*2)
                ics <- lapply(length(est$par), function(x) cbind(est$par - qnorm(.975) * ep, est$par + qnorm(.975) * ep))
                result <- data.frame("Coefficient" = cbind(est$par), "Std.Err" = ep, "z" = estZ[[1]], 
                                     "P-value" = p_valor[[1]], "Lower limit" = ics[[1]][,1], "Upper limit" = ics[[1]][,2])
                lns2u <- est$par[c(k2 + 1):k3]
                lns2w <- est$par[c(k3 + 1):k4]
                if (length(lns2u) == 1 & length(lns2w) == 1){
                        names(lns2u) <- c('lns2u')
                        cov.u <- ginv(est$hessian)[c(k2 + 1):k3,c(k2 + 1):k3]
                        S2U <- deltaMethod(lns2u, "exp(lns2u)", vcov = cov.u)
                        zS2U <- S2U$Estimate/S2U$SE
                        pS2U <- (1 - pnorm(abs(zS2U)))*2
                        icS2U <- lapply(length(S2U$Estimate), function(x) cbind(S2U$Estimate - qnorm(.975) * S2U$SE, S2U$Estimate + qnorm(.975) * S2U$SE))
                        S2U0 <- data.frame(S2U$Estimate, S2U$SE, zS2U, pS2U, icS2U[[1]][,1], icS2U[[1]][,2])
                        names(S2U0) <- names(result)
                        rownames(S2U0) <- "s2u"
                        
                        names(lns2w) <- c('lns2w')
                        cov.w <- ginv(est$hessian)[c(k3 + 1):k4,c(k3 + 1):k4]
                        S2W <- deltaMethod(lns2w, "exp(lns2w)", vcov = cov.w)
                        zS2W <- S2W$Estimate/S2W$SE
                        pS2W <- (1 - pnorm(abs(zS2W)))*2
                        icS2W <- lapply(length(S2W$Estimate), function(x) cbind(S2W$Estimate - qnorm(.975) * S2W$SE, S2W$Estimate + qnorm(.975) * S2W$SE))
                        S2W0 <- data.frame(S2W$Estimate, S2W$SE, zS2W, pS2W, icS2W[[1]][,1], icS2W[[1]][,2])
                        names(S2W0) <- names(result)
                        rownames(S2W0) <- "s2w"
                        
                        s20 <- c(S2U$Estimate, S2W$Estimate)
                        names(s20) <- c('s2u','s2w')
                        cov.s <- ginv(est$hessian)[c(k2 + 1):k4,c(k2 + 1):k2]
                        S2 <- deltaMethod(s20, "s2u + s2w", vcov = cov.s)
                        zS2 <- S2$Estimate/S2$SE
                        pS2 <- (1 - pnorm(abs(zS2)))*2
                        icS2 <- lapply(length(S2$Estimate), function(x) cbind(S2$Estimate - qnorm(.975) * S2$SE, S2$Estimate + qnorm(.975) * S2$SE))
                        S20 <- data.frame(S2$Estimate, S2$SE, zS2, pS2, icS2[[1]][,1], icS2[[1]][,2])
                        names(S20) <- names(result)
                        rownames(S20) <- "s2"
                        
                        Lambda <- deltaMethod(s20, "sqrt(s2u/s2w)", vcov = cov.s)
                        zLam <- Lambda$Estimate/Lambda$SE
                        pLam <- (1 - pnorm(abs(zLam)))*2
                        icLam <- lapply(length(Lambda$Estimate), function(x) cbind(Lambda$Estimate - qnorm(.975) * Lambda$SE, Lambda$Estimate + qnorm(.975) * Lambda$SE))
                        LAM <- data.frame(Lambda$Estimate, Lambda$SE, zLam, pLam, icLam[[1]][,1], icLam[[1]][,2])
                        names(LAM) <- names(result)
                        rownames(LAM) <- "lambda"
                        
                        RESP <- data.frame(rbind(result, S2U0, S2W0, S20, LAM))
                } else if (length(lns2u) == 1){
                        lns2u <- est$par[c(k2 + 1):k3]
                        names(lns2u) <- c('lns2u')
                        cov.u <- ginv(est$hessian)[c(k2 + 1):k3,c(k2 + 1):k3]
                        S2U <- deltaMethod(lns2u, "exp(lns2u)", vcov = cov.u)
                        zS2U <- S2U$Estimate/S2U$SE
                        pS2U <- (1 - pnorm(abs(zS2U)))*2
                        icS2U <- lapply(length(S2U$Estimate), function(x) cbind(S2U$Estimate - qnorm(.975) * S2U$SE, S2U$Estimate + qnorm(.975) * S2U$SE))
                        S2U0 <- data.frame(S2U$Estimate, S2U$SE, zS2U, pS2U, icS2U[[1]][,1], icS2U[[1]][,2])
                        names(S2U0) <- names(result)
                        rownames(S2U0) <- "s2u"
                        
                        RESP <- data.frame(rbind(result, S2U0))
                } else if (length(lns2w) == 1){
                        lns2w <- est$par[c(k3 + 1):k4]
                        names(lns2w) <- c('lns2w')
                        cov.w <- ginv(est$hessian)[c(k3 + 1):k4,c(k3 + 1):k4]
                        S2W <- deltaMethod(lns2w, "exp(lns2w)", vcov = cov.w)
                        zS2W <- S2W$Estimate/S2W$SE
                        pS2W <- (1 - pnorm(abs(zS2W)))*2
                        icS2W <- lapply(length(S2W$Estimate), function(x) cbind(S2W$Estimate - qnorm(.975) * S2W$SE, S2W$Estimate + qnorm(.975) * S2W$SE))
                        S2W0 <- data.frame(S2W$Estimate, S2W$SE, zS2W, pS2W, icS2W[[1]][,1], icS2W[[1]][,2])
                        names(S2W0) <- names(result)
                        rownames(S2W0) <- "s2w"
                        
                        RESP <- data.frame(rbind(result, S2W0))
                } else {
                        RESP <- result
                }
                
                Z <- xend
                delta <- matrix(z[1:k1], ncol = p)
                residual <- yend - Z %*% delta
                Y_est <- X_fr %*% est$par[c(k1 + 1):k2]
                s2u <- exp(sigma_u %*% est$par[c(k2 + 1):k3])
                s2w <- exp(sigma_w %*% est$par[c(k3 + 1):k4])
                s2cw <- exp(est$par[c(k3 + 1):k4][1])
                etas <- est$par[c(k4 + 1):kk]
                yest <- Y_est + sqrt(s2w/s2cw) * residual %*% etas
                erro <- Y_fr - yest
                s2 <- s2u + s2w
                
                # Efficiency estimate by E(exp(-u)|e)
                mu.mod <- - erro * s2u / s2
                s.mod <- sqrt(s2w * s2u/s2)
                uf <- mu.mod + s.mod * ((dnorm(-mu.mod / s.mod)) / (pnorm(mu.mod / s.mod)))
                ef <- ((pnorm(-s.mod + mu.mod / s.mod)) / (pnorm(mu.mod / s.mod))) * exp(- mu.mod + 0.5 * s.mod^2)
                
                # LR test for endogeneity
                lnL1 <- -ll(est$par)
                lnL00 <- function(z) {
                        z_p <- z[1:k1]
                        y_p <- z[c(k1 + 1):k2]
                        u_p <- z[c(k2 + 1):k3]
                        w_p <- z[c(k3 + 1):k4]
                        eta_p <- rep(0,p)
                        
                        # Part x
                        Z <- xend
                        delta <- matrix(z_p, ncol = p)
                        residual <- yend - Z %*% delta
                        n <- nrow(residual)
                        V <- crossprod(residual) / n
                        
                        ## Likelihood of x
                        lx <- n * 0.5 * (- p * log(2 * pi) - log(det(V))) - 0.5 * sum(mahalanobis(residual, 0, V))
                        
                        ## Part y|x
                        Y <- X_fr %*% y_p
                        U <- sigma_u %*% u_p
                        W <- sigma_w %*% w_p
                        W0 <- sigma_w[,1] %*% as.matrix(w_p[1])
                        Eta <- eta_p
                        
                        Sigma_u <- exp(U)
                        Sigma_w <- exp(W)
                        lambda <- sqrt(Sigma_u/Sigma_w)
                        sigma <- Sigma_u + Sigma_w
                        sigma_cw <- exp(W0)
                        cw <- sqrt(Sigma_w/sigma_cw)
                        omega <- residual %*% Eta
                        
                        e <- Y_fr - Y - cw * omega
                        
                        ## Likelihood of y|x, where x is the endogenous variable
                        zz <- - e * lambda / sqrt(sigma)
                        pz <- pmax(pnorm(zz), 9.88131291682493e-324)
                        ly.x.half <- 0.5 * log(2 / pi) - 0.5 * log(sigma) + log(pz) - 0.5 * e^2 / sigma
                        # Sum of the two likelihoods
                        val <- sum(ly.x.half) + lx
                        
                        return(-val)
                }
                lnL0 <- -lnL00(est$par)
                LRtest <- 2 * (lnL1 - lnL0)
                LRp <- pchisq(LRtest, df = p, lower.tail = FALSE)
                
                # Wald test for endogeneity
                etas <- est$par[c(k4 + 1):kk]
                cov_etas <- ginv(est$hessian)[c(k4 + 1):kk, c(k4 + 1):kk]
                waldtest <- t(etas) %*% ginv(cov_etas) %*% etas
                waldp <- pchisq(waldtest, df = p, lower.tail = FALSE)
                
                delta_est <- matrix(est$par[1:k1], ncol = p)
                yend_est <- Z %*% delta_est
                if(ncol(yend_est) > 1){colnames(yend_est) <- paste(colnames(yend), "IV", sep = ".")}
                cor_p.iv <- cor(yend, yend_est, method = "pearson")
                cor_s.iv <- cor(yend, yend_est, method = "spearman")
                
                bias <- mean(Y_est) - mean(Y_fr)
                RMSE <- sqrt(var(Y_est) + bias^2)
                pearson <- cor(Y_est, Y_fr, method = "pearson")
                spearman <- cor(Y_est, Y_fr, method = "spearman")
                
                n <- length(Y_fr)
                K <- length(est$par)
                lnL <- -est$value
                AIC <- - 2 * lnL + 2 * K
                BIC <- - 2 * lnL + log(n) * K
        }
        lista <- list("efficiency" = ef, "error" = erro, "fitted.y_without.correction" = Y_est,
                      "fitted.y_with.correction" = yest, "table" = RESP, 
                      "summary.ef" = summary(ef), "sd.ef" = sd(ef, na.rm = TRUE),
                      "value" = -est$value, "AIC" = AIC, "BIC" = BIC, 
                      "cor pearson IV" = cor_p.iv, "cor spearman IV" = cor_s.iv,
                      "cor.pearson" = c(pearson), "cor.spearman" = c(spearman),
                      "LR chisq test" = matrix(cbind(LRtest, LRp), 1, 2, dimnames = list(c(),c("statistic", "P-value"))),
                      "Wald chisq test" = matrix(cbind(waldtest, waldp), 1, 2, dimnames = list(c(),c("statistic", "P-value"))),
                      "bias" = c(bias), "RMSE" = c(RMSE), 
                      "sample size" = n, "estimated parameters" = K)
        return(lista)
}

SF.exp1S <- function(fr.form, end.form, s2u.form, s2w.form, data = sys.parent()) {
        {# Frontier
                frontier <- fr.form
                mf <- model.frame(frontier, data)
                X_fr <- model.matrix(frontier, mf)
                Y_fr <- model.response(mf, "numeric")
                
                # Endogenous variables
                End <- end.form
                m <- model.frame(End, data)
                xend <- model.matrix(End, m)
                yend <- model.response(m, "numeric")
                p <- ncol(as.matrix(yend))
                
                # Sigma_u and Sigma_w
                Sigma_u <- s2u.form
                sig_u <- model.frame(Sigma_u, data)
                sigma_u <- model.matrix(Sigma_u, sig_u)
                
                Sigma_w <- s2w.form
                sig_w <- model.frame(Sigma_w, data)
                sigma_w <- model.matrix(Sigma_w, sig_w)
                
                # Initial values
                ols <- lm(End, data); if(length(names(ols$coefficients)) > 1){names(ols$coefficients)[1] <- c("constant")} else{rownames(ols$coefficients)[1] <- c("constant")}
                fr <- lm(frontier, data = data); names(fr$coefficients)[1] <- c("constant")
                u.p <- c("constant" = 1, rep(0, ncol(sigma_u)-1)); if(length(names(u.p)) > 1){names(u.p)[2:length(u.p)] <- colnames(sigma_u)[2:length(u.p)]}
                w.p <- c("constant" = 1, rep(0, ncol(sigma_w)-1)); if(length(names(w.p)) > 1){names(w.p)[2:length(w.p)] <- colnames(sigma_w)[2:length(w.p)]}
                eta.p <- rep(0, p); names(eta.p) <- paste("eta", seq(1:p), sep = ".")
                
                ols.p <- c(ols$coefficients); names(ols.p) <- rep(rownames(ols$coefficients), p); if(length(names(ols$coefficients)) > 1){names(ols.p) <- names(ols$coefficients)}
                z <- c(ols.p, fr$coefficients, u.p, w.p, eta.p)
                names(z) <- make.names(names(z), unique=T)
                
                kk <- length(z)
                k1 <- length(coef(ols))
                k2 <- k1 + length(coef(fr))
                k3 <- k2 + length(u.p)
                k4 <- k3 + length(w.p)
        }
        ll <- function(z){
                z_p <- z[1:k1]
                y_p <- z[c(k1 + 1):k2]
                u_p <- z[c(k2 + 1):k3]
                w_p <- z[c(k3 + 1):k4]
                eta_p <- z[c(k4 + 1):kk]
                
                # Part x
                Z <- xend
                delta <- matrix(z_p, ncol = p)
                residual <- yend - Z %*% delta
                n <- nrow(residual)
                V <- crossprod(residual) / n
                
                ## Likelihood of x
                lx <- n * 0.5 * (- p * log(2 * pi) - log(det(V))) - 0.5 * sum(mahalanobis(residual, 0, V))
                
                ## Part y|x
                Y <- X_fr %*% y_p
                U <- sigma_u %*% u_p
                W <- sigma_w %*% w_p
                W0 <- sigma_w[,1] %*% as.matrix(w_p[1])
                Eta <- eta_p
                
                Sigma_u <- exp(U)
                Sigma_w <- exp(W)
                lambda <- sqrt(Sigma_u/Sigma_w)
                sigma <- Sigma_u + Sigma_w
                sigma_cw <- exp(W0)
                cw <- sqrt(Sigma_w/sigma_cw)
                omega <- residual %*% Eta
                
                e <- Y_fr - Y - cw * omega
                
                ## Likelihood of y|x, where x is the endogenous variable
                zz <- (-e - Sigma_w / sqrt(Sigma_u)) / sqrt(Sigma_w)
                pz <- pnorm(zz, log = TRUE)
                ly.x.exp <- -0.5 * log(Sigma_u) + 0.5 * Sigma_w / Sigma_u + pz + e / sqrt(Sigma_u)
                
                # Sum of the two likelihoods
                val <- sum(ly.x.exp) + lx
                
                return(-val)
        }
        G <- function(z){
                z_p <- z[1:k1]
                y_p <- z[c(k1 + 1):k2]
                u_p <- z[c(k2 + 1):k3]
                w_p <- z[c(k3 + 1):k4]
                eta_p <- z[c(k4 + 1):kk]
                
                # Part x
                Z <- xend
                delta <- matrix(z_p, ncol = p)
                residual <- yend - Z %*% delta
                n <- nrow(residual)
                V <- crossprod(residual) / n
                
                ## Part y|x
                Y <- X_fr %*% y_p
                U <- sigma_u %*% u_p
                W <- sigma_w %*% w_p
                W0 <- sigma_w[,1] %*% as.matrix(w_p[1])
                Eta <- eta_p
                
                Sigma_u <- as.vector(exp(U))
                Sigma_w <- as.vector(exp(W))
                lambda <- sqrt(Sigma_u/Sigma_w)
                sigma <- Sigma_u + Sigma_w
                sigma_cw <- as.vector(exp(W0))
                cw <- sqrt(Sigma_w/sigma_cw)
                omega <- as.vector(residual %*% Eta)
                e <- as.vector(Y_fr - Y - cw * omega)
                
                zz <- (-e - Sigma_w / sqrt(Sigma_u)) / sqrt(Sigma_w)
                dz <- pmax(dnorm(zz), 9.88131291682493e-324)
                pz <- pmax(pnorm(zz), 9.88131291682493e-324)
                fdp <- dz; cdf <- pz; fdp_cdf <- fdp / cdf
                
                valor <- (fdp_cdf / sqrt(Sigma_w) - 1 / sqrt(Sigma_u))
                Part1 <- t(Z) %*% residual %*% ginv(V)
                Part2 <- t(- Eta %*% t(cw * valor) %*% Z)
                
                g.end <- Part1 + Part2
                
                g.b <- t(X_fr) %*% valor
                
                Partu <- fdp_cdf * sqrt(Sigma_w / Sigma_u) - Sigma_w / Sigma_u - e / sqrt(Sigma_u) - 1
                g.u <- t(sigma_u) %*% ((0.5 / Sigma_u * Partu) * Sigma_u)
                
                sigma_w0 <- sigma_w
                sigma_w0[,1] <- sigma_w[,1] - c(1)
                g.w <- t(sigma_w) %*% ((0.5 * (1 / Sigma_u + fdp_cdf / sqrt(Sigma_w) * (e / Sigma_w - 1 / sqrt(Sigma_u)))) * Sigma_w) +
                        t(sigma_w0) %*% (0.5 * omega * cw * valor)
                
                g.eta <- t(residual) %*% (cw * valor)
                
                g <- c(g.end, g.b, g.u, g.w, g.eta)
                
                return(-g)
        }
        est <- optim(par = z, fn = ll, gr = G, hessian = TRUE, method = "Nelder-Mead", 
                     control = list(fnscale = 1, trace = TRUE, REPORT = 1, maxit = 150000))
        #round(data.frame("numerical" = grad(ll,est$par), "analytical" = G(est$par)), 4)
        {
                ep <- sqrt(diag(ginv(est$hessian)))
                estZ <- lapply(length(est$par), function(x) est$par/ep)
                p_valor <- lapply(length(est$par), function(x) (1 - pnorm(abs(est$par/ep)))*2)
                ics <- lapply(length(est$par), function(x) cbind(est$par - qnorm(.975) * ep, est$par + qnorm(.975) * ep))
                result <- data.frame("Coefficient" = cbind(est$par), "Std.Err" = ep, "z" = estZ[[1]], 
                                     "P-value" = p_valor[[1]], "Lower limit" = ics[[1]][,1], "Upper limit" = ics[[1]][,2])
                lns2u <- est$par[c(k2 + 1):k3]
                lns2w <- est$par[c(k3 + 1):k4]
                if (length(lns2u) == 1 & length(lns2w) == 1){
                        names(lns2u) <- c('lns2u')
                        cov.u <- ginv(est$hessian)[c(k2 + 1):k3,c(k2 + 1):k3]
                        S2U <- deltaMethod(lns2u, "exp(lns2u)", vcov = cov.u)
                        zS2U <- S2U$Estimate/S2U$SE
                        pS2U <- (1 - pnorm(abs(zS2U)))*2
                        icS2U <- lapply(length(S2U$Estimate), function(x) cbind(S2U$Estimate - qnorm(.975) * S2U$SE, S2U$Estimate + qnorm(.975) * S2U$SE))
                        S2U0 <- data.frame(S2U$Estimate, S2U$SE, zS2U, pS2U, icS2U[[1]][,1], icS2U[[1]][,2])
                        names(S2U0) <- names(result)
                        rownames(S2U0) <- "s2u"
                        
                        names(lns2w) <- c('lns2w')
                        cov.w <- ginv(est$hessian)[c(k3 + 1):k4,c(k3 + 1):k4]
                        S2W <- deltaMethod(lns2w, "exp(lns2w)", vcov = cov.w)
                        zS2W <- S2W$Estimate/S2W$SE
                        pS2W <- (1 - pnorm(abs(zS2W)))*2
                        icS2W <- lapply(length(S2W$Estimate), function(x) cbind(S2W$Estimate - qnorm(.975) * S2W$SE, S2W$Estimate + qnorm(.975) * S2W$SE))
                        S2W0 <- data.frame(S2W$Estimate, S2W$SE, zS2W, pS2W, icS2W[[1]][,1], icS2W[[1]][,2])
                        names(S2W0) <- names(result)
                        rownames(S2W0) <- "s2w"
                        
                        s20 <- c(S2U$Estimate, S2W$Estimate)
                        names(s20) <- c('s2u','s2w')
                        cov.s <- ginv(est$hessian)[c(k2 + 1):k4,c(k2 + 1):k2]
                        S2 <- deltaMethod(s20, "s2u + s2w", vcov = cov.s)
                        zS2 <- S2$Estimate/S2$SE
                        pS2 <- (1 - pnorm(abs(zS2)))*2
                        icS2 <- lapply(length(S2$Estimate), function(x) cbind(S2$Estimate - qnorm(.975) * S2$SE, S2$Estimate + qnorm(.975) * S2$SE))
                        S20 <- data.frame(S2$Estimate, S2$SE, zS2, pS2, icS2[[1]][,1], icS2[[1]][,2])
                        names(S20) <- names(result)
                        rownames(S20) <- "s2"
                        
                        Lambda <- deltaMethod(s20, "sqrt(s2u/s2w)", vcov = cov.s)
                        zLam <- Lambda$Estimate/Lambda$SE
                        pLam <- (1 - pnorm(abs(zLam)))*2
                        icLam <- lapply(length(Lambda$Estimate), function(x) cbind(Lambda$Estimate - qnorm(.975) * Lambda$SE, Lambda$Estimate + qnorm(.975) * Lambda$SE))
                        LAM <- data.frame(Lambda$Estimate, Lambda$SE, zLam, pLam, icLam[[1]][,1], icLam[[1]][,2])
                        names(LAM) <- names(result)
                        rownames(LAM) <- "lambda"
                        
                        RESP <- data.frame(rbind(result, S2U0, S2W0, S20, LAM))
                } else if (length(lns2u) == 1){
                        lns2u <- est$par[c(k2 + 1):k3]
                        names(lns2u) <- c('lns2u')
                        cov.u <- ginv(est$hessian)[c(k2 + 1):k3,c(k2 + 1):k3]
                        S2U <- deltaMethod(lns2u, "exp(lns2u)", vcov = cov.u)
                        zS2U <- S2U$Estimate/S2U$SE
                        pS2U <- (1 - pnorm(abs(zS2U)))*2
                        icS2U <- lapply(length(S2U$Estimate), function(x) cbind(S2U$Estimate - qnorm(.975) * S2U$SE, S2U$Estimate + qnorm(.975) * S2U$SE))
                        S2U0 <- data.frame(S2U$Estimate, S2U$SE, zS2U, pS2U, icS2U[[1]][,1], icS2U[[1]][,2])
                        names(S2U0) <- names(result)
                        rownames(S2U0) <- "s2u"
                        
                        RESP <- data.frame(rbind(result, S2U0))
                } else if (length(lns2w) == 1){
                        lns2w <- est$par[c(k3 + 1):k4]
                        names(lns2w) <- c('lns2w')
                        cov.w <- ginv(est$hessian)[c(k3 + 1):k4,c(k3 + 1):k4]
                        S2W <- deltaMethod(lns2w, "exp(lns2w)", vcov = cov.w)
                        zS2W <- S2W$Estimate/S2W$SE
                        pS2W <- (1 - pnorm(abs(zS2W)))*2
                        icS2W <- lapply(length(S2W$Estimate), function(x) cbind(S2W$Estimate - qnorm(.975) * S2W$SE, S2W$Estimate + qnorm(.975) * S2W$SE))
                        S2W0 <- data.frame(S2W$Estimate, S2W$SE, zS2W, pS2W, icS2W[[1]][,1], icS2W[[1]][,2])
                        names(S2W0) <- names(result)
                        rownames(S2W0) <- "s2w"
                        
                        RESP <- data.frame(rbind(result, S2W0))
                } else {
                        RESP <- result
                }
                
                Z <- xend
                delta <- matrix(z[1:k1], ncol = p)
                residual <- yend - Z %*% delta
                Y_est <- X_fr %*% est$par[c(k1 + 1):k2]
                s2u <- exp(sigma_u %*% est$par[c(k2 + 1):k3])
                s2w <- exp(sigma_w %*% est$par[c(k3 + 1):k4])
                s2cw <- exp(est$par[c(k3 + 1):k4][1])
                etas <- est$par[c(k4 + 1):kk]
                yest <- Y_est + sqrt(s2w/s2cw) * residual %*% etas
                erro <- Y_fr - yest
                s2 <- s2u + s2w
                
                # Efficiency estimate by E(exp(-u)|e)
                mu.mod <- - erro - (s2w/sqrt(s2u))
                s.mod <- sqrt(s2w)
                uf <- mu.mod + s.mod * ((dnorm(-mu.mod / s.mod)) / (pnorm(mu.mod / s.mod)))
                ef <- ((pnorm(-s.mod + mu.mod / s.mod)) / (pnorm(mu.mod / s.mod))) * exp(- mu.mod + 0.5 * s.mod^2)
                
                # LR test for endogeneity
                lnL1 <- -ll(est$par)
                lnL00 <- function(z) {
                        z_p <- z[1:k1]
                        y_p <- z[c(k1 + 1):k2]
                        u_p <- z[c(k2 + 1):k3]
                        w_p <- z[c(k3 + 1):k4]
                        eta_p <- rep(0,p)
                        
                        # Part x
                        Z <- xend
                        delta <- matrix(z_p, ncol = p)
                        residual <- yend - Z %*% delta
                        n <- nrow(residual)
                        V <- crossprod(residual) / n
                        
                        ## Likelihood of x
                        lx <- n * 0.5 * (- p * log(2 * pi) - log(det(V))) - 0.5 * sum(mahalanobis(residual, 0, V))
                        
                        ## Part y|x
                        Y <- X_fr %*% y_p
                        U <- sigma_u %*% u_p
                        W <- sigma_w %*% w_p
                        W0 <- sigma_w[,1] %*% as.matrix(w_p[1])
                        Eta <- eta_p
                        
                        Sigma_u <- exp(U)
                        Sigma_w <- exp(W)
                        lambda <- sqrt(Sigma_u/Sigma_w)
                        sigma <- Sigma_u + Sigma_w
                        sigma_cw <- exp(W0)
                        cw <- sqrt(Sigma_w/sigma_cw)
                        omega <- residual %*% Eta
                        
                        e <- Y_fr - Y - cw * omega
                        
                        ## Likelihood of y|x, where x is the endogenous variable
                        zz <- (-e - Sigma_w / sqrt(Sigma_u)) / sqrt(Sigma_w)
                        pz <- pnorm(zz, log = TRUE)
                        ly.x.exp <- -0.5 * log(Sigma_u) + 0.5 * Sigma_w / Sigma_u + pz + e / sqrt(Sigma_u)
                        
                        # Sum of the two likelihoods
                        val <- sum(ly.x.exp) + lx
                        
                        return(-val)
                }
                lnL0 <- -lnL00(est$par)
                LRtest <- 2 * (lnL1 - lnL0)
                LRp <- pchisq(LRtest, df = p, lower.tail = FALSE)
                
                # Wald test for endogeneity
                etas <- est$par[c(k4 + 1):kk]
                cov_etas <- ginv(est$hessian)[c(k4 + 1):kk, c(k4 + 1):kk]
                waldtest <- t(etas) %*% ginv(cov_etas) %*% etas
                waldp <- pchisq(waldtest, df = p, lower.tail = FALSE)
                
                delta_est <- matrix(est$par[1:k1], ncol = p)
                yend_est <- Z %*% delta_est
                if(ncol(yend_est) > 1){colnames(yend_est) <- paste(colnames(yend), "IV", sep = ".")}
                cor_p.iv <- cor(yend, yend_est, method = "pearson")
                cor_s.iv <- cor(yend, yend_est, method = "spearman")
                
                bias <- mean(Y_est) - mean(Y_fr)
                RMSE <- sqrt(var(Y_est) + bias^2)
                pearson <- cor(Y_est, Y_fr, method = "pearson")
                spearman <- cor(Y_est, Y_fr, method = "spearman")
                
                n <- length(Y_fr)
                K <- length(est$par)
                lnL <- -est$value
                AIC <- - 2 * lnL + 2 * K
                BIC <- - 2 * lnL + log(n) * K
        }
        lista <- list("efficiency" = ef, "error" = erro, "fitted.y_without.correction" = Y_est,
                      "fitted.y_with.correction" = yest, "table" = RESP,
                      "summary.ef" = summary(ef), "sd.ef" = sd(ef, na.rm = TRUE),
                      "value" = -est$value, "AIC" = AIC, "BIC" = BIC, 
                      "cor pearson IV" = cor_p.iv, "cor spearman IV" = cor_s.iv,
                      "cor.pearson" = c(pearson), "cor.spearman" = c(spearman),
                      "LR chisq test" = matrix(cbind(LRtest, LRp), 1, 2, dimnames = list(c(),c("statistic", "P-value"))),
                      "Wald chisq test" = matrix(cbind(waldtest, waldp), 1, 2, dimnames = list(c(),c("statistic", "P-value"))),
                      "bias" = c(bias), "RMSE" = c(RMSE), 
                      "sample size" = n, "estimated parameters" = K)
        return(lista)
}

SF.trunc1S <- function(fr.form, end.form, mu.form, s2u.form = ~1, s2w.form = ~1, data = sys.parent()) {
        {# Frontier
                frontier <- fr.form
                mf <- model.frame(frontier, data)
                X_fr <- model.matrix(frontier, mf)
                Y_fr <- model.response(mf, "numeric")
                
                # Endogenous variables
                End <- end.form
                m <- model.frame(End, data)
                xend <- model.matrix(End, m)
                yend <- model.response(m, "numeric")
                p <- ncol(as.matrix(yend))
                
                Mu <- mu.form
                Mi <- model.frame(Mu, data)
                mi <- model.matrix(Mu, Mi)
                
                # Sigma_u and Sigma_w
                Sigma_u <- s2u.form
                sig_u <- model.frame(Sigma_u, data)
                sigma_u <- model.matrix(Sigma_u, sig_u)
                
                Sigma_w <- s2w.form
                sig_w <- model.frame(Sigma_w, data)
                sigma_w <- model.matrix(Sigma_w, sig_w)
                
                # Initial values
                ols <- lm(End, data); if(length(names(ols$coefficients)) > 1){names(ols$coefficients)[1] <- c("constant")} else{rownames(ols$coefficients)[1] <- c("constant")}
                fr <- lm(frontier, data = data); names(fr$coefficients)[1] <- c("constant")
                mu.p <- c("constant" = 1, rep(0, ncol(mi) - 1)); if(length(names(mu.p)) > 1){names(mu.p)[2:length(mu.p)] <- colnames(mi)[2:length(mu.p)]}
                u.p <- c('constant' = 1)
                w.p <- c('constant' = 1)
                eta.p <- rep(0, p); names(eta.p) <- paste("eta", seq(1:p), sep = ".")
                
                ols.p <- c(ols$coefficients); names(ols.p) <- rep(rownames(ols$coefficients), p); if(length(names(ols$coefficients)) > 1){names(ols.p) <- names(ols$coefficients)}
                z <- c(ols.p, fr$coefficients, mu.p, u.p, w.p, eta.p)
                names(z) <- make.names(names(z), unique=T)
                
                kk <- length(z)
                k1 <- length(coef(ols))
                k2 <- k1 + length(coef(fr))
                k3 <- k2 + length(mu.p)
                k4 <- k3 + length(u.p)
                k5 <- k4 + length(w.p)
        }
        ll <- function(z){
                z_p <- z[1:k1]
                y_p <- z[c(k1 + 1):k2]
                mu_p <- z[c(k2 + 1):k3]
                u_p <- z[c(k3 + 1):k4]
                w_p <- z[c(k4 + 1):k5]
                eta_p <- z[c(k5 + 1):kk]
                
                # Part x
                Z <- xend
                delta <- matrix(z_p, ncol = p)
                residual <- yend - Z %*% delta
                n <- nrow(residual)
                V <- crossprod(residual) / n
                
                ## Likelihood of x
                lx <- n * 0.5 * (- p * log(2 * pi) - log(det(V))) - 0.5 * sum(mahalanobis(residual, 0, V))
                
                ## Part y|x
                Y <- X_fr %*% y_p
                mu <- mi %*% mu_p
                U <- sigma_u %*% u_p
                W <- sigma_w %*% w_p
                W0 <- sigma_w[,1] %*% as.matrix(w_p[1])
                Eta <- eta_p
                
                Sigma_u <- exp(U)
                Sigma_w <- exp(W)
                lambda <- sqrt(Sigma_u/Sigma_w)
                sigma <- Sigma_u + Sigma_w
                sigma_cw <- exp(W0)
                cw <- sqrt(Sigma_w/sigma_cw)
                omega <- residual %*% Eta
                
                e <- Y_fr - Y - cw * omega
                
                aa <- mu / sqrt(sigma) * sqrt(1 + lambda^-2)
                pa <- pmax(pnorm(aa), 9.88131291682493e-324)
                bb <- (mu / lambda - e * lambda) / sqrt(sigma)
                pb <- pmax(pnorm(bb), 9.88131291682493e-324)
                cc <- (e + mu)^2 / sigma
                
                ## Likelihood of y|x, where x is the endogenous variable
                ly.x.trun <- -0.5 * log(2 * pi) - 0.5 * log(sigma) - log(pa) + log(pb) - 0.5 * cc
                
                # Sum of the two likelihoods
                val <- sum(ly.x.trun) + lx
                
                return(-val)
        }
        G <- function(z){
                z_p <- z[1:k1]
                y_p <- z[c(k1 + 1):k2]
                mu_p <- z[c(k2 + 1):k3]
                u_p <- z[c(k3 + 1):k4]
                w_p <- z[c(k4 + 1):k5]
                eta_p <- z[c(k5 + 1):kk]
                
                # Part x
                Z <- xend
                delta <- matrix(z_p, ncol = p)
                residual <- yend - Z %*% delta
                n <- nrow(residual)
                V <- crossprod(residual) / n
                
                ## Part y|x
                Y <- X_fr %*% y_p
                mu <- mi %*% mu_p
                U <- sigma_u %*% u_p
                W <- sigma_w %*% w_p
                W0 <- sigma_w[,1] %*% as.matrix(w_p[1])
                Eta <- eta_p
                
                Sigma_u <- as.vector(exp(U))
                Sigma_w <- as.vector(exp(W))
                lambda <- sqrt(Sigma_u/Sigma_w)
                sigma <- Sigma_u + Sigma_w
                sigma_cw <- as.vector(exp(W0))
                cw <- sqrt(Sigma_w/sigma_cw)
                omega <- as.vector(residual %*% Eta)
                e <- as.vector(Y_fr - Y - cw * omega)
                
                aa <- mu / sqrt(sigma) * sqrt(1 + lambda^-2)
                da <- pmax(dnorm(aa), 9.88131291682493e-324)
                pa <- pmax(pnorm(aa), 9.88131291682493e-324)
                
                bb <- as.vector((mu / lambda - e * lambda) / sqrt(sigma))
                db <- pmax(dnorm(bb), 9.88131291682493e-324)
                pb <- pmax(pnorm(bb), 9.88131291682493e-324)
                
                cc <- as.vector((e + mu)^2 / sigma)
                da_pa <- da/pa; db_pb <- db/pb
                
                valor <- as.vector((e + mu) / sigma + (lambda / sqrt(sigma)) * db_pb)
                
                Part1 <- t(Z) %*% residual %*% ginv(V)
                Part2 <- - t(Eta %*% t(cw * valor) %*% Z)
                g.end <- Part1 + Part2
                
                g.b <- t(X_fr) %*% valor
                
                g.mu <- t(mi) %*% (db_pb / (lambda * sqrt(sigma)) - sqrt((1 + lambda^(-2)) / sigma) * da_pa - (e + mu) / sigma)
                
                Partu <- db_pb / (lambda * sqrt(sigma)) * (2 * mu + e + mu / lambda^2)
                g.u <- t(sigma_u) %*% ((0.5 * mu / Sigma_u^1.5 * da_pa + 0.5 / sigma * (cc - Partu - 1)) * Sigma_u)
                
                Partw <- as.vector(lambda / sqrt(sigma) * db_pb * (mu + 2 * e + e * lambda^2))
                g.w <- t(Sigma_w) %*% (sigma_w * (0.5 / sigma * (cc + Partw - 1)))
                
                g.eta <- t(residual) %*% (cw * valor)
                
                g <- c(g.end, g.b, g.mu, g.u, g.w, g.eta)
                
                return(-g)
        }
        est <- optim(par = z, fn = ll, gr = G, hessian = TRUE, method = "BFGS", 
                     control = list(fnscale = 1, trace = TRUE, REPORT = 1, maxit = 150000))
        #round(data.frame("numerical" = grad(ll,est$par), "analytical" = G(est$par)), 4)
        {
                ep <- sqrt(diag(ginv(est$hessian)))
                estZ <- lapply(length(est$par), function(x) est$par/ep)
                p_valor <- lapply(length(est$par), function(x) (1 - pnorm(abs(est$par/ep)))*2)
                ics <- lapply(length(est$par), function(x) cbind(est$par - qnorm(.975) * ep, est$par + qnorm(.975) * ep))
                result <- data.frame("Coefficient" = cbind(est$par), "Std.Err" = ep, "z" = estZ[[1]], 
                                     "P-value" = p_valor[[1]], "Lower limit" = ics[[1]][,1], "Upper limit" = ics[[1]][,2])
                
                lns2u <- est$par[c(k3 + 1):k4]
                names(lns2u) <- c('lns2u')
                cov.u <- ginv(est$hessian)[c(k3 + 1):k4,c(k3 + 1):k4]
                S2U <- deltaMethod(lns2u, "exp(lns2u)", vcov = cov.u)
                zS2U <- S2U$Estimate/S2U$SE
                pS2U <- (1 - pnorm(abs(zS2U)))*2
                icS2U <- lapply(length(S2U$Estimate), function(x) cbind(S2U$Estimate - qnorm(.975) * S2U$SE, S2U$Estimate + qnorm(.975) * S2U$SE))
                S2U0 <- data.frame(S2U$Estimate, S2U$SE, zS2U, pS2U, icS2U[[1]][,1], icS2U[[1]][,2])
                names(S2U0) <- names(result)
                rownames(S2U0) <- "s2u"
                
                lns2w <- est$par[c(k4 + 1):k5]
                names(lns2w) <- c('lns2w')
                cov.w <- ginv(est$hessian)[c(k4 + 1):k5,c(k4 + 1):k5]
                S2W <- deltaMethod(lns2w, "exp(lns2w)", vcov = cov.w)
                zS2W <- S2W$Estimate/S2W$SE
                pS2W <- (1 - pnorm(abs(zS2W)))*2
                icS2W <- lapply(length(S2W$Estimate), function(x) cbind(S2W$Estimate - qnorm(.975) * S2W$SE, S2W$Estimate + qnorm(.975) * S2W$SE))
                S2W0 <- data.frame(S2W$Estimate, S2W$SE, zS2W, pS2W, icS2W[[1]][,1], icS2W[[1]][,2])
                names(S2W0) <- names(result)
                rownames(S2W0) <- "s2w"
                
                s20 <- c(S2U$Estimate, S2W$Estimate)
                names(s20) <- c('s2u','s2w')
                cov.s <- ginv(est$hessian)[c(k3 + 1):k5,c(k3 + 1):k5]
                S2 <- deltaMethod(s20, "s2u + s2w", vcov = cov.s)
                zS2 <- S2$Estimate/S2$SE
                pS2 <- (1 - pnorm(abs(zS2)))*2
                icS2 <- lapply(length(S2$Estimate), function(x) cbind(S2$Estimate - qnorm(.975) * S2$SE, S2$Estimate + qnorm(.975) * S2$SE))
                S20 <- data.frame(S2$Estimate, S2$SE, zS2, pS2, icS2[[1]][,1], icS2[[1]][,2])
                names(S20) <- names(result)
                rownames(S20) <- "s2"
                
                Lambda <- deltaMethod(s20, "sqrt(s2u/s2w)", vcov = cov.s)
                zLam <- Lambda$Estimate/Lambda$SE
                pLam <- (1 - pnorm(abs(zLam)))*2
                icLam <- lapply(length(Lambda$Estimate), function(x) cbind(Lambda$Estimate - qnorm(.975) * Lambda$SE, Lambda$Estimate + qnorm(.975) * Lambda$SE))
                LAM <- data.frame(Lambda$Estimate, Lambda$SE, zLam, pLam, icLam[[1]][,1], icLam[[1]][,2])
                names(LAM) <- names(result)
                rownames(LAM) <- "lambda"
                
                RESP <- data.frame(rbind(result, S2U0, S2W0, S20, LAM))
                
                Z <- xend
                delta <- matrix(z[1:k1], ncol = p)
                residual <- yend - Z %*% delta
                Y_est <- X_fr %*% est$par[c(k1 + 1):k2]
                mui <- mi %*% est$par[c(k2 + 1):k3]
                s2u <- exp(sigma_u %*% est$par[c(k3 + 1):k4])
                s2w <- exp(sigma_w %*% est$par[c(k4 + 1):k5])
                s2cw <- exp(est$par[c(k4 + 1):k5][1])
                etas <- est$par[c(k5 + 1):kk]
                yest <- Y_est + sqrt(s2w/s2cw) * residual %*% etas
                erro <- Y_fr - yest
                s2 <- s2u + s2w
                
                # Efficiency estimate by E(exp(-u)|e)
                mu.mod <- (- erro * s2u + mui * s2w)/ s2
                s.mod <- sqrt((s2u * s2w)/s2)
                uf <- mu.mod + s.mod * ((dnorm(-mu.mod / s.mod)) / (pnorm(mu.mod / s.mod)))
                ef <- ((pnorm(-s.mod + mu.mod / s.mod)) / (pnorm(mu.mod / s.mod))) * exp(- mu.mod + 0.5 * s.mod^2)
                
                # LR test for endogeneity
                lnL1 <- -ll(est$par)
                lnL00 <- function(z) {
                        z_p <- z[1:k1]
                        y_p <- z[c(k1 + 1):k2]
                        mu_p <- z[c(k2 + 1):k3]
                        u_p <- z[c(k3 + 1):k4]
                        w_p <- z[c(k4 + 1):k5]
                        eta_p <- rep(0,p)
                        
                        # Part x
                        Z <- xend
                        delta <- matrix(z_p, ncol = p)
                        residual <- yend - Z %*% delta
                        n <- nrow(residual)
                        V <- crossprod(residual) / n
                        
                        ## Likelihood of x
                        lx <- n * 0.5 * (- p * log(2 * pi) - log(det(V))) - 0.5 * sum(mahalanobis(residual, 0, V))
                        
                        ## Part y|x
                        Y <- X_fr %*% y_p
                        mu <- mi %*% mu_p
                        U <- sigma_u %*% u_p
                        W <- sigma_w %*% w_p
                        W0 <- sigma_w[,1] %*% as.matrix(w_p[1])
                        Eta <- eta_p
                        
                        Sigma_u <- exp(U)
                        Sigma_w <- exp(W)
                        lambda <- sqrt(Sigma_u/Sigma_w)
                        sigma <- Sigma_u + Sigma_w
                        sigma_cw <- exp(W0)
                        cw <- sqrt(Sigma_w/sigma_cw)
                        omega <- residual %*% Eta
                        
                        e <- Y_fr - Y - cw * omega
                        
                        aa <- mu / sqrt(sigma) * sqrt(1 + lambda^-2)
                        pa <- pmax(pnorm(aa), 9.88131291682493e-324)
                        bb <- (mu / lambda - e * lambda) / sqrt(sigma)
                        pb <- pmax(pnorm(bb), 9.88131291682493e-324)
                        cc <- (e + mu)^2 / sigma
                        
                        ## Likelihood of y|x, where x is the endogenous variable
                        ly.x.trun <- -0.5 * log(2 * pi) - 0.5 * log(sigma) - log(pa) + log(pb) - 0.5 * cc
                        
                        # Sum of the two likelihoods
                        val <- sum(ly.x.trun) + lx
                        
                        return(-val)
                }
                lnL0 <- -lnL00(est$par)
                LRtest <- 2 * (lnL1 - lnL0)
                LRp <- pchisq(LRtest, df = p, lower.tail = FALSE)
                
                # Wald test for endogeneity
                etas <- est$par[c(k5 + 1):kk]
                cov_etas <- ginv(est$hessian)[c(k5 + 1):kk, c(k5 + 1):kk]
                waldtest <- t(etas) %*% ginv(cov_etas) %*% etas
                waldp <- pchisq(waldtest, df = p, lower.tail = FALSE)
                
                delta_est <- matrix(est$par[1:k1], ncol = p)
                yend_est <- Z %*% delta_est
                if(ncol(yend_est) > 1){colnames(yend_est) <- paste(colnames(yend), "IV", sep = ".")}
                cor_p.iv <- cor(yend, yend_est, method = "pearson")
                cor_s.iv <- cor(yend, yend_est, method = "spearman")
                
                bias <- mean(Y_est) - mean(Y_fr)
                RMSE <- sqrt(var(Y_est) + bias^2)
                pearson <- cor(Y_est, Y_fr, method = "pearson")
                spearman <- cor(Y_est, Y_fr, method = "spearman")
                
                n <- length(Y_fr)
                K <- length(est$par)
                lnL <- -est$value
                AIC <- - 2 * lnL + 2 * K
                BIC <- - 2 * lnL + log(n) * K
        }
        lista <- list("efficiency" = ef, "error" = erro, "fitted.y_without.correction" = Y_est,
                      "fitted.y_with.correction" = yest, "table" = RESP,
                      "summary.ef" = summary(ef), "sd.ef" = sd(ef, na.rm = TRUE),
                      "value" = -est$value, "AIC" = AIC, "BIC" = BIC, 
                      "cor pearson IV" = cor_p.iv, "cor spearman IV" = cor_s.iv,
                      "cor.pearson" = c(pearson), "cor.spearman" = c(spearman),
                      "LR chisq test" = matrix(cbind(LRtest, LRp), 1, 2, dimnames = list(c(),c("statistic", "P-value"))),
                      "Wald chisq test" = matrix(cbind(waldtest, waldp), 1, 2, dimnames = list(c(),c("statistic", "P-value"))),
                      "bias" = c(bias), "RMSE" = c(RMSE), 
                      "sample size" = n, "estimated parameters" = K)
        return(lista)
}


# Stochastic production frontier with endogenous variables in two-step

SF.half2S <- function(fr.form, end.form, s2u.form, s2w.form, data = sys.parent()) {
        {# Frontier
                frontier <- fr.form
                mf <- model.frame(frontier, data)
                X_fr <- model.matrix(frontier, mf)
                Y_fr <- model.response(mf, "numeric")
                
                # Endogenous variables
                End <- end.form
                m <- model.frame(End, data)
                xend <- model.matrix(End, m)
                yend <- model.response(m, "numeric")
                p <- ncol(as.matrix(yend))
                
                # Sigma_u and Sigma_w
                Sigma_u <- s2u.form
                sig_u <- model.frame(Sigma_u, data)
                sigma_u <- model.matrix(Sigma_u, sig_u)
                
                Sigma_w <- s2w.form
                sig_w <- model.frame(Sigma_w, data)
                sigma_w <- model.matrix(Sigma_w, sig_w)
                
                # Step 1 - OLS of the Endogenous Variables
                ols <- lm(End, data); if(length(names(ols$coefficients)) > 1){names(ols$coefficients)[1] <- c("constant")} else{rownames(ols$coefficients)[1] <- c("constant")}
                residual <- as.matrix(residuals(ols))
                
                yend_est <- as.matrix(fitted(ols))
                if(ncol(yend_est) > 1){colnames(yend_est) <- paste(colnames(yend), "IV", sep = ".")}
                cor_p.iv <- cor(yend, yend_est, method = "pearson")
                cor_s.iv <- cor(yend, yend_est, method = "spearman")
                
                result_IV <- summary(ols)
                
                # Initial values
                fr <- lm(frontier, data = data); names(fr$coefficients)[1] <- c("constant")
                u.p <- c("constant" = 1, rep(0, ncol(sigma_u)-1)); if(length(names(u.p)) > 1){names(u.p)[2:length(u.p)] <- colnames(sigma_u)[2:length(u.p)]}
                w.p <- c("constant" = 1, rep(0, ncol(sigma_w)-1)); if(length(names(w.p)) > 1){names(w.p)[2:length(w.p)] <- colnames(sigma_w)[2:length(w.p)]}
                res.p <- rep(0, p); names(res.p) <- paste("residual", seq(1:p), sep = ".")
                
                z <- c(fr$coefficients, u.p, w.p, res.p)
                names(z) <- make.names(names(z), unique=T)
                
                kk <- length(z)
                k1 <- length(coef(fr))
                k2 <- k1 + length(u.p)
                k3 <- k2 + length(w.p)
        }
        ll <- function(z){
                x_p <- z[1:k1]
                u_p <- z[c(k1 + 1):k2]
                w_p <- z[c(k2 + 1):k3]
                res_p <- z[c(k3 + 1):kk]
                
                Y <- X_fr %*% x_p
                U <- sigma_u %*% u_p
                W <- sigma_w %*% w_p
                Residual <- residual %*% res_p
                
                Sigma_u <- exp(U)
                Sigma_w <- exp(W)
                lambda <- sqrt(Sigma_u/Sigma_w)
                sigma <- Sigma_u + Sigma_w
                
                e <- Y_fr - Y - Residual
                
                ly.x.half <- 0.5 * log(2 / pi) - 0.5 * log(sigma) + pnorm(-lambda * e / sqrt(sigma), log = TRUE) - 0.5 * e^2 / sigma
                
                val <- sum(ly.x.half)
                
                return(-val)
        }
        G <- function(z){
                x_p <- z[1:k1]
                u_p <- z[c(k1 + 1):k2]
                w_p <- z[c(k2 + 1):k3]
                res_p <- z[c(k3 + 1):kk]
                
                Y <- X_fr %*% x_p
                U <- sigma_u %*% u_p
                W <- sigma_w %*% w_p
                Residual <- residual %*% res_p
                
                Sigma_u <- exp(U)
                Sigma_w <- exp(W)
                lambda <- sqrt(Sigma_u/Sigma_w)
                sigma <- Sigma_u + Sigma_w
                
                e <- Y_fr - Y - Residual
                
                zz <- - e * lambda / sqrt(sigma)
                dz <- pmax(dnorm(zz), 9.88131291682493e-324)
                pz <- pmax(pnorm(zz), 9.88131291682493e-324)
                fdp <- dz; cdf <- pz; fdp_cdf <- fdp / cdf
                valor <- as.vector(e / sigma + lambda / sqrt(sigma) * fdp_cdf)
                
                g.b <- t(X_fr) %*% valor
                g.u <- t(sigma_u) %*% as.vector((0.5 / sigma * (e^2 / sigma - e / (lambda * sqrt(sigma)) * fdp_cdf - 1)) * Sigma_u)
                g.w <- t(sigma_w) %*% as.vector((0.5 / sigma * (e^2 / sigma + e * lambda * (2 + lambda^2) * fdp_cdf / sqrt(sigma) - 1)) * Sigma_w)
                g.eta <- t(residual) %*% valor
                
                g <- c(g.b, g.u, g.w, g.eta)
                
                return(-g)
        }
        est <- optim(par = z, fn = ll, gr = G, hessian = TRUE, method = "BFGS",
                     control = list(fnscale = 1, trace = TRUE, REPORT = 1, maxit = 150000))
        #round(data.frame("numerical" = grad(ll,est$par), "analytical" = G(est$par)), 4)
        {
                Z <- xend
                g11 <- lapply(1:p, function(x) cbind(-2 * Z * residual[,x]))
                G11 <- matrix(unlist(g11), nrow = dim(Z)[1L], ncol = dim(Z)[2L]*p)
                
                b <- est$par
                Y_est <- X_fr %*% b[1:k1]
                s2u <- exp(sigma_u %*% b[c(k1 + 1):k2])
                s2w <- exp(sigma_w %*% b[c(k2 + 1):k3])
                Residuals <- residual %*% b[c(k3 + 1):kk]
                sigma <- s2u + s2w
                lambda <- sqrt(s2u/s2w)
                yest <- Y_est + Residuals
                erro <- Y_fr - yest
                
                G2 <- function(z){
                        x_p <- z[1:k1]
                        u_p <- z[c(k1 + 1):k2]
                        w_p <- z[c(k2 + 1):k3]
                        res_p <- z[c(k3 + 1):kk]
                        
                        Y <- X_fr %*% x_p
                        U <- sigma_u %*% u_p
                        W <- sigma_w %*% w_p
                        Residual <- residual %*% res_p
                        
                        Sigma_u <- exp(U)
                        Sigma_w <- exp(W)
                        lambda <- sqrt(Sigma_u / Sigma_w)
                        sigma <- Sigma_u + Sigma_w
                        
                        e <- Y_fr - Y - Residual
                        
                        zz <- - e * lambda / sqrt(sigma)
                        dz <- pmax(dnorm(zz), 9.88131291682493e-324)
                        pz <- pmax(pnorm(zz), 9.88131291682493e-324)
                        fdp <- dz; cdf <- pz; fdp_cdf <- fdp / cdf
                        valor <- as.vector(e / sigma + lambda / sqrt(sigma) * fdp_cdf)
                        
                        g.b <- X_fr * valor
                        g.u <- sigma_u * as.vector((0.5 / sigma * (e ^ 2 / sigma - e / (lambda * sqrt(sigma)) * fdp_cdf - 1)) * Sigma_u)
                        g.w <- sigma_w * as.vector((0.5 / sigma * (e ^ 2 / sigma + e * lambda * (2 + lambda ^ 2) * fdp_cdf / sqrt(sigma) - 1)) * Sigma_w)
                        g.eta <- residual * valor # c(0,0)
                        
                        g <- cbind(g.b, g.u, g.w, g.eta)
                        
                        return(g)
                }
                G22 <- G2(b)
                
                g21 <- function(z){
                        zz <- - erro * lambda / sqrt(sigma)
                        dz <- pmax(dnorm(zz), 9.88131291682493e-324)
                        pz <- pmax(pnorm(zz), 9.88131291682493e-324)
                        fdp <- dz; cdf <- pz; fdp_cdf <- fdp / cdf
                        valor <- as.vector(erro / sigma + lambda / sqrt(sigma) * fdp_cdf)
                        eta <- z[c(k3 + 1):kk]
                        
                        g0 <- lapply(1:p, function(x) cbind(- Z * valor * eta[x]))
                        g <- matrix(unlist(g0), nrow = dim(Z)[1L], ncol = dim(Z)[2L]*p)
                        return(g)
                }
                G21 <- g21(b)
                
                V1 <- vcov(ols); V2 <- ginv(est$hessian); C <- t(G22) %*% G21; R <- t(G22) %*% G11
                V2_adj <- V2 + V2 %*% ((C %*% V1 %*% t(C)) - (R %*% V1 %*% t(C)) - (C %*% V1 %*% t(R))) %*% V2
                
                ep <- sqrt(diag(V2_adj))
                estZ <- lapply(length(b), function(x) b/ep)
                p_valor <- lapply(length(b), function(x) (1 - pnorm(abs(b/ep)))*2)
                ics <- lapply(length(b), function(x) cbind(b - qnorm(.975) * ep, b + qnorm(.975) * ep))
                result <- data.frame("Coefficient" = cbind(est$par), "Std.Err" = ep, "z" = estZ[[1]], 
                                     "P-value" = p_valor[[1]], "Lower limit" = ics[[1]][,1], "Upper limit" = ics[[1]][,2])
                
                mu.mod <- - erro * s2u / sigma
                s.mod <- sqrt(s2w * s2u / sigma)
                uf <- mu.mod + s.mod * ((dnorm(-mu.mod / s.mod)) / (pnorm(mu.mod / s.mod)))
                ef <- ((pnorm(-s.mod + mu.mod / s.mod)) / (pnorm(mu.mod / s.mod))) * exp(- mu.mod + 0.5 * s.mod^2)
                
                lns2u <- est$par[c(k1 + 1):k2]
                lns2w <- est$par[c(k2 + 1):k3]
                if (length(lns2u) == 1 & length(lns2w) == 1){
                        names(lns2u) <- c('lns2u')
                        cov.u <- ginv(est$hessian)[c(k1 + 1):k2,c(k1 + 1):k2]
                        S2U <- deltaMethod(lns2u, "exp(lns2u)", vcov = cov.u)
                        zS2U <- S2U$Estimate/S2U$SE
                        pS2U <- (1 - pnorm(abs(zS2U)))*2
                        icS2U <- lapply(length(S2U$Estimate), function(x) cbind(S2U$Estimate - qnorm(.975) * S2U$SE, S2U$Estimate + qnorm(.975) * S2U$SE))
                        S2U0 <- data.frame(S2U$Estimate, S2U$SE, zS2U, pS2U, icS2U[[1]][,1], icS2U[[1]][,2])
                        names(S2U0) <- names(result)
                        rownames(S2U0) <- "s2u"
                        
                        names(lns2w) <- c('lns2w')
                        cov.w <- ginv(est$hessian)[c(k2 + 1):k3,c(k2 + 1):k3]
                        S2W <- deltaMethod(lns2w, "exp(lns2w)", vcov = cov.w)
                        zS2W <- S2W$Estimate/S2W$SE
                        pS2W <- (1 - pnorm(abs(zS2W)))*2
                        icS2W <- lapply(length(S2W$Estimate), function(x) cbind(S2W$Estimate - qnorm(.975) * S2W$SE, S2W$Estimate + qnorm(.975) * S2W$SE))
                        S2W0 <- data.frame(S2W$Estimate, S2W$SE, zS2W, pS2W, icS2W[[1]][,1], icS2W[[1]][,2])
                        names(S2W0) <- names(result)
                        rownames(S2W0) <- "s2w"
                        
                        s20 <- c(S2U$Estimate, S2W$Estimate)
                        names(s20) <- c('s2u','s2w')
                        cov.s <- ginv(est$hessian)[c(k1 + 1):k3,c(k1 + 1):k3]
                        S2 <- deltaMethod(s20, "s2u + s2w", vcov = cov.s)
                        zS2 <- S2$Estimate/S2$SE
                        pS2 <- (1 - pnorm(abs(zS2)))*2
                        icS2 <- lapply(length(S2$Estimate), function(x) cbind(S2$Estimate - qnorm(.975) * S2$SE, S2$Estimate + qnorm(.975) * S2$SE))
                        S20 <- data.frame(S2$Estimate, S2$SE, zS2, pS2, icS2[[1]][,1], icS2[[1]][,2])
                        names(S20) <- names(result)
                        rownames(S20) <- "s2"
                        
                        Lambda <- deltaMethod(s20, "sqrt(s2u/s2w)", vcov = cov.s)
                        zLam <- Lambda$Estimate/Lambda$SE
                        pLam <- (1 - pnorm(abs(zLam)))*2
                        icLam <- lapply(length(Lambda$Estimate), function(x) cbind(Lambda$Estimate - qnorm(.975) * Lambda$SE, Lambda$Estimate + qnorm(.975) * Lambda$SE))
                        LAM <- data.frame(Lambda$Estimate, Lambda$SE, zLam, pLam, icLam[[1]][,1], icLam[[1]][,2])
                        names(LAM) <- names(result)
                        rownames(LAM) <- "lambda"
                        
                        RESP <- data.frame(rbind(result, S2U0, S2W0, S20, LAM))
                } else if (length(lns2u) == 1){
                        lns2u <- est$par[c(k1 + 1):k2]
                        names(lns2u) <- c('lns2u')
                        cov.u <- ginv(est$hessian)[c(k1 + 1):k2,c(k1 + 1):k2]
                        S2U <- deltaMethod(lns2u, "exp(lns2u)", vcov = cov.u)
                        zS2U <- S2U$Estimate/S2U$SE
                        pS2U <- (1 - pnorm(abs(zS2U)))*2
                        icS2U <- lapply(length(S2U$Estimate), function(x) cbind(S2U$Estimate - qnorm(.975) * S2U$SE, S2U$Estimate + qnorm(.975) * S2U$SE))
                        S2U0 <- data.frame(S2U$Estimate, S2U$SE, zS2U, pS2U, icS2U[[1]][,1], icS2U[[1]][,2])
                        names(S2U0) <- names(result)
                        rownames(S2U0) <- "s2u"
                        
                        RESP <- data.frame(rbind(result, S2U0))
                } else if (length(lns2w) == 1){
                        lns2w <- est$par[c(k2 + 1):k3]
                        names(lns2w) <- c('lns2w')
                        cov.w <- ginv(est$hessian)[c(k2 + 1):k3,c(k2 + 1):k3]
                        S2W <- deltaMethod(lns2w, "exp(lns2w)", vcov = cov.w)
                        zS2W <- S2W$Estimate/S2W$SE
                        pS2W <- (1 - pnorm(abs(zS2W)))*2
                        icS2W <- lapply(length(S2W$Estimate), function(x) cbind(S2W$Estimate - qnorm(.975) * S2W$SE, S2W$Estimate + qnorm(.975) * S2W$SE))
                        S2W0 <- data.frame(S2W$Estimate, S2W$SE, zS2W, pS2W, icS2W[[1]][,1], icS2W[[1]][,2])
                        names(S2W0) <- names(result)
                        rownames(S2W0) <- "s2w"
                        
                        RESP <- data.frame(rbind(result, S2W0))
                } else {
                        RESP <- result
                }
                
                # LR test for endogeneity
                lnL1 <- -ll(est$par)
                lnL00 <- function(z) {
                        x_p <- z[1:k1]
                        u_p <- z[c(k1 + 1):k2]
                        w_p <- z[c(k2 + 1):k3]
                        res_p <- rep(0,p)
                        
                        Y <- X_fr %*% x_p
                        U <- sigma_u %*% u_p
                        W <- sigma_w %*% w_p
                        Residual <- residual %*% res_p
                        
                        Sigma_u <- exp(U)
                        Sigma_w <- exp(W)
                        lambda <- sqrt(Sigma_u/Sigma_w)
                        sigma <- Sigma_u + Sigma_w
                        
                        e <- Y_fr - Y - Residual
                        
                        ly.x.half <- 0.5 * log(2 / pi) - 0.5 * log(sigma) + pnorm(-lambda * e / sqrt(sigma), log = TRUE) - 0.5 * e^2 / sigma
                        
                        val <- sum(ly.x.half)
                        
                        return(-val)
                }
                lnL0 <- -lnL00(est$par)
                LRtest <- 2 * (lnL1 - lnL0)
                LRp <- pchisq(LRtest, df = p, lower.tail = FALSE)
                
                # Wald test for endogeneity
                etas <- b[c(k3 + 1):kk]
                cov_etas <- V2_adj[c(k3 + 1):kk, c(k3 + 1):kk]
                waldtest <- t(etas) %*% ginv(cov_etas) %*% etas
                waldp <- pchisq(waldtest, df = p, lower.tail = FALSE)
                
                bias <- mean(Y_est) - mean(Y_fr)
                RMSE <- sqrt(var(Y_est) + bias^2)
                pearson <- cor(Y_est, Y_fr, method = "pearson")
                spearman <- cor(Y_est, Y_fr, method = "spearman")
                
                n <- length(Y_fr)
                W <- crossprod(residual) / n
                ## Likelihood of x
                lx <- n * 0.5 * (- p * log(2 * pi) - log(det(W))) - 0.5 * sum(mahalanobis(residual, 0, W))
                loglik <- lx - est$value
                lnL <- loglik
                K <- length(b)
                AIC <- - 2 * lnL + 2 * K
                BIC <- - 2 * lnL + log(n) * K
        }
        lista <- list("efficiency" = ef, "error" = erro, "fitted.y_without.correction" = Y_est,
                      "fitted.y_with.correction" = yest, "reg_IV" = result_IV, "table" = RESP,
                      "summary.ef" = summary(ef), "sd.ef" = sd(ef, na.rm = TRUE),
                      "value" = loglik, "AIC" = AIC, "BIC" = BIC, 
                      "cor pearson IV" = cor_p.iv, "cor spearman IV" = cor_s.iv,
                      "cor.pearson" = c(pearson), "cor.spearman" = c(spearman),
                      "LR chisq test" = matrix(cbind(LRtest, LRp), 1, 2, dimnames = list(c(),c("statistic", "P-value"))),
                      "Wald chisq test" = matrix(cbind(waldtest, waldp), 1, 2, dimnames = list(c(),c("statistic", "P-value"))),
                      "bias" = c(bias), "RMSE" = c(RMSE), 
                      "sample size" = n, "estimated parameters" = K)
        return(lista)
}

SF.exp2S <- function(fr.form, end.form, s2u.form, s2w.form, data = sys.parent()) {
        {# Frontier
                frontier <- fr.form
                mf <- model.frame(frontier, data)
                X_fr <- model.matrix(frontier, mf)
                Y_fr <- model.response(mf, "numeric")
                
                # Endogenous variables
                End <- end.form
                m <- model.frame(End, data)
                xend <- model.matrix(End, m)
                yend <- model.response(m, "numeric")
                p <- ncol(as.matrix(yend))
                
                # Sigma_u and Sigma_w
                Sigma_u <- s2u.form
                sig_u <- model.frame(Sigma_u, data)
                sigma_u <- model.matrix(Sigma_u, sig_u)
                
                Sigma_w <- s2w.form
                sig_w <- model.frame(Sigma_w, data)
                sigma_w <- model.matrix(Sigma_w, sig_w)
                
                # Step 1 - OLS of the Endogenous Variables
                ols <- lm(End, data); if(length(names(ols$coefficients)) > 1){names(ols$coefficients)[1] <- c("constant")} else{rownames(ols$coefficients)[1] <- c("constant")}
                residual <- as.matrix(residuals(ols))
                
                yend_est <- as.matrix(fitted(ols))
                if(ncol(yend_est) > 1){colnames(yend_est) <- paste(colnames(yend), "IV", sep = ".")}
                cor_p.iv <- cor(yend, yend_est, method = "pearson")
                cor_s.iv <- cor(yend, yend_est, method = "spearman")
                
                result_IV <- summary(ols)
                
                # Initial values
                fr <- lm(frontier, data = data); names(fr$coefficients)[1] <- c("constant")
                u.p <- c("constant" = 1, rep(0, ncol(sigma_u)-1)); if(length(names(u.p)) > 1){names(u.p)[2:length(u.p)] <- colnames(sigma_u)[2:length(u.p)]}
                w.p <- c("constant" = 1, rep(0, ncol(sigma_w)-1)); if(length(names(w.p)) > 1){names(w.p)[2:length(w.p)] <- colnames(sigma_w)[2:length(w.p)]}
                res.p <- rep(0, p); names(res.p) <- paste("residual", seq(1:p), sep = ".")
                
                z <- c(fr$coefficients, u.p, w.p, res.p)
                names(z) <- make.names(names(z), unique=T)
                
                kk <- length(z)
                k1 <- length(coef(fr))
                k2 <- k1 + length(u.p)
                k3 <- k2 + length(w.p)
        }
        ll <- function(z){
                x_p <- z[1:k1]
                u_p <- z[c(k1 + 1):k2]
                w_p <- z[c(k2 + 1):k3]
                res_p <- z[c(k3 + 1):kk]
                
                Y <- X_fr %*% x_p
                U <- sigma_u %*% u_p
                W <- sigma_w %*% w_p
                Residual <- residual %*% res_p
                
                Sigma_u <- exp(U)
                Sigma_w <- exp(W)
                lambda <- sqrt(Sigma_u/Sigma_w)
                sigma <- Sigma_u + Sigma_w
                
                e <- Y_fr - Y - Residual
                
                zz <- (-e - Sigma_w / sqrt(Sigma_u)) / sqrt(Sigma_w)
                pz <- pnorm(zz, log = TRUE)
                ly.x.exp <- -0.5 * log(Sigma_u) + 0.5 * Sigma_w / Sigma_u + pz + e / sqrt(Sigma_u)
                
                val <- sum(ly.x.exp)
                
                return(-val)
        }
        G <- function(z){
                x_p <- z[1:k1]
                u_p <- z[c(k1 + 1):k2]
                w_p <- z[c(k2 + 1):k3]
                res_p <- z[c(k3 + 1):kk]
                
                Y <- X_fr %*% x_p
                U <- sigma_u %*% u_p
                W <- sigma_w %*% w_p
                Residual <- residual %*% res_p
                
                Sigma_u <- exp(U)
                Sigma_w <- exp(W)
                lambda <- sqrt(Sigma_u/Sigma_w)
                sigma <- Sigma_u + Sigma_w
                
                e <- Y_fr - Y - Residual
                
                zz <- (-e - Sigma_w / sqrt(Sigma_u)) / sqrt(Sigma_w)
                dz <- pmax(dnorm(zz), 9.88131291682493e-324)
                pz <- pmax(pnorm(zz), 9.88131291682493e-324)
                fdp <- dz; cdf <- pz; fdp_cdf <- fdp / cdf
                valor <- as.vector(fdp_cdf / sqrt(Sigma_w) - 1 / sqrt(Sigma_u))
                
                g.b <- t(X_fr) %*% valor
                g.u <- t(sigma_u) %*% as.vector((0.5 / Sigma_u * (fdp_cdf * sqrt(Sigma_w / Sigma_u) - Sigma_w / Sigma_u - e / sqrt(Sigma_u) - 1)) * Sigma_u)
                g.w <- t(sigma_w) %*% as.vector((0.5 * (1 / Sigma_u + fdp_cdf / sqrt(Sigma_w) * (e / Sigma_w - 1 / sqrt(Sigma_u)))) * Sigma_w)
                g.eta <- t(residual) %*% valor
                
                g <- c(g.b, g.u, g.w, g.eta)
                
                return(-g)
        }
        est <- optim(par = z, fn = ll, gr = G, hessian = TRUE, method = "Nelder-Mead", 
                     control = list(fnscale = 1, trace = TRUE, REPORT = 1, maxit = 150000))
        #round(data.frame("numerical" = grad(ll,est$par), "analytical" = G(est$par)), 4)
        {
                Z <- xend
                g11 <- lapply(1:p, function(x) cbind(-2 * Z * residual[,x]))
                G11 <- matrix(unlist(g11), nrow = dim(Z)[1L], ncol = dim(Z)[2L]*p)
                
                b <- est$par
                Y_est <- X_fr %*% b[1:k1]
                s2u <- exp(sigma_u %*% b[c(k1 + 1):k2])
                s2w <- exp(sigma_w %*% b[c(k2 + 1):k3])
                Residuals <- residual %*% b[c(k3 + 1):kk]
                sigma <- s2u + s2w
                lambda <- sqrt(s2u/s2w)
                yest <- Y_est + Residuals
                erro <- Y_fr - yest
                
                G2 <- function(z){
                        x_p <- z[1:k1]
                        u_p <- z[c(k1 + 1):k2]
                        w_p <- z[c(k2 + 1):k3]
                        res_p <- z[c(k3 + 1):kk]
                        
                        Y <- X_fr %*% x_p
                        U <- sigma_u %*% u_p
                        W <- sigma_w %*% w_p
                        Residual <- residual %*% res_p
                        
                        Sigma_u <- exp(U)
                        Sigma_w <- exp(W)
                        lambda <- sqrt(Sigma_u/Sigma_w)
                        sigma <- Sigma_u + Sigma_w
                        
                        e <- Y_fr - Y - Residual
                        
                        zz <- (-e - Sigma_w / sqrt(Sigma_u)) / sqrt(Sigma_w)
                        dz <- pmax(dnorm(zz), 9.88131291682493e-324)
                        pz <- pmax(pnorm(zz), 9.88131291682493e-324)
                        fdp <- dz; cdf <- pz; fdp_cdf <- fdp / cdf
                        valor <- as.vector(fdp_cdf / sqrt(Sigma_w) - 1 / sqrt(Sigma_u))
                        
                        g.b <- X_fr * valor
                        g.u <- sigma_u * as.vector((0.5 / Sigma_u * (fdp_cdf * sqrt(Sigma_w / Sigma_u) - Sigma_w / Sigma_u - e / sqrt(Sigma_u) - 1)) * Sigma_u)
                        g.w <- sigma_w * as.vector((0.5 * (1 / Sigma_u + fdp_cdf / sqrt(Sigma_w) * (e / Sigma_w - 1 / sqrt(Sigma_u)))) * Sigma_w)
                        g.eta <- residual * valor
                        
                        g <- cbind(g.b, g.u, g.w, g.eta)
                        
                        return(g)
                }
                G22 <- G2(b)
                
                g21 <- function(z){
                        zz <- (-erro - s2w / sqrt(s2u)) / sqrt(s2w)
                        dz <- pmax(dnorm(zz), 9.88131291682493e-324)
                        pz <- pmax(pnorm(zz), 9.88131291682493e-324)
                        fdp <- dz; cdf <- pz; fdp_cdf <- fdp / cdf
                        valor <- as.vector(fdp_cdf / sqrt(s2w) - 1 / sqrt(s2u))
                        eta <- z[c(k3 + 1):kk]
                        
                        g0 <- lapply(1:p, function(x) cbind(- Z * valor * eta[x]))
                        g <- matrix(unlist(g0), nrow = dim(Z)[1L], ncol = dim(Z)[2L]*p)
                        return(g)
                }
                G21 <- g21(b)
                
                V1 <- vcov(ols); V2 <- ginv(est$hessian); C <- t(G22) %*% G21; R <- t(G22) %*% G11
                V2_adj <- V2 + V2 %*% ((C %*% V1 %*% t(C)) - (R %*% V1 %*% t(C)) - (C %*% V1 %*% t(R))) %*% V2
                
                ep <- sqrt(diag(V2_adj))
                estZ <- lapply(length(b), function(x) b/ep)
                p_valor <- lapply(length(b), function(x) (1 - pnorm(abs(b/ep)))*2)
                ics <- lapply(length(b), function(x) cbind(b - qnorm(.975) * ep, b + qnorm(.975) * ep))
                result <- data.frame("Coefficient" = cbind(est$par), "Std.Err" = ep, "z" = estZ[[1]], 
                                     "P-value" = p_valor[[1]], "Lower limit" = ics[[1]][,1], "Upper limit" = ics[[1]][,2])
                
                mu.mod <- - erro - (s2w/sqrt(s2u))
                s.mod <- sqrt(s2w)
                uf <- mu.mod + s.mod * ((dnorm(-mu.mod / s.mod)) / (pnorm(mu.mod / s.mod)))
                ef <- ((pnorm(-s.mod + mu.mod / s.mod)) / (pnorm(mu.mod / s.mod))) * exp(- mu.mod + 0.5 * s.mod^2)
                
                lns2u <- est$par[c(k1 + 1):k2]
                lns2w <- est$par[c(k2 + 1):k3]
                if (length(lns2u) == 1 & length(lns2w) == 1){
                        names(lns2u) <- c('lns2u')
                        cov.u <- ginv(est$hessian)[c(k1 + 1):k2,c(k1 + 1):k2]
                        S2U <- deltaMethod(lns2u, "exp(lns2u)", vcov = cov.u)
                        zS2U <- S2U$Estimate/S2U$SE
                        pS2U <- (1 - pnorm(abs(zS2U)))*2
                        icS2U <- lapply(length(S2U$Estimate), function(x) cbind(S2U$Estimate - qnorm(.975) * S2U$SE, S2U$Estimate + qnorm(.975) * S2U$SE))
                        S2U0 <- data.frame(S2U$Estimate, S2U$SE, zS2U, pS2U, icS2U[[1]][,1], icS2U[[1]][,2])
                        names(S2U0) <- names(result)
                        rownames(S2U0) <- "s2u"
                        
                        names(lns2w) <- c('lns2w')
                        cov.w <- ginv(est$hessian)[c(k2 + 1):k3,c(k2 + 1):k3]
                        S2W <- deltaMethod(lns2w, "exp(lns2w)", vcov = cov.w)
                        zS2W <- S2W$Estimate/S2W$SE
                        pS2W <- (1 - pnorm(abs(zS2W)))*2
                        icS2W <- lapply(length(S2W$Estimate), function(x) cbind(S2W$Estimate - qnorm(.975) * S2W$SE, S2W$Estimate + qnorm(.975) * S2W$SE))
                        S2W0 <- data.frame(S2W$Estimate, S2W$SE, zS2W, pS2W, icS2W[[1]][,1], icS2W[[1]][,2])
                        names(S2W0) <- names(result)
                        rownames(S2W0) <- "s2w"
                        
                        s20 <- c(S2U$Estimate, S2W$Estimate)
                        names(s20) <- c('s2u','s2w')
                        cov.s <- ginv(est$hessian)[c(k1 + 1):k3,c(k1 + 1):k3]
                        S2 <- deltaMethod(s20, "s2u + s2w", vcov = cov.s)
                        zS2 <- S2$Estimate/S2$SE
                        pS2 <- (1 - pnorm(abs(zS2)))*2
                        icS2 <- lapply(length(S2$Estimate), function(x) cbind(S2$Estimate - qnorm(.975) * S2$SE, S2$Estimate + qnorm(.975) * S2$SE))
                        S20 <- data.frame(S2$Estimate, S2$SE, zS2, pS2, icS2[[1]][,1], icS2[[1]][,2])
                        names(S20) <- names(result)
                        rownames(S20) <- "s2"
                        
                        Lambda <- deltaMethod(s20, "sqrt(s2u/s2w)", vcov = cov.s)
                        zLam <- Lambda$Estimate/Lambda$SE
                        pLam <- (1 - pnorm(abs(zLam)))*2
                        icLam <- lapply(length(Lambda$Estimate), function(x) cbind(Lambda$Estimate - qnorm(.975) * Lambda$SE, Lambda$Estimate + qnorm(.975) * Lambda$SE))
                        LAM <- data.frame(Lambda$Estimate, Lambda$SE, zLam, pLam, icLam[[1]][,1], icLam[[1]][,2])
                        names(LAM) <- names(result)
                        rownames(LAM) <- "lambda"
                        
                        RESP <- data.frame(rbind(result, S2U0, S2W0, S20, LAM))
                } else if (length(lns2u) == 1){
                        lns2u <- est$par[c(k1 + 1):k2]
                        names(lns2u) <- c('lns2u')
                        cov.u <- ginv(est$hessian)[c(k1 + 1):k2,c(k1 + 1):k2]
                        S2U <- deltaMethod(lns2u, "exp(lns2u)", vcov = cov.u)
                        zS2U <- S2U$Estimate/S2U$SE
                        pS2U <- (1 - pnorm(abs(zS2U)))*2
                        icS2U <- lapply(length(S2U$Estimate), function(x) cbind(S2U$Estimate - qnorm(.975) * S2U$SE, S2U$Estimate + qnorm(.975) * S2U$SE))
                        S2U0 <- data.frame(S2U$Estimate, S2U$SE, zS2U, pS2U, icS2U[[1]][,1], icS2U[[1]][,2])
                        names(S2U0) <- names(result)
                        rownames(S2U0) <- "s2u"
                        
                        RESP <- data.frame(rbind(result, S2U0))
                } else if (length(lns2w) == 1){
                        lns2w <- est$par[c(k2 + 1):k3]
                        names(lns2w) <- c('lns2w')
                        cov.w <- ginv(est$hessian)[c(k2 + 1):k3,c(k2 + 1):k3]
                        S2W <- deltaMethod(lns2w, "exp(lns2w)", vcov = cov.w)
                        zS2W <- S2W$Estimate/S2W$SE
                        pS2W <- (1 - pnorm(abs(zS2W)))*2
                        icS2W <- lapply(length(S2W$Estimate), function(x) cbind(S2W$Estimate - qnorm(.975) * S2W$SE, S2W$Estimate + qnorm(.975) * S2W$SE))
                        S2W0 <- data.frame(S2W$Estimate, S2W$SE, zS2W, pS2W, icS2W[[1]][,1], icS2W[[1]][,2])
                        names(S2W0) <- names(result)
                        rownames(S2W0) <- "s2w"
                        
                        RESP <- data.frame(rbind(result, S2W0))
                } else {
                        RESP <- result
                }
                
                # LR test for endogeneity
                lnL1 <- -ll(est$par)
                lnL00 <- function(z) {
                        x_p <- z[1:k1]
                        u_p <- z[c(k1 + 1):k2]
                        w_p <- z[c(k2 + 1):k3]
                        res_p <- rep(0,p)
                        
                        Y <- X_fr %*% x_p
                        U <- sigma_u %*% u_p
                        W <- sigma_w %*% w_p
                        Residual <- residual %*% res_p
                        
                        Sigma_u <- exp(U)
                        Sigma_w <- exp(W)
                        lambda <- sqrt(Sigma_u/Sigma_w)
                        sigma <- Sigma_u + Sigma_w
                        
                        e <- Y_fr - Y - Residual
                        
                        zz <- (-e - Sigma_w / sqrt(Sigma_u)) / sqrt(Sigma_w)
                        pz <- pnorm(zz, log = TRUE)
                        ly.x.exp <- -0.5 * log(Sigma_u) + 0.5 * Sigma_w / Sigma_u + pz + e / sqrt(Sigma_u)
                        
                        val <- sum(ly.x.exp)
                        
                        return(-val)
                }
                lnL0 <- -lnL00(est$par)
                LRtest <- 2 * (lnL1 - lnL0)
                LRp <- pchisq(LRtest, df = p, lower.tail = FALSE)
                
                # Wald test for endogeneity
                etas <- b[c(k3 + 1):kk]
                cov_etas <- V2_adj[c(k3 + 1):kk, c(k3 + 1):kk]
                waldtest <- t(etas) %*% ginv(cov_etas) %*% etas
                waldp <- pchisq(waldtest, df = p, lower.tail = FALSE)
                
                bias <- mean(Y_est) - mean(Y_fr)
                RMSE <- sqrt(var(Y_est) + bias^2)
                pearson <- cor(Y_est, Y_fr, method = "pearson")
                spearman <- cor(Y_est, Y_fr, method = "spearman")
                
                n <- length(Y_fr)
                W <- crossprod(residual) / n
                ## Likelihood of x
                lx <- n * 0.5 * (- p * log(2 * pi) - log(det(W))) - 0.5 * sum(mahalanobis(residual, 0, W))
                loglik <- lx - est$value
                lnL <- loglik
                K <- length(b)
                AIC <- - 2 * lnL + 2 * K
                BIC <- - 2 * lnL + log(n) * K
        }
        lista <- list("efficiency" = ef, "error" = erro, "fitted.y_without.correction" = Y_est,
                      "fitted.y_with.correction" = yest, "reg_IV" = result_IV, "table" = RESP,
                      "summary.ef" = summary(ef), "sd.ef" = sd(ef, na.rm = TRUE),
                      "value" = loglik, "AIC" = AIC, "BIC" = BIC, 
                      "cor pearson IV" = cor_p.iv, "cor spearman IV" = cor_s.iv,
                      "cor.pearson" = c(pearson), "cor.spearman" = c(spearman),
                      "LR chisq test" = matrix(cbind(LRtest, LRp), 1, 2, dimnames = list(c(),c("statistic", "P-value"))),
                      "Wald chisq test" = matrix(cbind(waldtest, waldp), 1, 2, dimnames = list(c(),c("statistic", "P-value"))),
                      "bias" = c(bias), "RMSE" = c(RMSE), 
                      "sample size" = n, "estimated parameters" = K)
        return(lista)
}

SF.trunc2S <- function(fr.form, end.form, mu.form, s2u.form = ~1, s2w.form = ~1, data = sys.parent()) {
        {# Frontier
                frontier <- fr.form
                mf <- model.frame(frontier, data)
                X_fr <- model.matrix(frontier, mf)
                Y_fr <- model.response(mf, "numeric")
                
                # Endogenous variables
                End <- end.form
                m <- model.frame(End, data)
                xend <- model.matrix(End, m)
                yend <- model.response(m, "numeric")
                p <- ncol(as.matrix(yend))
                
                Mu <- mu.form
                Mi <- model.frame(Mu, data)
                mi <- model.matrix(Mu, Mi)
                
                # Sigma_u and Sigma_w
                Sigma_u <- s2u.form
                sig_u <- model.frame(Sigma_u, data)
                sigma_u <- model.matrix(Sigma_u, sig_u)
                
                Sigma_w <- s2w.form
                sig_w <- model.frame(Sigma_w, data)
                sigma_w <- model.matrix(Sigma_w, sig_w)
                
                # Step 1 - OLS of the Endogenous Variables
                ols <- lm(End, data); if(length(names(ols$coefficients)) > 1){names(ols$coefficients)[1] <- c("constant")} else{rownames(ols$coefficients)[1] <- c("constant")}
                residual <- as.matrix(residuals(ols))
                
                yend_est <- as.matrix(fitted(ols))
                if(ncol(yend_est) > 1){colnames(yend_est) <- paste(colnames(yend), "IV", sep = ".")}
                cor_p.iv <- cor(yend, yend_est, method = "pearson")
                cor_s.iv <- cor(yend, yend_est, method = "spearman")
                
                result_IV <- summary(ols)
                
                # Initial values
                fr <- lm(frontier, data = data); names(fr$coefficients)[1] <- c("constant")
                mu.p <- c("constant" = 1, rep(0, ncol(mi) - 1)); if(length(names(mu.p)) > 1){names(mu.p)[2:length(mu.p)] <- colnames(mi)[2:length(mu.p)]}
                u.p <- c("constant" = 1)
                w.p <- c("constant" = 1)
                res.p <- rep(0, p); names(res.p) <- paste("residual", seq(1:p), sep = ".")
                
                z <- c(fr$coefficients, mu.p, u.p, w.p, res.p)
                names(z) <- make.names(names(z), unique=T)
                
                kk <- length(z)
                k1 <- length(coef(fr))
                k2 <- k1 + length(mu.p)
                k3 <- k2 + length(u.p)
                k4 <- k3 + length(w.p)
        }
        ll <- function(z){
                x_p <- z[1:k1]
                mu_p <- z[c(k1 + 1):k2]
                u_p <- z[c(k2 + 1):k3]
                w_p <- z[c(k3 + 1):k4]
                res_p <- z[c(k4 + 1):kk]
                
                Y <- X_fr %*% x_p
                mu <- mi %*% mu_p
                U <- sigma_u %*% u_p
                W <- sigma_w %*% w_p
                Residual <- residual %*% res_p
                
                Sigma_u <- exp(U)
                Sigma_w <- exp(W)
                lambda <- sqrt(Sigma_u/Sigma_w)
                sigma <- Sigma_u + Sigma_w
                
                e <- Y_fr - Y - Residual
                
                aa <- mu / sqrt(sigma) * sqrt(1 + lambda^-2)
                pa <- pmax(pnorm(aa), 9.88131291682493e-324)
                bb <- (mu / lambda - e * lambda) / sqrt(sigma)
                pb <- pmax(pnorm(bb), 9.88131291682493e-324)
                cc <- (e + mu)^2 / sigma
                
                ly.x.trun <- -0.5 * log(2 * pi) - 0.5 * log(sigma) - log(pa) + log(pb) - 0.5 * cc
                
                val <- sum(ly.x.trun)
                
                return(-val)
        }
        G <- function(z){
                x_p <- z[1:k1]
                mu_p <- z[c(k1 + 1):k2]
                u_p <- z[c(k2 + 1):k3]
                w_p <- z[c(k3 + 1):k4]
                res_p <- z[c(k4 + 1):kk]
                
                Y <- X_fr %*% x_p
                mu <- mi %*% mu_p
                U <- sigma_u %*% u_p
                W <- sigma_w %*% w_p
                Residual <- residual %*% res_p
                
                Sigma_u <- exp(U)
                Sigma_w <- exp(W)
                lambda <- sqrt(Sigma_u/Sigma_w)
                sigma <- Sigma_u + Sigma_w
                
                e <- Y_fr - Y - Residual
                
                aa <- mu / sqrt(sigma) * sqrt(1 + lambda^-2)
                da <- pmax(dnorm(aa), 9.88131291682493e-324)
                pa <- pmax(pnorm(aa), 9.88131291682493e-324)
                
                bb <- (mu / lambda - e * lambda) / sqrt(sigma)
                db <- pmax(dnorm(bb), 9.88131291682493e-324)
                pb <- pmax(pnorm(bb), 9.88131291682493e-324)
                
                cc <- (e + mu)^2 / sigma
                da_pa <- da/pa; db_pb <- db/pb
                
                valor <- as.vector((e + mu) / sigma + db_pb * lambda / sqrt(sigma))
                
                g.b <- t(X_fr) %*% valor
                
                g.mu <- t(mi) %*% as.vector(1 / sqrt(sigma) * (db_pb / lambda - da_pa * sqrt(1 + lambda^(-2)) - (e + mu) / sqrt(sigma)))
                
                Partu <- (2 * mu + e + mu / lambda^2) * db_pb / (lambda * sqrt(sigma))
                g.u <- t(sigma_u) %*% as.vector((0.5 * mu / Sigma_u^1.5 * da_pa + 0.5 / sigma * (cc - Partu - 1)) * Sigma_u)
                
                Partw <- lambda / sqrt(sigma) * db_pb * (mu + 2 * e + e * lambda^2)
                g.w <- t(sigma_w) %*% as.vector((0.5 / sigma * (cc + Partw - 1)) * Sigma_w)
                
                g.eta <- t(residual) %*% valor
                
                g <- c(g.b, g.mu, g.u, g.w, g.eta)
                
                return(-g)
        }
        est <- optim(par = z, fn = ll, gr = G, hessian = TRUE, method = "BFGS", 
                     control = list(fnscale = 1, trace = TRUE, REPORT = 1, maxit = 150000))
        #round(data.frame("numerical" = grad(ll,est$par), "analytical" = G(est$par)), 4)
        {        
                Z <- xend
                g11 <- lapply(1:p, function(x) cbind(-2 * Z * residual[,x]))
                G11 <- matrix(unlist(g11), nrow = dim(Z)[1L], ncol = dim(Z)[2L]*p)
                
                b <- est$par
                Y_est <- X_fr %*% b[1:k1]
                mu <- mi %*% b[c(k1 + 1):k2]
                s2u <- exp(sigma_u %*% b[c(k2 + 1):k3])
                s2w <- exp(sigma_w %*% b[c(k3 + 1):k4])
                Residuals <- residual %*% b[c(k4 + 1):kk]
                sigma <- s2u + s2w
                lambda <- sqrt(s2u/s2w)
                yest <- Y_est + Residuals
                erro <- Y_fr - yest
                
                G2 <- function(z){
                        x_p <- z[1:k1]
                        mu_p <- z[c(k1 + 1):k2]
                        u_p <- z[c(k2 + 1):k3]
                        w_p <- z[c(k3 + 1):k4]
                        res_p <- z[c(k4 + 1):kk]
                        
                        Y <- X_fr %*% x_p
                        mu <- mi %*% mu_p
                        U <- sigma_u %*% u_p
                        W <- sigma_w %*% w_p
                        Residual <- residual %*% res_p
                        
                        Sigma_u <- exp(U)
                        Sigma_w <- exp(W)
                        lambda <- sqrt(Sigma_u/Sigma_w)
                        sigma <- Sigma_u + Sigma_w
                        
                        e <- Y_fr - Y - Residual
                        
                        aa <- mu / sqrt(sigma) * sqrt(1 + lambda^-2)
                        da <- pmax(dnorm(aa), 9.88131291682493e-324)
                        pa <- pmax(pnorm(aa), 9.88131291682493e-324)
                        
                        bb <- (mu / lambda - e * lambda) / sqrt(sigma)
                        db <- pmax(dnorm(bb), 9.88131291682493e-324)
                        pb <- pmax(pnorm(bb), 9.88131291682493e-324)
                        
                        cc <- (e + mu)^2 / sigma
                        da_pa <- da/pa; db_pb <- db/pb
                        
                        valor <- as.vector((e + mu) / sigma + db_pb * lambda / sqrt(sigma))
                        
                        g.b <- X_fr * valor
                        
                        g.mu <- mi * as.vector(1 / sqrt(sigma) * (db_pb / lambda - da_pa * sqrt(1 + lambda^(-2)) - (e + mu) / sqrt(sigma)))
                        
                        Partu <- (2 * mu + e + mu / lambda^2) * db_pb / (lambda * sqrt(sigma))
                        g.u <- sigma_u * as.vector((0.5 * mu / Sigma_u^1.5 * da_pa + 0.5 / sigma * (cc - Partu - 1)) * Sigma_u)
                        
                        Partw <- lambda / sqrt(sigma) * db_pb * (mu + 2 * e + e * lambda^2)
                        g.w <- sigma_w * as.vector((0.5 / sigma * (cc + Partw - 1)) * Sigma_w)
                        
                        g.eta <- residual * valor
                        
                        g <- cbind(g.b, g.mu, g.u, g.w, g.eta)
                        
                        return(g)
                }
                G22 <- G2(b)
                
                g21 <- function(z){
                        zz <- (mu / lambda - erro * lambda) / sqrt(sigma)
                        dz <- pmax(dnorm(zz), 9.88131291682493e-324)
                        pz <- pmax(pnorm(zz), 9.88131291682493e-324)
                        fdp <- dz; cdf <- pz
                        valor <- as.vector((erro + mu) / sigma + fdp / cdf * lambda / sqrt(sigma))
                        eta <- z[c(k4 + 1):kk]
                        
                        g0 <- lapply(1:p, function(x) cbind(- Z * valor * eta[x]))
                        g <- matrix(unlist(g0), nrow = dim(Z)[1L], ncol = dim(Z)[2L]*p)
                        return(g)
                }
                G21 <- g21(b)
                
                V1 <- vcov(ols); V2 <- ginv(est$hessian); C <- t(G22) %*% G21; R <- t(G22) %*% G11
                V2_adj <- V2 + V2 %*% ((C %*% V1 %*% t(C)) - (R %*% V1 %*% t(C)) - (C %*% V1 %*% t(R))) %*% V2
                
                ep <- sqrt(diag(V2_adj))
                estZ <- lapply(length(b), function(x) b/ep)
                p_valor <- lapply(length(b), function(x) (1 - pnorm(abs(b/ep)))*2)
                ics <- lapply(length(b), function(x) cbind(b - qnorm(.975) * ep, b + qnorm(.975) * ep))
                result <- data.frame("Coefficient" = cbind(est$par), "Std.Err" = ep, "z" = estZ[[1]], 
                                     "P-value" = p_valor[[1]], "Lower limit" = ics[[1]][,1], "Upper limit" = ics[[1]][,2])
                
                mu.mod <- (- erro * s2u + mu * s2w)/ sigma
                s.mod <- sqrt((s2u * s2w)/sigma)
                uf <- mu.mod + s.mod * ((dnorm(-mu.mod / s.mod)) / (pnorm(mu.mod / s.mod)))
                ef <- ((pnorm(-s.mod + mu.mod / s.mod)) / (pnorm(mu.mod / s.mod))) * exp(- mu.mod + 0.5 * s.mod^2)
                
                lns2u <- est$par[c(k2 + 1):k3]
                names(lns2u) <- c('lns2u')
                cov.u <- ginv(est$hessian)[c(k2 + 1):k3,c(k2 + 1):k3]
                S2U <- deltaMethod(lns2u, "exp(lns2u)", vcov = cov.u)
                zS2U <- S2U$Estimate/S2U$SE
                pS2U <- (1 - pnorm(abs(zS2U)))*2
                icS2U <- lapply(length(S2U$Estimate), function(x) cbind(S2U$Estimate - qnorm(.975) * S2U$SE, S2U$Estimate + qnorm(.975) * S2U$SE))
                S2U0 <- data.frame(S2U$Estimate, S2U$SE, zS2U, pS2U, icS2U[[1]][,1], icS2U[[1]][,2])
                names(S2U0) <- names(result)
                rownames(S2U0) <- "s2u"
                
                lns2w <- est$par[c(k3 + 1):k4]
                names(lns2w) <- c('lns2w')
                cov.w <- ginv(est$hessian)[c(k3 + 1):k4,c(k3 + 1):k4]
                S2W <- deltaMethod(lns2w, "exp(lns2w)", vcov = cov.w)
                zS2W <- S2W$Estimate/S2W$SE
                pS2W <- (1 - pnorm(abs(zS2W)))*2
                icS2W <- lapply(length(S2W$Estimate), function(x) cbind(S2W$Estimate - qnorm(.975) * S2W$SE, S2W$Estimate + qnorm(.975) * S2W$SE))
                S2W0 <- data.frame(S2W$Estimate, S2W$SE, zS2W, pS2W, icS2W[[1]][,1], icS2W[[1]][,2])
                names(S2W0) <- names(result)
                rownames(S2W0) <- "s2w"
                
                s20 <- c(S2U$Estimate, S2W$Estimate)
                names(s20) <- c('s2u','s2w')
                cov.s <- ginv(est$hessian)[c(k2 + 1):k4,c(k2 + 1):k4]
                S2 <- deltaMethod(s20, "s2u + s2w", vcov = cov.s)
                zS2 <- S2$Estimate/S2$SE
                pS2 <- (1 - pnorm(abs(zS2)))*2
                icS2 <- lapply(length(S2$Estimate), function(x) cbind(S2$Estimate - qnorm(.975) * S2$SE, S2$Estimate + qnorm(.975) * S2$SE))
                S20 <- data.frame(S2$Estimate, S2$SE, zS2, pS2, icS2[[1]][,1], icS2[[1]][,2])
                names(S20) <- names(result)
                rownames(S20) <- "s2"
                
                Lambda <- deltaMethod(s20, "sqrt(s2u/s2w)", vcov = cov.s)
                zLam <- Lambda$Estimate/Lambda$SE
                pLam <- (1 - pnorm(abs(zLam)))*2
                icLam <- lapply(length(Lambda$Estimate), function(x) cbind(Lambda$Estimate - qnorm(.975) * Lambda$SE, Lambda$Estimate + qnorm(.975) * Lambda$SE))
                LAM <- data.frame(Lambda$Estimate, Lambda$SE, zLam, pLam, icLam[[1]][,1], icLam[[1]][,2])
                names(LAM) <- names(result)
                rownames(LAM) <- "lambda"
                
                RESP <- data.frame(rbind(result, S2U0, S2W0, S20, LAM))
                
                # LR test for endogeneity
                lnL1 <- -ll(est$par)
                lnL00 <- function(z) {
                        x_p <- z[1:k1]
                        mu_p <- z[c(k1 + 1):k2]
                        u_p <- z[c(k2 + 1):k3]
                        w_p <- z[c(k3 + 1):k4]
                        res_p <- rep(0,p)
                        
                        Y <- X_fr %*% x_p
                        mu <- mi %*% mu_p
                        U <- sigma_u %*% u_p
                        W <- sigma_w %*% w_p
                        Residual <- residual %*% res_p
                        
                        Sigma_u <- exp(U)
                        Sigma_w <- exp(W)
                        lambda <- sqrt(Sigma_u/Sigma_w)
                        sigma <- Sigma_u + Sigma_w
                        
                        e <- Y_fr - Y - Residual
                        
                        aa <- mu / sqrt(sigma) * sqrt(1 + lambda^-2)
                        pa <- pmax(pnorm(aa), 9.88131291682493e-324)
                        bb <- (mu / lambda - e * lambda) / sqrt(sigma)
                        pb <- pmax(pnorm(bb), 9.88131291682493e-324)
                        cc <- (e + mu)^2 / sigma
                        
                        ly.x.trun <- -0.5 * log(2 * pi) - 0.5 * log(sigma) - log(pa) + log(pb) - 0.5 * cc
                        
                        val <- sum(ly.x.trun)
                        
                        return(-val)
                }
                lnL0 <- -lnL00(est$par)
                LRtest <- 2 * (lnL1 - lnL0)
                LRp <- pchisq(LRtest, df = p, lower.tail = FALSE)
                
                # Wald test for endogeneity
                etas <- b[c(k4 + 1):kk]
                cov_etas <- V2_adj[c(k4 + 1):kk, c(k4 + 1):kk]
                waldtest <- t(etas) %*% ginv(cov_etas) %*% etas
                waldp <- pchisq(waldtest, df = p, lower.tail = FALSE)
                
                bias <- mean(Y_est) - mean(Y_fr)
                RMSE <- sqrt(var(Y_est) + bias^2)
                pearson <- cor(Y_est, Y_fr, method = "pearson")
                spearman <- cor(Y_est, Y_fr, method = "spearman")
                
                n <- length(Y_fr)
                W <- crossprod(residual) / n
                ## Likelihood of x
                lx <- n * 0.5 * (- p * log(2 * pi) - log(det(W))) - 0.5 * sum(mahalanobis(residual, 0, W))
                loglik <- lx - est$value
                lnL <- loglik
                K <- length(b)
                AIC <- - 2 * lnL + 2 * K
                BIC <- - 2 * lnL + log(n) * K
        }
        lista <- list("efficiency" = ef, "error" = erro, "fitted.y_without.correction" = Y_est,
                      "fitted.y_with.correction" = yest, "reg_IV" = result_IV, "table" = RESP,
                      "summary.ef" = summary(ef), "sd.ef" = sd(ef, na.rm = TRUE),
                      "value" = loglik, "AIC" = AIC, "BIC" = BIC, 
                      "cor pearson IV" = cor_p.iv, "cor spearman IV" = cor_s.iv,
                      "cor.pearson" = c(pearson), "cor.spearman" = c(spearman),
                      "LR chisq test" = matrix(cbind(LRtest, LRp), 1, 2, dimnames = list(c(),c("statistic", "P-value"))),
                      "Wald chisq test" = matrix(cbind(waldtest, waldp), 1, 2, dimnames = list(c(),c("statistic", "P-value"))),
                      "bias" = c(bias), "RMSE" = c(RMSE), 
                      "sample size" = n, "estimated parameters" = K)
        return(lista)
}



SF.prod <- function(model, fr.form, s2u.form, s2w.form, end.form, mu.form, data)
{
        required.packages <- c("Matrix","MASS","car")
        new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
        if(length(new.packages)) install.packages(new.packages)
        library(Matrix); library(MASS); library(car)
        if(is.character(model)) {
                distname <- tolower(model)
                model <- switch(distname, half = "half", exp = "exp", trunc = "trunc", 
                                half1s = "half1s", half2s = "half2s", 
                                exp1s = "exp1s", exp2s = "exp2s", 
                                trunc1s = "trunc1s", trunc2s = "trunc2s", NULL)
                if(is.null(model))
                        stop("unsupported distribution")
                if(distname == "half") {
                        return(SF.half(fr.form = fr.form, s2u.form = s2u.form, s2w.form = s2w.form, data = data))
                }
                if(distname == "exp") {
                        return(SF.exp(fr.form = fr.form, s2u.form = s2u.form, s2w.form = s2w.form, data = data))
                }
                if(distname == "trunc") {
                        return(SF.trunc(fr.form = fr.form, mu.form = mu.form, s2u.form = ~1, s2w.form = ~1, data = data))
                }
                if(distname == "half1s") {
                        return(SF.half1S(fr.form = fr.form, end.form = end.form, s2u.form = s2u.form, s2w.form = s2w.form, data = data))
                }
                if(distname == "half2s") {
                        return(SF.half2S(fr.form = fr.form, end.form = end.form, s2u.form = s2u.form, s2w.form = s2w.form, data = data))
                }
                if(distname == "exp1s") {
                        return(SF.exp1S(fr.form = fr.form, end.form = end.form, s2u.form = s2u.form, s2w.form = s2w.form, data = data))
                }
                if(distname == "exp2s") {
                        return(SF.exp2S(fr.form = fr.form, end.form = end.form, s2u.form = s2u.form, s2w.form = s2w.form, data = data))
                }
                if(distname == "trunc1s") {
                        return(SF.trunc1S(fr.form = fr.form, end.form = end.form, mu.form = mu.form, s2u.form = ~1, s2w.form = ~1, data = data))
                }
                if(distname == "trunc2s") {
                        return(SF.trunc2S(fr.form = fr.form, end.form = end.form, mu.form = mu.form, s2u.form = ~1, s2w.form = ~1, data = data))
                }
        }
        return(model)
}
