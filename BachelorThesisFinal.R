library("tidyverse")
library("reshape2")
library("ggpubr")
library("glmnet")
library("pls")


###########################################
###### Simulation of the latent Model ##### 
###########################################

latent.model <- function(n = 500,
                         m = 5,
                         p = 150,
                         sigma_X = 0.1,
                         sigma_Y = 0.1,
                         sigma_P = 1/sqrt(m),
                         alpha.true = NULL,
                         P = NULL,
                         seed = 123) {
  
  set.seed(seed)
  
  N <- matrix(rnorm(n*m, 0, 1), nrow = n, ncol = m)
  if (is.null(alpha.true)) {
    alpha.true <- c(-2*m, m:1)
  }
  
  f <- rnorm(n, 0, 1)
  Y <- N %*% alpha.true[-1] + alpha.true[1] + (sigma_Y * f) 
  
  if (is.null(P)) P <- sigma_P * matrix(rnorm(p*m, 0, 1), nrow = p, ncol = m)
  
  E <- matrix(rnorm(n*p, 0, 1), nrow = n, ncol = p)
  X <- N %*% t(P) + (sigma_X * E)
  
  beta.true <- c(alpha.true[1], P %*% solve(t(P) %*% P) %*% alpha.true[-1])
  
  return(list(X, Y, beta.true))
  
}



##################################
#######      SCRIPTS       ####### 
##################################

#### PCR as function ####


pcr <- function(m = 30, ...) {
  
  # Scaling
  
  mu.X <- apply(X.train, 2, mean)
  sd.X <- apply(X.train, 2, sd)
  
  mu.Y <- apply(Y.train, 2, mean)
  sd.Y <- apply(Y.train, 2, sd)
  
  X.norm.train <- scale(X.train)
  Y.norm.train <- scale(Y.train)
  
  # PCR
  
  m.pca.seq <- 2:m
  pca.results.MSE.y <- Inf
  pca.results.MSE.beta <- Inf
  m.pca.best <- NA
  beta.pca.best <- NA
  cor.pca <- rep(NA, m)
  
  # Global PCA on normalized training data
  
  fit.pca.train <- prcomp(X.norm.train, center=FALSE, scale.=FALSE, rank.=max(m.pca.seq)) 
  
  # PCR model selection
  
  for (m.pca in m.pca.seq) {
    
    # Reduced data set
    X.norm.train.pca <- fit.pca.train$x[, 1:m.pca]                                # prcomp gives you projected data already
    R.pca <- fit.pca.train$rotation[, 1:m.pca]                                    # with prcomp we only use projection matrix for coeffs
    
    data.pca.comp <- cbind(Y.norm.train, X.norm.train.pca)
    data.pca.comp <- as.data.frame(data.pca.comp)
    colnames(data.pca.comp) <- c("Y", paste0("PC", 1:m.pca))
    
    # PCR on reduced data set
    fit.pcr <- lm(Y ~ .,  data=data.pca.comp)                                     # we fit lm to get all coeffs
    beta.norm.pca.proj <- fit.pcr$coefficients                                    # we extract lm coeffs for norm/proj data
    
    # Coeffs for non-reduced normalized data set
    beta.norm.pca <- c(beta.norm.pca.proj[1], R.pca %*% beta.norm.pca.proj[-1])
    
    # Coeffs for original data set
    beta.pca <- beta.norm.pca*0
    
    beta.pca[1] <- sd.Y[1]*beta.norm.pca[1] + mu.Y[1] - sd.Y[1]*(mu.X/sd.X)%*%beta.norm.pca[-1]
    beta.pca[-1] <- sd.Y[1]*beta.norm.pca[-1]/sd.X
    
    # Predict original dataset
    Y.hat.pca <- X.test %*% beta.pca[-1] + beta.pca[1]
    
    # Score prediction
    MSE.y.candidate.pca <- mean((Y.test - Y.hat.pca)^2)
    MSE.beta.candidate.pca <- mean((beta.true - beta.pca)^2)
    cor.pca[m.pca] <- cor(Y.test, Y.hat.pca)
    
    if (MSE.y.candidate.pca < pca.results.MSE.y) {
      
      pca.results.MSE.y <- MSE.y.candidate.pca
      pca.results.MSE.beta <- MSE.beta.candidate.pca
      m.pca.best <- m.pca
      beta.pca.best <- beta.pca
      CondNumbPCA <- kappa((t(X.train) %*% X.train) %*% R.pca)
        
      
    } else {
      next;
    }
    
  }
  
  return(list(pca.results.MSE.y, pca.results.MSE.beta, m.pca.best, cor.pca, CondNumbPCA))
  
}


#### PLS as function ####

pls <- function(m = 30, ...) {
  
  
  # Scaling
  
  mu.X <- apply(X.train, 2, mean)
  sd.X <- apply(X.train, 2, sd)
  
  mu.Y <- apply(Y.train, 2, mean)
  sd.Y <- apply(Y.train, 2, sd)
  
  X.norm.train <- scale(X.train)
  Y.norm.train <- scale(Y.train)
  
  # PLS
  
  m.pls.seq <- 2:m
  pls.results.MSE.y <- Inf
  pls.results.MSE.beta <- Inf
  m.pls.best <- NA
  beta.pls.best <- NA
  cor.pls <- rep(NA, m)
  
  # Global PLS on normalized training data
  
  fit.pls.train <- plsr(Y.norm.train~X.norm.train, center=FALSE, scale=FALSE, ncomp=max(m.pls.seq))
  
  # PLS model selection
  
  for (m.pls in m.pls.seq) {
    
    # Reduced data set
    R.pls <- fit.pls.train$projection[, 1:m.pls]                                  # plsr requires to compute projection matrix first
    X.norm.train.pls <- X.norm.train%*%R.pls                                      # we project data ourselves to be sure
    
    data.pls.comp <- cbind(Y.norm.train, X.norm.train.pls)
    data.pls.comp <- as.data.frame(data.pls.comp)
    colnames(data.pls.comp) <- c("Y", paste0("PLS", 1:m.pls))
    
    # PLS on reduced data set
    fit.plsr <- lm(Y ~ .,  data=data.pls.comp)                                    # we fit lm to get all coeffs
    beta.norm.pls.proj <- fit.plsr$coefficients                                   # we extract lm coeffs for norm/proj data
    
    # Coeffs for non-reduced normalized data set
    beta.norm.pls <- c(beta.norm.pls.proj[1], R.pls%*%beta.norm.pls.proj[-1])
    
    # Coeffs for original data set
    beta.pls <- beta.norm.pls*0
    
    beta.pls[1] <- sd.Y[1]*beta.norm.pls[1] + mu.Y[1] - sd.Y[1]*(mu.X/sd.X)%*%beta.norm.pls[-1]
    beta.pls[-1] <- sd.Y[1]*beta.norm.pls[-1]/sd.X
    
    # Predict original dataset
    Y.hat.pls <- X.test%*%beta.pls[-1] + beta.pls[1]
    
    # Score prediction
    MSE.y.candidate.pls <- mean((Y.test - Y.hat.pls)^2)
    MSE.beta.candidate.pls <- mean((beta.true - beta.pls)^2)
    cor.pls[m.pls] <- cor(Y.test, Y.hat.pls)
    
    if (MSE.y.candidate.pls < pls.results.MSE.y) {
      
      pls.results.MSE.y <- MSE.y.candidate.pls
      pls.results.MSE.beta <- MSE.beta.candidate.pls
      m.pls.best <- m.pls
      beta.pls.best <- beta.pls
      CondNumbPLS <- kappa((t(X.train) %*% X.train) %*% R.pls)
      
    } else {
      next;
    }
    
  }
  
  return(list(pls.results.MSE.y, pls.results.MSE.beta, m.pls.best, cor.pls, CondNumbPLS))
  
}




#### RIDGE as function and without CV ####

ridge <- function(lambda.seq = seq(1, 0.01, -0.01), ...) {
  
  # Scaling
  
  mu.X <- apply(X.train, 2, mean)
  sd.X <- apply(X.train, 2, sd)
  
  mu.Y <- apply(Y.train, 2, mean)
  sd.Y <- apply(Y.train, 2, sd)
  
  X.norm.train <- scale(X.train)
  Y.norm.train <- scale(Y.train)
  
  # RIDGE
  
  ridge.results.MSE.y <- Inf
  ridge.results.MSE.beta <- Inf
  ridge.lambdaopt <- NA
  beta.ridge.best <- NA
  cor.ridge <- rep(NA, length(lambda.seq))
  
  # Global RIDGE on normalized training data
  
  fit.ridge <- glmnet(X.norm.train, Y.norm.train, alpha=0, lambda=lambda.seq)
  
  # RIDGE model selection
  
  for (lambda.ridge in lambda.seq) {
    
    beta.ridge <- coef(fit.ridge, s=lambda.ridge)
    Y.hat.ridge <- X.test %*% beta.ridge[-1] + beta.ridge[1]
    
    MSE.y.candidate.ridge <- mean((Y.test - Y.hat.ridge)^2)
    MSE.beta.candidate.ridge <- mean((beta.true - beta.ridge)^2)
    
    cor.ridge[which(lambda.ridge == lambda.seq)] <- cor(Y.test, Y.hat.ridge)
    
    if (MSE.y.candidate.ridge < ridge.results.MSE.y) {
      
      ridge.results.MSE.y <- MSE.y.candidate.ridge
      ridge.results.MSE.beta <- MSE.beta.candidate.ridge
      ridge.lambdaopt <- lambda.ridge
      beta.ridge.best <- beta.ridge
      CondNumbRidge <- kappa(t(X.train) %*% X.train + (ridge.lambdaopt * diag(ncol(X.train))))
      
    } else {
      next;
    }
    
  }
  
  cor.ridge[which(is.na(cor.ridge))] <- 0
  return(list(ridge.results.MSE.y, ridge.results.MSE.beta, ridge.lambdaopt, cor.ridge, CondNumbRidge))
  
}


#### LASSO as a function and without CV ####


lasso <- function(lambda.seq = seq(1, 0.01, -0.01), ...) {
  
  # Scaling
  
  mu.X <- apply(X.train, 2, mean)
  sd.X <- apply(X.train, 2, sd)
  
  mu.Y <- apply(Y.train, 2, mean)
  sd.Y <- apply(Y.train, 2, sd)
  
  X.norm.train <- scale(X.train)
  Y.norm.train <- scale(Y.train)
  
  # LASSO
  
  lasso.results.MSE.y <- Inf
  lasso.results.MSE.beta <- Inf
  lasso.lambdaopt <- NA
  beta.lasso.best <- NA
  cor.lasso <- rep(NA, length(lambda.seq))
  
  # Global LASSO on normalized training data
  
  fit.lasso <- glmnet(X.norm.train, Y.norm.train, alpha=1, lambda=lambda.seq)
  
  # LASSO model selection
  
  for (lambda.lasso in lambda.seq) {
    
    beta.lasso <- coef(fit.lasso, s=lambda.lasso)
    Y.hat.lasso <- X.test %*% beta.lasso[-1] + beta.lasso[1]
    
    MSE.y.candidate.lasso <- mean((Y.test - Y.hat.lasso)^2)
    MSE.beta.candidate.lasso <- mean((beta.true - beta.lasso)^2)
    
    cor.lasso[which(lambda.lasso == lambda.seq)] <- cor(Y.test, Y.hat.lasso)
    
    if (MSE.y.candidate.lasso < lasso.results.MSE.y) {
      
      lasso.results.MSE.y <- MSE.y.candidate.lasso
      lasso.results.MSE.beta <- MSE.beta.candidate.lasso
      lasso.lambdaopt <- lambda.lasso
      beta.lasso.best <- beta.lasso
      
    } else {
      next;
    }
    
  }
  
  cor.lasso[which(is.na(cor.lasso))] <- 0
  return(list(lasso.results.MSE.y, lasso.results.MSE.beta, lasso.lambdaopt, cor.lasso))
  
}



####### Simulation study ########

k <- 100

m.pcr.max <- 40
m.pls.max <- 40

lambda.ridge.seq <- seq(120, 0.01, length.out = 40)
lambda.lasso.seq <- seq(120, 0.01, length.out = 40)

# Creating result object
  
results.MSE <- matrix(rep(Inf, k * 8), ncol = 8)
results.MSE <- as.data.frame(results.MSE)
colnames(results.MSE) <- c("MSE Y PCR", "MSE beta PCR", 
                             "MSE Y PLS", "MSE beta PLS",
                             "MSE Y Ridge", "MSE beta Ridge",
                             "MSE Y Lasso", "MSE beta Lasso")
  
Parameters.opt <- matrix(rep(NA, k * 4), ncol = 4)
colnames(Parameters.opt) <- c("Princ. Comps.", "Dim. PLS", "Ridge Lambda", "Lasso Lambda")
CondNumb <- matrix(rep(Inf, 4 * k), ncol = 4)
colnames(CondNumb) <- c("Non-transf.", "PCA", "PLS", "Ridge")
CorrPCR <- matrix(rep(Inf, m.pcr.max * k), nrow = k)
CorrPLS <- matrix(rep(Inf, m.pls.max * k), nrow = k)
CorrRidge <- matrix(rep(Inf, length(lambda.ridge.seq) * k), nrow = k)
CorrLasso <- matrix(rep(Inf, length(lambda.lasso.seq) * k), nrow = k)
runtimes <- matrix(rep(NA, k * 4), ncol = 4)
colnames(runtimes) <- c("PCR", "PLS", "Ridge", "Lasso")

### Model parameters ###

n <- 500
m <- 5
p <- 150
sigma_X <- 0.1
sigma_Y <- 0.1
sigma_P <- 1/sqrt(m)

set.seed(0)
alpha.true <- c(-2*m, m:1)
P.true <- sigma_P * matrix(rnorm(p*m, 0, 1), nrow = p, ncol = m)
  
for (i in 1:k) {
  
  print(paste("Experiment",i))
  
  ### Model Creation ###
  
  LatMod <- latent.model(n=n, m=m, p=p, 
                         sigma_X=sigma_X, sigma_Y=sigma_Y, sigma_P=sigma_P, 
                         alpha.true=alpha.true, P=P.true,
                         seed=i)
  X <- LatMod[[1]]
  Y <- LatMod[[2]]
  beta.true <- LatMod[[3]]
  
  ### Dividing Train/Test Data ###
  
  train <- sample(n, n * 0.7, replace = FALSE)
  
  X.train <- as.matrix(X[train, ])
  Y.train <- as.matrix(Y[train])
  
  X.test <- as.matrix(X[-train, ])
  Y.test <- as.matrix(Y[-train])
  
  ### Calculate Condition Number ###
  
  CondNumb[i, 1] <- kappa(t(X.train) %*% X.train)
  
  ### PCA / PCR ###
  
  startpcr <- Sys.time()
  pca.iter <- pcr(m = m.pcr.max, X.train, Y.train, X.test, Y.test, beta.true)
  results.MSE[i, 1] <- pca.iter[[1]]
  results.MSE[i, 2] <- pca.iter[[2]]
  Parameters.opt[i, 1] <- pca.iter[[3]]
  CorrPCR[i, ] <- pca.iter[[4]]
  CondNumb[i, 2] <- pca.iter[[5]]
  endpcr <- Sys.time()
  runtimes[i, 1] <- difftime(endpcr, startpcr, units = "secs")
  
  ### PLS ###
  
  startpls <- Sys.time()
  pls.iter <- pls(m = m.pls.max, X.train, Y.train, X.test, Y.test, beta.true)
  results.MSE[i, 3] <- pls.iter[[1]]
  results.MSE[i, 4] <- pls.iter[[2]]
  Parameters.opt[i, 2] <- pls.iter[[3]]
  CorrPLS[i, ] <- pls.iter[[4]]
  CondNumb[i, 3] <- pls.iter[[5]]
  endpls <- Sys.time()
  runtimes[i, 2] <- difftime(endpls, startpls, units = "secs")
  
  ### RIDGE ###
  
  startridge <- Sys.time()
  ridge.iter <- ridge(lambda = lambda.ridge.seq, X.train, Y.train, X.test, Y.test, beta.true)
  results.MSE[i, 5] <- ridge.iter[[1]]
  results.MSE[i, 6] <- ridge.iter[[2]]
  Parameters.opt[i, 3] <- ridge.iter[[3]]
  CorrRidge[i, ] <- ridge.iter[[4]]
  CondNumb[i, 4] <- ridge.iter[[5]]
  endridge <- Sys.time()
  runtimes[i, 3] <- difftime(endridge, startridge, units = "secs")
  
  ### LASSO ###
  
  startlasso <- Sys.time()
  lasso.iter <- lasso(lambda = lambda.lasso.seq, X.train, Y.train, X.test, Y.test, beta.true)
  results.MSE[i, 7] <- lasso.iter[[1]]
  results.MSE[i, 8] <- lasso.iter[[2]]
  Parameters.opt[i, 4] <- lasso.iter[[3]]
  CorrLasso[i, ] <- lasso.iter[[4]]
  endlasso <- Sys.time()
  runtimes[i, 4] <- difftime(endlasso, startlasso, units = "secs")
    
}

colMeans(runtimes)


### Plots ###

#### MSE PLOTS ####

# PCR = 506, PLS = 44, RIDGE = 49, LASSO = 468

col <- colors()[c(506, 44, 49, 468)]

# MSE Y PLS vs PCR

MSEYdimred.sigma10 <- ggplot(melt(results.MSE[c(1, 3)])) + geom_boxplot(aes(y = value, x = variable, fill = variable)) +
  scale_fill_manual(name = "Methods", values = col[1:2]) +
  theme(axis.title.x=element_blank()) +
  theme(legend.title=element_text(size = 15),
        legend.text=element_text(size=15)) +
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=14), 
        axis.text.x=element_text(size = 10))


# MSE Y RIDGE vs LASSO

MSEYpen.sigma10 <- ggplot(melt(results.MSE[c(5, 7)])) + geom_boxplot(aes(y = value, x = variable, fill = variable)) +
  scale_fill_manual(name = "Methods", values = col[3:4]) +
  xlab(c("Ridge", "Lasso")) +
  theme(axis.title.x=element_blank()) +
  theme(legend.title=element_text(size = 15),
         legend.text=element_text(size=15)) +
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=14), 
        axis.text.x=element_text(size = 10))


MSEY0.1 <- ggarrange(MSEYdimred.sigma0.1, MSEYpen.sigma0.1)
MSEY0.1 <- annotate_figure(MSEY0.1, top = text_grob(paste("MSE Y, ", "Sigma_X = ", sigma_X), 
                                        color = "red", face = "bold", size = 14))

MSEY1 <- ggarrange(MSEYdimred.sigma1, MSEYpen.sigma1)
MSEY1 <- annotate_figure(MSEY1, top = text_grob(paste("MSE Y, ", "Sigma_X = ", sigma_X), 
                                        color = "red", face = "bold", size = 14))
MSEY10 <- ggarrange(MSEYdimred.sigma10, MSEYpen.sigma10)
MSEY10 <- annotate_figure(MSEY10, top = text_grob(paste("MSE Y, ", "Sigma_X = ", sigma_X), 
                                      color = "red", face = "bold", size = 14))




# MSE beta PLS vs PCR

MSEbdimred.sigma10 <- ggplot(melt(results.MSE[c(2, 4)])) + geom_boxplot(aes(y = value, x = variable, fill = variable)) +
  scale_fill_manual(name = "Methods", values = col[1:2]) +
  theme(axis.title.x=element_blank()) +
  theme(legend.title=element_text(size = 15),
        legend.text=element_text(size=15)) +
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=14), 
        axis.text.x=element_text(size = 7))

# MSE beta RIDGE vs LASSO

MSEbpen.sigma10 <- ggplot(melt(results.MSE[c(6, 8)])) + geom_boxplot(aes(y = value, x = variable, fill = variable)) +
  scale_fill_manual(name = "Methods", values = col[3:4]) +
  theme(axis.title.x=element_blank()) +
  theme(legend.title=element_text(size = 15),
        legend.text=element_text(size=15)) +
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=14), 
        axis.text.x=element_text(size = 7))


MSEb0.1 <- ggarrange(MSEbdimred.sigma0.1, MSEbpen.sigma0.1)
MSEb0.1 <- annotate_figure(MSEb0.1, top = text_grob(paste("MSE beta, ", "Sigma_X = ", sigma_X), 
                                                    color = "red", face = "bold", size = 14))
MSEb1 <- ggarrange(MSEbdimred.sigma1, MSEbpen.sigma1)
MSEb1 <- annotate_figure(MSEb1, top = text_grob(paste("MSE beta, ", "Sigma_X = ", sigma_X), 
                                                color = "red", face = "bold", size = 14))
MSEb10 <- ggarrange(MSEbdimred.sigma10, MSEbpen.sigma10)
MSEb10 <- annotate_figure(MSEb10, top = text_grob(paste("MSE beta, ", "Sigma_X = ", sigma_X), 
                                                 color = "red", face = "bold", size = 14))




#### Parameteres plot ####


Parameters.opt.plt <- melt(Parameters.opt)  
Parameters.opt.plt[(k * 2 + 1):(k * 4), 3] <- p

# Parameters PLS vs PCR vs RIDGE vs LASSO

param.sigma0.1 <- ggplot(Parameters.opt.plt[1:(nrow(Parameters.opt.plt) / 2), ]) + 
 geom_boxplot(aes(y = value, x = Var2, fill = Var2)) + 
 geom_hline(yintercept = m, size = 1.5) +
 scale_fill_manual(name = "Methods", values = col,
                   labels = c("Principal Components", "Dimension PLS")) +
  theme(axis.title.x=element_blank()) +
  theme(legend.title=element_text(size = 15),
        legend.text=element_text(size=15)) +
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=14), 
        axis.text.x=element_text(size = 10)) +
  ggtitle(label = "Sigma_X = 0.1")

param.sigma1 <- ggplot(Parameters.opt.plt[1:(nrow(Parameters.opt.plt) / 2), ]) + 
  geom_boxplot(aes(y = value, x = Var2, fill = Var2)) + 
  geom_hline(yintercept = m, size = 1.5) +
  scale_fill_manual(name = "Methods", values = col,
                    labels = c("Principal Components", "Dimension PLS")) +
  theme(axis.title.x=element_blank()) +
  theme(legend.title=element_text(size = 15),
        legend.text=element_text(size=15)) +
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=14), 
        axis.text.x=element_text(size = 10)) +
  ggtitle(label = "Sigma_X = 1")

param.sigma10 <- ggplot(Parameters.opt.plt[1:(nrow(Parameters.opt.plt) / 2), ]) + 
  geom_boxplot(aes(y = value, x = Var2, fill = Var2)) + 
  geom_hline(yintercept = m, size = 1.5) +
  scale_fill_manual(name = "Methods", values = col,
                    labels = c("Principal Components", "Dimension PLS")) +
  theme(axis.title.x=element_blank()) +
  theme(legend.title=element_text(size = 15),
        legend.text=element_text(size=15)) +
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=14), 
        axis.text.x=element_text(size = 10)) +
  ggtitle(label = "Sigma_X = 10")

param <- ggarrange(param.sigma0.1, param.sigma1, param.sigma10, common.legend = TRUE, ncol = 3, legend = "bottom")
param <- annotate_figure(param, top = text_grob("Real Latent dimension (black) and dimensions for PCR/PLS", 
                                                  color = "red", face = "bold", size = 14))



# Mean over the columns for the correlation and line plots

CorrPCR0.1 <- colMeans(CorrPCR, na.rm = TRUE)
CorrPCR0.1[1] <- 0
CorrPLS0.1 <- colMeans(CorrPLS, na.rm = TRUE)
CorrPLS0.1[1] <- 0

CorrRidge0.1 <- colMeans(CorrRidge, na.rm = TRUE)
CorrLasso0.1 <- colMeans(CorrLasso, na.rm = TRUE)

CorrPCR1 <- colMeans(CorrPCR, na.rm = TRUE)
CorrPCR1[1] <- 0
CorrPLS1 <- colMeans(CorrPLS, na.rm = TRUE)
CorrPLS1[1] <- 0

CorrRidge1 <- colMeans(CorrRidge, na.rm = TRUE)
CorrLasso1 <- colMeans(CorrLasso, na.rm = TRUE)

CorrPCR10 <- colMeans(CorrPCR, na.rm = TRUE)
CorrPCR10[1] <- 0
CorrPLS10 <- colMeans(CorrPLS, na.rm = TRUE)
CorrPLS10[1] <- 0

CorrRidge10 <- colMeans(CorrRidge, na.rm = TRUE)
CorrLasso10 <- colMeans(CorrLasso, na.rm = TRUE)

### SIGMA = 0.1 ###


corr.sigma0.1.zoom <- ggplot() + geom_line(aes(x = 1:m.pcr.max, y = CorrPCR0.1,  color = col[2]), size = 1.5, alpha = 0.7) +
  geom_line(aes(x = 1:m.pls.max, y = CorrPLS0.1, color = col[3]), size = 1.5, alpha = 0.7) + 
  geom_line(aes(x = lambda.ridge.seq / 3, y = CorrRidge0.1, color = col[4]), size = 1.5, alpha = 0.7) + 
  geom_line(aes(x = lambda.lasso.seq / 3, y = CorrLasso0.1, color = col[1]), size = 1.5, alpha = 0.7) +
  scale_color_manual(name = "Methods", values = col, 
                     labels = c("PCR", "PLS", "Ridge", "Lasso"),
                     ) +
  scale_x_continuous(sec.axis = sec_axis(~. * 3, name = "Lambda"), name = "Dimensions") +
  ylab("Correlation") +
  coord_cartesian(ylim=c(0.99,1)) +
  ggtitle("sigma_X = 0.1, zoom")


### SIGMA = 1 ###

corr.sigma1.zoom <- ggplot() + geom_line(aes(x = 1:m.pcr.max, y = CorrPCR1,  color = col[2]), size = 1.5, alpha = 0.7) +
  geom_line(aes(x = 1:m.pls.max, y = CorrPLS1, color = col[3]), size = 1.5, alpha = 0.7) + 
  geom_line(aes(x = lambda.ridge.seq / 3, y = CorrRidge1, color = col[4]), size = 1.5, alpha = 0.7) + 
  geom_line(aes(x = lambda.lasso.seq / 3, y = CorrLasso1, color = col[1]), size = 1.5, alpha = 0.7) +
  scale_color_manual(name = "Methods", values = col, 
                     labels = c("PCR", "PLS", "Ridge", "Lasso"),
  ) +
  scale_x_continuous(sec.axis = sec_axis(~ . * 3, name = "Lambda"), name = "Dimensions") +
  ylab("Correlation") +
  coord_cartesian(ylim=c(0.95,0.99)) +
  ggtitle("sigma_X = 1, zoom")



### SIGMA = 10 ###

corr.sigma10.zoom <- ggplot() + geom_line(aes(x = 1:m.pcr.max, y = CorrPCR10,  color = col[2]), size = 1.5, alpha = 0.7) +
  geom_line(aes(x = 1:m.pls.max, y = CorrPLS10, color = col[3]), size = 1.5, alpha = 0.7) + 
  geom_line(aes(x = lambda.ridge.seq / 3, y = CorrRidge10, color = col[4]), size = 1.5, alpha = 0.7) + 
  geom_line(aes(x = lambda.lasso.seq / 3, y = CorrLasso10, color = col[1]), size = 1.5, alpha = 0.7) +
  scale_color_manual(name = "Methods", values = col, 
                     labels = c("PCR", "PLS", "Ridge", "Lasso"),
  ) +
  scale_x_continuous(sec.axis = sec_axis(~ . * 3, name = "Lambda"), name = "Dimensions") +
  ylab("Correlation") +
  coord_cartesian(ylim=c(0.2,0.35)) +
  ggtitle("sigma_X = 10, zoom")


rm(Corr.plt)
Corr.plt  <- ggarrange(corr.sigma0.1.zoom, corr.sigma1.zoom, corr.sigma10.zoom, common.legend = TRUE, legend = "bottom", ncol = 3)
Corr.plt <- annotate_figure(Corr.plt, top = text_grob(paste("Correlation for all methods"), 
                                                                     color = "red", face = "bold", size = 14))



# Condition Numbers

### SIGMA = 0.1 ###

condnumb0.1 <- ggplot(melt(CondNumb)) + 
  geom_boxplot(aes(x = Var2, y = value, fill = Var2)) +
  scale_fill_manual(name = "Condition Numbers", values = c("darkgrey", col[1:3])) +
  ggtitle(label = paste("Sigma_X = ", sigma_X)) + 
  theme(axis.title.x=element_blank(), #remove x axis labels
  )

### SIGMA = 1 ###

condnumb1 <- ggplot(melt(CondNumb)) + 
  geom_boxplot(aes(x = Var2, y = value, fill = Var2)) +
  scale_fill_manual(name = "Condition Numbers", values = c("darkgrey", col[1:3])) +
  ggtitle(label = paste("Sigma_X = ", sigma_X)) + 
  theme(axis.title.x=element_blank(), #remove x axis labels
  )

### SIGMA = 10 ###

condnumb10 <- ggplot(melt(CondNumb)) + 
  geom_boxplot(aes(x = Var2, y = value, fill = Var2)) +
  scale_fill_manual(name = "Condition Numbers", values = c("darkgrey", col[1:3])) +
  ggtitle(label = paste("Sigma_X = ", sigma_X)) + 
  theme(axis.title.x=element_blank(), #remove x axis labels
  )


condnum <- ggarrange(condnumb0.1, condnumb1, condnumb10, common.legend = TRUE, legend = "bottom", ncol = 3)
condnum <- annotate_figure(condnum, top = text_grob("Condition Numbers for X.train and after PCR/PLS/Ridge", 
                                                color = "red", face = "bold", size = 14))


