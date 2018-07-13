## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = FALSE)
# install.packages(c("fda","png"))
## Load packages
library("fda")
data(growth)

## ---- fig.height=4, fig.align='center'-----------------------------------
library("png")
library("grid")
FDA_Book_Pic <- png::readPNG("./images/IntroFDA_Book.png")
grid::grid.raster(FDA_Book_Pic)

## ----FD1, out.width = "100%"---------------------------------------------
knitr::include_graphics("./images/boygrowth.png")

## ----FD2, out.width = "100%"---------------------------------------------
knitr::include_graphics("./images/Fig_POI.png")

## ----FD4, out.width = "100%"---------------------------------------------
knitr::include_graphics("./images/Fig_Electr.png")

## ----FD5, out.width = "100%"---------------------------------------------
knitr::include_graphics("./images/Fig_PM10_1.png")

## ----FD6, out.width = "100%"---------------------------------------------
knitr::include_graphics("./images/Fig_PM10_2.png")

## ---- fig.width=8, fig.height=3, fig.align='center'----------------------
par(mfrow=c(1,2), mar=c(5.1, 4.1, 2.1, 2.1))
matplot(x=growth$age, y=growth$hgtf[,1:5], type="p", lty=1, pch=21, 
        xlab="age", ylab="height (cm)", cex.lab=1.2,col=1:5,bg=1:5,
        main=expression(x[i](t[j])))
matplot(x=growth$age, y=growth$hgtf[,1:5], type="l", lty=1, pch=1, 
        xlab="age", ylab="height (cm)", cex.lab=1.2,
        main=expression(X[i](t)))
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))

## ---- echo=TRUE, out.width = "80%", fig.align='center'-------------------
bspl_bf <- create.bspline.basis(rangeval=c(0,1), norder=3, 
                                     breaks=seq(0,1,len=7))
plot(bspl_bf, main="B-spline Basis Functions", xlab="[a,b]=[0,1]")

## ---- echo=TRUE, out.width = "80%", fig.align='center'-------------------
fourier_bf <- create.fourier.basis(rangeval=c(0,1), nbasis=5)
plot(fourier_bf, main="Fourier Basis Functions", xlab="[a,b]=[0,1]")

## ---- fig.align='center', out.width = "80%", echo=TRUE-------------------
matplot(x=growth$age, y=growth$hgtf[,1:5], type="p", lty=1, pch=21, 
        xlab="age", ylab="height (cm)", cex.lab=1.2,col=1:5,bg=1:5,
        main="5 Girls in Berkeley Growth Study")

## ---- fig.align='center', out.width = "80%", echo=TRUE-------------------
SmObj <- smooth.basisPar(growth$age,y=growth$hgtf[,1:5],lam=0.1)
plot(SmObj$fd, xlab="age", ylab="height (cm)", cex.lab=1.2,
     main="5 Girls in Berkeley Growth Study", lty=1)

## ---- fig.align='center', out.width = "80%", echo=TRUE-------------------
plot(deriv(SmObj$fd, 1), xlab="age", 
     ylab="growth rate (cm / year)", cex.lab=1.2,
     main="5 Girls in Berkeley Growth Study (1st derivative)", lty=1)

## ---- fig.align='center', out.width = "80%", echo=TRUE-------------------
plot(deriv(SmObj$fd, 2), xlab="age", cex.lab=1.2,
     ylab="growth acceleration (cm / year^2)",
     main="5 Girls in Berkeley Growth Study (2nd derivative)", lty=1)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
## BM.Mean <- rowMeans(BM.mat)
## BM.SD   <- apply(BM.mat, 1, sd)

## ---- echo=FALSE, fig.align='center', fig.width=8, fig.height=4----------
set.seed(109)
n       <- 50
J       <- 500
BM.mat  <- matrix(0, ncol=n, nrow=J)
for(i in 1:n){BM.mat[,i] <- cumsum(rnorm(J, sd = 1/100))}
BM.Mean <- rowMeans(BM.mat)
BM.SD   <- apply(BM.mat,1,sd)
xx          <- seq(0,1,len=J) 
par(mfrow=c(1,2))
matplot(x=xx, y=BM.mat, xlab="", ylab="", type="l", col=gray(.7), lty=1, main="Brownian Motions")
lines(x=xx, y=BM.Mean)
lines(x=xx, y=BM.SD, lty=2)
legend("topleft", lty = c(1,2), legend = c("Sample Mean","Sample SD"))
matplot(x=xx, y=BM.mat, xlab="", ylab="", type="l", col=gray(.7), lty=1, main="Brownian Motion")
lines(x=c(0,1), y=c(0,0), lty=1)
lines(x=c(0,1), y=c(0,sqrt(J*(1/100)^2)), lty=2)
legend("topleft", lty = c(1,2), legend = c("True Mean","True SD"))
par(mfrow=c(1,1))

## ---- eval=FALSE, echo=TRUE----------------------------------------------
## BM.Cov <- var(BM.mat)

## ---- echo=FALSE, fig.align='center', out.width="100%"-------------------
BM.cov <- var(t(BM.mat))
slct   <- c(seq.int(1,500,by=20),500)
par(mfrow=c(1,2), mar=c(1, 1.1, 1.2, 1.1))
persp(xx[slct], xx[slct], BM.cov[slct,slct], xlab="s", ylab="t", zlab="",main="Sample Covariance Function",
      theta = -40, phi = 20, expand = .75, col = "blue", shade = 1.05)
x <- seq(0, 1, length= 30); y <- x
f <- function(x, y){min(x,y)}
f <- Vectorize(f)
z <- outer(x, y, f)
persp(x, y, z, xlab="s", ylab="t", zlab="",
      main="True Covariance Function",
      theta = -40, phi = 20, expand = .75, col = "blue", shade = 1.05)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))

## ----echo=FALSE, fig.align='center', fig.height=3------------------------
BSPL.basis <- create.bspline.basis(rangeval=c(0,1), nbasis=200)
BM.fd      <- smooth.basis(argvals=xx, y=BM.mat, fdParobj=BSPL.basis)
BM.pca     <- pca.fd(BM.fd$fd, nharm=3)
par(mfrow=c(1,3), mar=c(2.1, 1.1, 4.1, 1.1))
persp(xx[slct], xx[slct], BM.cov[slct,slct], xlab="s", ylab="t", zlab="",
      main="Sample Covariance Function", theta = -40, phi = 20, expand = .75, col = "blue", shade = 1.05)
invisible(plot(BM.pca$values[1:3], type="o", ylab="", main="Estimated Eigenvalues (p=3)"))
invisible(plot(BM.pca$harmonics, lwd=3, ylab="", main="Estimated FPC's (p=3)"))
par(mfrow=c(1,1))

## ----echo=TRUE, echo=FALSE, fig.align='center', fig.height=4-------------
BM.pca <- pca.fd(BM.fd$fd, nharm=15)
v_hat_mat  <- eval.fd(xx, BM.pca$harmonics)
xi_hat_mat <- BM.pca$scores

# FPCA-approx
X_fpca_fit <- BM.Mean + (xi_hat_mat %*% t(v_hat_mat)) 

# plot
par(mfrow=c(1,1), mar=c(5.1, 4.1, 2.1, 2.1))
invisible(plot(BM.fd$fd[1], lwd=1.3, main="Approximation of Brownian Motion (p=15)"))
lines(y=X_fpca_fit[1,], x=xx, col="red", lwd=1.3)
legend("topright", lty = c(1,1), col=c("black","red"), legend = c("True","Approx"))

## ---- echo=FALSE, eval=TRUE----------------------------------------------
# Data
# BOA        <- read.csv(file = "https://raw.githubusercontent.com/lidom/Teaching_Repo/master/stock_prices.csv",
#                        header = TRUE, sep = " ", dec = ",")
BOA        <- read.csv(file = "stock_prices.csv", header = TRUE, sep = " ", dec = ",")
Dates      <- BOA$date
BOA        <- BOA[,-1]
BOA        <- data.matrix(BOA)
BOA        <- BOA[-which(Dates=="08/26/2004"),]# Outlier
n          <- dim(BOA)[1]
J          <- dim(BOA)[2]
Times      <- seq(0,6.5,length=J)

# Cumulative log-return functions (raw data)
log_BOA    <- log(BOA) - matrix(log(BOA)[,1], nrow=n, ncol=J)

# B-spline basis functions
bspl_basis <- create.bspline.basis(rangeval=c(0,6.5),nord=4,nb=200)

# Cumulative log-return functions (with basis functions)
log_BOA_fd <- Data2fd(Times,t(log_BOA),basisobj = bspl_basis)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
## # Plot functional data
## plot(log_BOA_fd, xlab="Trading Hours",ylab="",lwd=1, col=gray(.5),
##      main="Cumulative Log-Return Functions")
## lines(log_BOA_fd[1:10],lwd=1.5, lty=1)

## ---- eval=TRUE, echo=FALSE, fig.align='center', fig.height=4------------
par(mfrow=c(1,1), mar=c(5.1, 4.1, 2.1, 2.1))
invisible(plot(log_BOA_fd, xlab="Trading Hours",ylab="",lwd=1, col=gray(.5), 
               main="Cumulative Log-Return Functions"))
lines(log_BOA_fd[1:10],lwd=1.5, lty=1)

## ----echo=FALSE, fig.align='center', fig.cap="Plot of mean function for BOA cumulative returns. Point-wise 95% confidence intervals are included in red."----
muhat          <- mean.fd(log_BOA_fd)
sdhat          <- sd.fd(log_BOA_fd)
SE_hat_U       <- fd(basisobj=bspl_basis) # create upper CI bound
SE_hat_L       <- fd(basisobj=bspl_basis) # create lower CI bound
SE_hat_U$coefs <-  2*sdhat$coefs/sqrt(n) + muhat$coefs
SE_hat_L$coefs <- -2*sdhat$coefs/sqrt(n) + muhat$coefs 
# plot
invisible(plot.fd(SE_hat_U,ylim=c(-0.002,0.002),col='red',lty=2,
                  xlab="Trading Hours",ylab="Returns", main="Mean Function with 95% CI"))
invisible(plot.fd(SE_hat_L,add=TRUE,col='red',lty=2))
invisible(plot.fd(muhat,add=TRUE))

## ----echo=FALSE, fig.cap="EFPC's of BOA cumulative returns versus the theoretical eigenfunctions of Brownian Motions.", fig.align='center'----
log_BOA_fpca  <- pca.fd(log_BOA_fd, nharm=4)
# Eigenfunctions of BM process
efwp_fun      <- function(x, k=1){sqrt(2)*sin((k-1/2)*pi*x)}
# integrate(function(x){(efwp_fun(x,k=1))^2}, lower = 0, upper = 1)
EFWP_mat      <- cbind(efwp_fun(x=xx,k=1), efwp_fun(x=xx,k=2), efwp_fun(x=xx,k=3), efwp_fun(x=xx,k=4))
# plot
par(mfrow=c(1,2))
invisible(plot(log_BOA_fpca$harmonics, main="Estimated Eigenfunctions of\nCumulative Return Functions",
               xlab="Trading Hours", ylab=""))
matplot(x = xx, y=EFWP_mat, type="l",  main="True Eigenfunctions of\nBrownian Motion",ylab="",xlab="")

## ---- echo=TRUE----------------------------------------------------------
set.seed(9000)
n    <- 1000
grid <- seq(0, 1, length = 101)
##
beta <- sin(grid * 2 * pi)
X    <- matrix(0, nrow=n, ncol=length(grid))
for(i2 in 1:n){
 X[i2,] <- X[i2,]+rnorm(length(grid), 0, 1)
 X[i2,] <- X[i2,]+runif(1, 0, 5)
 X[i2,] <- X[i2,]+rnorm(1, 1, 0.2)*grid
 for(j2 in 1:10){
  e      <- rnorm(2, 0, 1/j2^(2))
  X[i2,] <- X[i2,]+e[1]*sin((2*pi)*grid*j2)
  X[i2,] <- X[i2,]+e[2]*cos((2*pi)*grid*j2)
}}
## Note the Integral-Approximation using '* .01':
Y = X %*% beta * .01 + rnorm(n, 0, .4)

## ---- echo=TRUE----------------------------------------------------------
matplot(x = grid, y=t(X), type = "l", lty=1)

## ---- echo=TRUE----------------------------------------------------------
library("refund")

## Roughness Penalty
fit.RP   <- pfr(Y ~ lf(X, bs = "ps", k = 50))

## FPCA
fit.FPCA <- pfr(Y ~ fpc(X, pve = 0.99))

## ---- echo=FALSE, fig.align="center"-------------------------------------
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
##
coefs <-  data.frame(grid      = grid, 
                     FPCA      = coef(fit.FPCA)$value, 
                     Penalized = coef(fit.RP)$value,
                     Truth     = beta)
coefs.m           <- melt(coefs, id = "grid")
colnames(coefs.m) <- c("grid", "Method", "Value")
##
# dev.new(width=6/2,height=3/2)
ggplot(coefs.m, aes(x=grid,y=Value,color=Method,group=Method), 
       width=12,height=6) + geom_path() + theme_bw()

## ---- echo=TRUE----------------------------------------------------------
##
beta <- -1 * dnorm(grid, mean=.20, sd=.03) + 
         3 * dnorm(grid, mean=.50, sd=.04) + 
         1 * dnorm(grid, mean=.75, sd=.05)

## Note the Integral-Approximation using '* .01':
Y = X %*% beta * .01 + rnorm(n, 0, .4)

## ---- echo=FALSE, fig.align="center"-------------------------------------
## Roughness Penalty
fit.RP   <- pfr(Y ~ lf(X, bs = "ps", k = 50))
## FPCA
fit.FPCA <- pfr(Y ~ fpc(X, pve = 0.99))
##
coefs <-  data.frame(grid      = grid, 
                     FPCA      = coef(fit.FPCA)$value, 
                     Penalized = coef(fit.RP)$value,
                     Truth     = beta)
coefs.m           <- melt(coefs, id = "grid")
colnames(coefs.m) <- c("grid", "Method", "Value")
##
ggplot(coefs.m, aes(x=grid,y=Value,color=Method,group=Method), 
       width=12,height=6) + geom_path() + theme_bw()

