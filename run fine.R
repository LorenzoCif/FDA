rm(list=ls())

library(npreg)
library(R.matlab)
library(rlist)
library(SCBmeanfd)
library(fda)
library(refund)
library(fields)
library(lpSolve)
library(tidyverse)

wd_base <- "C:/Users/gaiom/Desktop/FUNZIONALI/progetto/dati sport/SportDB/RUN"

# IMPORT DEI DATI ---------------------------------------------------------

# covariate ---------------------------------------------------------------

covariate <- matrix(NA, ncol = 7, nrow = 10)
for(i in 1:10){
  filename <- paste(paste(paste(wd_base, "/S", sep = ""), 
                          i, sep = ""), "/CRD1/Dem.txt", sep = "")
  dati <- read.table(filename, header = T)
  covariate[i,] <- (as.matrix(dati))
}

covariate <- as.data.frame(covariate)
names(covariate) <- c("Sex", "Age", "Weight", "Height", "Smoker", "Alcool",
                      "Training_Rate")
covariate

# hr e br -----------------------------------------------------------------

hr <- list()
br <- list()

for(i in 1:10){
  filename <- paste(paste(paste(wd_base, "/S", sep = ""), 
                          i, sep = ""), "/CRD1/Data.mat", sep = "")
  hr <- list.append(hr, readMat(filename)$Data[[1]])
  br <- list.append(br, readMat(filename)$Data[[3]])
}

# PREPROCESSING -----------------------------------------------------------

# imputazione di valori mal registrati ------------------------------------

# correggiamo gli errori nei dati

# hr
plot(hr[[10]], type = "l")
sum(hr[[10]]==0)
# imputiamo la media dei due valori corretti più vicini
hr[[10]][1537:1556] <- mean(hr[[10]][1536], hr[[10]][1557])
hr[[10]][1896:1918] <- mean(hr[[10]][1895], hr[[10]][1919])
plot(hr[[10]], type = "l")

# serie di lunghezze diverse
HR.temp <- matrix(NA, nrow = max(unlist(lapply(hr, length))), ncol = 10)
dim(HR.temp)
for(i in 1:10){
  HR.temp[1:length(hr[[i]]),i] <- hr[[i]]
}
matplot(HR.temp, type = "l",col = "grey50", xlab = "tempo", ylab = "hr", cex.axis = 1.3,
        cex.lab = 1.5)
lunghezze <- unlist(lapply(hr, length))
hr.fine <- rep(NA, 10)
for(i in 1:10){
  points(lunghezze[i], hr[[i]][lunghezze[i]], pch = 19)
  hr.fine[i] <- hr[[i]][lunghezze[i]]
}

# in br
br.clean <- function(br){
  idx <- which(br>100)
  if (length(idx)>0){
    for (i in 1:length(idx)){
      br[idx[i]] <- mean(br[idx[i]-1], br[idx[i]+1])
    }
  }
  br
}
br <- sapply(br, br.clean)

# da dominio temporale a domninio spaziale --------------------------------
# modifichiamo le serie hr e br facendo sì che il loro dominio sia lo spazio (uguale
# per tutti i soggetti) e non più il tempo (serie di lunghezze diverse per ogni
# soggetto).

# tempi in cui cambia la pendenza del percorso (dai file con le note, presenti
# nel dataset per ciadcun corridore)
int <- matrix(NA, 10, 4)
int[1,] <- c(7*60+45, 18*60+32, 25*60+47, 43*60+54)
int[2,] <- c(5*60+7, 13*60+54, 19*60+13, 32*60+50)
int[3,] <- c(8*60+5, 17*60+55, 23*60+48, 38*60+55)
int[4,] <- c(5*60+16, 15*60+3, 20*60+35, 34*60+31)
int[5,] <- c(7*60+36, 15*60+58, 21*60+54, 37*60+7)
int[6,] <- c(6*60+41, 14*60+43, 22*60+4, 34*60+45)
int[7,] <- c(05*60+24, 13*60+35, 17*60+15, 29*60+40)
int[8,] <- c(6*60+55, 14*60+31, 19*60+38, 30*60+11)
int[9,] <- c(7*60+19, 17*60+35, 26*60+45, 38*60+46)
int[10,] <- c(5*60+24, 16*60+22, 20*60+28, 32*60+5)

# lunghezza in km di ciascun tratto del percorso
flessi = c(1.3, 1.2, 1, 2.6)

# tempo di percorrenza di ciascun tratto a pendenza costante
dt <- matrix(NA, 10, 4)
dt[,1] <- int[,1]
dt[,2:4] <- int[,2:4]-int[,1:3]
dt

# velocità medie di ciascun tratto a pendenza costante 
# (assumendo velocita costante)
vel = matrix(NA, 10, 4)
vel = t(apply(dt, 1, function(x) flessi/x)) # velocita' = spazio/tempo
vel

# velocità in ogni istante temporale
velfun = list()
for (i in 1:10){
  velfun[[i]] = c(rep(vel[i,1], dt[i,1]), 
                  rep(vel[i,2], dt[i,2]),
                  rep(vel[i,3], dt[i,3]),
                  rep(vel[i,4], dt[i,4]))
}

# lisciamento delle velocità
# splines di lisciamento in ogni curva
smoothspeed =list()
for (i in 1:10){
  s = ss(1:int[i,4], velfun[[i]], df=35, m = 1)
  smoothspeed[[i]] = predict(s, 1:int[i,4])[,2]
}
par(mfrow=c(2,5))
for(i in 1:10){
  plot(smoothspeed[[i]], type = "l", xlab = "time",
       ylab = "speed", main = paste("soggetto ", i))
}
par(mfrow = c(1,1))

# spazio = somma delle distanze percorse al secondo
space = list()
for (i in 1:10){
  space[[i]] = cumsum(velfun[[i]]) 
}
lapply(space, range) # okay
# le serie che erano più lunghe nel tempo avranno osservazioni
# in griglie piu fitte nello spazio, ma sempre entro il range (0, 6.1)

# scegliamo 200 punti equispaziati dove osservare hr e br
space_ob <- seq(0, 6.1, length = 200)

# in corrispondenza di ciascuno dei punti della griglia, mettiamo il valore della serie
# rilveato se presente o il valore rilevato più vicino

# hr
HR <- matrix(NA, ncol = 10, nrow = length(space_ob))
for(j in 1:10){
  for(i in 1:length(space_ob)){
    HR[i, j] <- hr[[j]][which.min(abs(space[[j]]-space_ob[i]))]
  }
}
rownames(HR) = space_ob
matplot(HR, type = "l", lty = 1, xlab = "spazio", ylab = "HR",cex.axis = 1.5,
        cex.lab = 1.3)

# br
BR <- matrix(NA, ncol = 10, nrow = length(space_ob))
for(j in 1:10){
  for(i in 1:length(space_ob)){
    BR[i, j] <- br[[j]][which.min(abs(space[[j]]-space_ob[i]))]
  }
}
matplot(BR, type = "l", lty = 1, xlab = "spazio", ylab = "BR",cex.axis = 1.5,
        cex.lab = 1.3)

# speed
SPEED <- matrix(NA, ncol = 10, nrow = length(space_ob))
for(j in 1:10){
  for(i in 1:length(space_ob)){
    SPEED[i, j] <- smoothspeed[[j]][which.min(abs(space[[j]]-space_ob[i]))]
  }
}
matplot(SPEED, type = "l", xlab = "spazio", ylab = "speed",cex.axis = 1.3,
        cex.lab = 1.5, lty = 1)

# SMOOTHING ---------------------------------------------------------------

lambda.grid <- seq(2.5, 15, length = 20)

# hr
#  ciclo per scegliere il parametro di lisciamento
SmoothStats = array(NA, dim=c(length(lambda.grid), 3),
                        dimnames=list(lambda.grid, c("log10.lambda", "df", "gcv") ) )
SmoothStats[, 1] = lambda.grid
for (ilam in 1:length(lambda.grid)) {
  Smooth = smooth.basisPar(space_ob, HR, 4, Lfdobj=int2Lfd(2), lambda=10^(-lambda.grid[ilam]))
  # penalizziamo la derivata quarta per avere anche acelerazione liscia
  SmoothStats[ilam, "df"]  = Smooth$df
  SmoothStats[ilam, "gcv"] = sum(Smooth$gcv)
}
par(mfrow = c(1,2))
plot(SmoothStats[, c(1, 3)], type="b", log="y",lwd=2)
plot(SmoothStats[, 1:2], type="b", log="y",lwd=2)
par(mfrow = c(1,1))
lambda.opt <- 10^-lambda.grid[which.min(SmoothStats[,"gcv"])]
hr.fd <- smooth.basisPar(space_ob, HR, 6, Lfdobj=int2Lfd(4),lambda=lambda.opt)$fd
plot(hr.fd)

# scegliamo a mano un lambda un po' piu' grande
hr.fd <- smooth.basisPar(space_ob, HR, 6, Lfdobj=int2Lfd(4),lambda=0.00001)$fd
# funzione
plot(hr.fd)
abline(v = cumsum(c(1.3, 1.2, 1, 2.6)), lty = 2, lwd = 2) # punti di cambio pendenza
# derivata prima
plot(hr.fd, Lfd = 1)
abline(v = cumsum(c(1.3, 1.2, 1, 2.6)), lty = 2, lwd = 2) # punti di cambio pendenza
# derivata seconda
plot(hr.fd, Lfd = 2)
abline(v = cumsum(c(1.3, 1.2, 1, 2.6)), lty = 2, lwd = 2) # punti di cambio pendenza

# br
#  ciclo per scegliere il parametro di lisciamento
SmoothStats = array(NA, dim=c(length(lambda.grid), 3),
                        dimnames=list(lambda.grid, c("log10.lambda", "df", "gcv") ) )
SmoothStats[, 1] = lambda.grid
for (ilam in 1:length(lambda.grid)) {
  Smooth = smooth.basisPar(space_ob, BR,6, Lfdobj=int2Lfd(4), lambda=10^(-lambda.grid[ilam]))
  # penalizziamo la derivata quarta per avere anche acelerazione liscia
  SmoothStats[ilam, "df"]  = Smooth$df
  SmoothStats[ilam, "gcv"] = sum(Smooth$gcv)
}
par(mfrow = c(1,2))
plot(SmoothStats[, c(1, 3)], type="b", log="y",lwd=2)
plot(SmoothStats[, 1:2], type="b", log="y",lwd=2)
par(mfrow = c(1,1))
lambda.opt <- 10^-lambda.grid[which.min(SmoothStats[,"gcv"])]
br.fd <- smooth.basisPar(space_ob, BR, 6, Lfdobj=int2Lfd(4),lambda=lambda.opt)$fd
plot(br.fd)
# con il parametro scelto via covalida incrociata generalizzata si liscia troppo poco

br.fd <- smooth.basisPar(space_ob, BR, 6, Lfdobj=int2Lfd(4),lambda=0.001)$fd
# funzione
plot(br.fd)
abline(v = cumsum(c(1.3, 1.2, 1, 2.6)), lty = 2, lwd = 2) # punti di cambio pendenza
# derivata prima
plot(br.fd, Lfd = 1)
abline(v = cumsum(c(1.3, 1.2, 1, 2.6)), lty = 2, lwd = 2) # punti di cambio pendenza
# derivata seconda
plot(br.fd, Lfd = 2)
abline(v = cumsum(c(1.3, 1.2, 1, 2.6)), lty = 2, lwd = 2) # punti di cambio pendenza

# ANLISI ESPLORATIVA HR ---------------------------------------------------

# medie e varianze --------------------------------------------------------

# hr
# MEDIA
mean.hr <- mean.fd(hr.fd)
plot(mean.hr)
abline(v = cumsum(c(1.3, 1.2, 1, 2.6)), lty = 2, lwd = 2) 

plot(hr.fd)
lines(mean.hr, lwd = 2)

# BOXPLOT
fbp <- boxplot(hr.fd)

# VARIANZA
std.hr <- std.fd(hr.fd)
plot(std.hr)

# COVARIANZA
cov.hr <- var.fd(hr.fd)
grid <- seq(0, 6.1, length = 100)
varcov.M <- eval.bifd(grid, grid, cov.hr)
persp(grid, grid, varcov.M, theta = -45, phi = 25, r  = 3, expand = 0.5,
      ticktype = "detailed", xlab = "spazio", ylab = "spazio", 
      zlab = "Covarianza hr")
image(grid, grid, varcov.M)
# all'inizio c'è un trend crescente molto forte
# se non consideriamo il pezzo iniziale vediamo meglio la variabilità nei km successivi
grid <- seq(0.5, 6.1, length = 100)
varcov.M <- eval.bifd(grid, grid, cov.hr)
persp(grid, grid, varcov.M, theta = -45, phi = 25, r  = 3, expand = 0.5,
      ticktype = "detailed", xlab = "spazio", ylab = "spazio", 
      zlab = "Covarianza hr")
image(grid, grid, varcov.M)

# componenti principali ---------------------------------------------------

# hr
ncomp <- 4 
hr.pca <- pca.fd(hr.fd, ncomp)

par(mfrow = c(2,2))
plot(hr.pca)
par(mfrow = c(1,1))
# prima:
# riflette l'andamento generale: varianza circa costante dopo il primo tratto
# ma non cattura la variabilità iniziale
# seconda:
# cattura la variabilita iniziale e un po' nell'ultima fase (quella pianeggiante)
# le curve + e - si invertono: chi ha battiti piu alti della media
# all'inizio, può averli piu bassi nell'ultima fase
# terza:
# cattura la variabilità dalla fase in cui inizia la salita e nell'ultimo tratto 
# pianeggiante

# MODELLI PER I TEMPI DI ARRIVO -------------------------------------------

# y ~ hr ------------------------------------------------------------------

y <- int[,4]
y <- c(y)

# rappresentazione in basi del coefficiente funzionale
beta.basis <- create.bspline.basis(rangeval = c(0, 6.1), nbasis = 25)

# aggiungiamo l'intercetta al modello:
cbasis <- create.constant.basis(rangeval = c(0, 6.1)) 
cfd <- fd(coef = matrix(1,1,10), basisobj = cbasis)

xfdlist <- list(cfd, hr.fd)

# tuning di lambda (parametro di lisciamento per beta)
set.seed(42)
lambda <- seq(-3.5, -2, length = 100)
gocv <- matrix(NA, length(lambda), 3)
for(j in 1:length(lambda)){
  betalist <- list(fdPar(fdobj = cbasis, Lfdobj = 0, lambda = 0),
                   fdPar(fdobj = beta.basis, Lfdobj = 2, lambda = 10^lambda[j]))
  hr.reg <- fRegress(y, xfdlist, betalist)
  gocv[j,] <- c(lambda[j], hr.reg$GCV, hr.reg$OCV)
  if (j%%10 == 0) cat(round(j/length(lambda)*100),"%\n")
}

plot(gocv[,1:2], xlab = "", ylab = "", type = "l", lwd = 2)
abline(v = lambda[which.min(gocv[,2])], lty = 2, lwd = 2, col = 2)
lambdastar_ocv = 10^lambda[which.min(gocv[,2])]

# modello con il lambda scelto
betalist <- list(fdPar(fdobj = cbasis, Lfdobj = 0, lambda = 0),
                 fdPar(fdobj = beta.basis, Lfdobj = 2, lambda = lambdastar_ocv))
hr.reg <- fRegress(y, xfdlist, betalist)

# gradi di liberta
hr.reg$df 

# test di permutazione
Fpermres = Fperm.fd(y,xfdlist, betalist, plotres = FALSE)
Fpermres$pval
hist(Fpermres$Fnullvals, nclass = 50, main ="", ylab = "Frequenza")
abline(v = Fpermres$Fobs, lwd = 2, lty = 2, col = "red4")
abline(v = Fpermres$qval, lwd = 2, lty = 2, col = 3)
legend("topright", c("Statistica ossservata", "Quantile di livello 0.05"), col = c("red4",3), pch = 19)

# stima della varianza dell'errore
sigmae <- sum((y - hr.reg$yhatfdobj)^2)/(10 - hr.reg$df); sigmae

# R^2
Rsq <- 1-sigmae/var(y)
Rsq

# errori standard punto a punto per i coefficienti
hr.reg.std <- fRegress.stderr(hr.reg, y2cMap = NULL, SigmaE = sigmae*diag(rep(1,10)))

# intervallo di confidenza per l'intercetta
alpha.hat <- hr.reg$betaestlist[[1]]$fd$coefs
alpha.hat # coefficiente stimato per la base relativa all'intercetta
alphahaterr <- hr.reg.std$betastderrlist[[1]]$coefs
alphaci <- c(alpha.hat -2*alphahaterr, alpha.hat, alpha.hat + 2*alphahaterr)
alphaci

# intervallo di confidenza per beta(t)
betaest <- hr.reg$betaestlist[[2]]$fd
plot(betaest, ylab = expression(beta), ylim = c(-320, 100), xlab = "")
title("Stima del coefficiente funzionale di regressione")
lines(betaest + 2*hr.reg.std$betastderrlist[[2]], col = "red4", lty = 1)
lines(betaest - 2*hr.reg.std$betastderrlist[[2]], col = "red4", lty = 1)
abline(v = cumsum(c(1.3, 1.2, 1, 2.6)), lty = "dotted", lwd = 2, col = "grey70") 

# se il beta è positivo, all'aumentare di hr il tempo aumenta (arriva dopo)
# se il beta è negativo, all'aumentare di hr il tempo diminuisce (arriva prima)

# in salita: all'aumentare di hr arriva dopo
# in discesa:  all'aumentare di hr arriva prima
# alla fine del persorso: se hr aumenta il tempo diminuisce molto
# (battiti altissimi -> sprint finale)

# y ~ hr pca --------------------------------------------------------------

hr.pca <- pca.fd(hr.fd, 8)
hr.mpc.step <- step(lm(y~., data = as.data.frame(hr.pca$scores)), direction = "both", trace = 1)
summary(hr.mpc.step)
# teniamo solo le cp significative
hr.mpc <- lm(y ~ hr.pca$scores[,c(4,5,6)])
summary(hr.mpc)

pcabeta = hr.mpc$coef[2]*hr.pca$harmonics[4]+
  hr.mpc$coef[3]*hr.pca$harmonics[5] + 
  hr.mpc$coef[4]*hr.pca$harmonics[6]

pcacoefvar = summary(hr.mpc)$coefficients[,2]^2
pcabetavar = pcacoefvar[2]*hr.pca$harmonics[4]^2+
  pcacoefvar[3]*hr.pca$harmonics[5]^2 +
  pcacoefvar[4]*hr.pca$harmonics[6]^2

pcabeta1<- predict(pcabeta, space_ob)
pcabetavar1<- predict(pcabetavar, space_ob)

# grafico del coefficiente stimato
plot(space_ob,pcabeta1, type="l", col=1, ylab=expression(beta), xlab="km")
lines(space_ob,pcabeta1+2*sqrt(pcabetavar1), col = "red4", lty = 2)
lines(space_ob,pcabeta1-2*sqrt(pcabetavar1), col = "red4", lty = 2)
abline(h=0)

# fitting del modello 
hr.pca.pred <- predict(hr.mpc)
plot(y, hr.pca.pred)
abline(0,1)
sigmae.hr_pca =  sum((y - hr.pca.pred)^2)/(10-hr.mpc$df)
sqrt(sigmae.hr_pca)
1 - sigmae.hr_pca/var(y)

# confronto con il beta stimato col modello funzionale
par(mfrow = c(1,2))
betaest <- hr.reg$betaestlist[[2]]$fd
plot(betaest, ylab = expression(beta), ylim=c(-100,100))
lines(betaest + 2*hr.reg.std$betastderrlist[[2]], col = "red4", lty = 1)
lines(betaest - 2*hr.reg.std$betastderrlist[[2]], col = "red4", lty = 1)
abline(v = cumsum(c(1.3, 1.2, 1, 2.6)), lty = "dotted", lwd = 2, col = "grey70")

plot(space_ob, pcabeta1, type="l",ylab=expression(beta), xlab="",ylim=c(-100,100))
abline(h = 0, lty = 2)
lines(space_ob,pcabeta1+2*sqrt(pcabetavar1), lty=1, col="red4", lwd=1)
lines(space_ob,pcabeta1-2*sqrt(pcabetavar1), lty=1, col="red4", lwd=1)
abline(v = cumsum(c(1.3, 1.2, 1, 2.6)), lty = "dotted", lwd = 2, col = "grey70")
par(mfrow = c(1,1))
# andameti piu o meno simili

# y ~ hr flirti -----------------------------------------------------------

setwd("C:/Users/gaiom/Desktop/FUNZIONALI/lab/5 x funzionale")
source("flrti.txt")

# regressione linerare funzionale interpretabile
m.flirti.cv <- flrti.cv(Y = matrix(y, 10, 1), X = scale(t(HR)),
                  sigma = seq(0.001, 0.003, length = 11), deriv = 2,
                  weight = seq(0.1, 1, length = 10), cv = 10)
m.flirti <- flrti(matrix(y, 10, 1), scale(t(HR)), sigma = 0.001, weight = 0.1, deriv = 3)

str(m.flirti)
plot(space_ob, m.flirti$beta, type = "l")
abline(v = cumsum(c(1.3, 1.2, 1, 2.6)), lty = "dotted", lwd = 2, col = "grey70")
abline(h = 0, col = 2, lty = 2)


# Modello simultaneo hr ~ br + speed --------------------------------------

temp_df = data.frame(smoke = as.factor(covariate$Smoker))
temp_df$BR = t(BR)
temp_df$SPEED = t(SPEED)
temp_df$HR = t(HR)
fit3 = pffr(HR~BR+SPEED, data=temp_df)
rboot <- coefboot.pffr(fit3, B=200)
plot(rboot$smterms$`BR(yindex)`$value, type = "l", ylim = c(-2,3), xlab="Spazio", ylab=expression(beta),
     col="red4", lwd=3)
lines(rboot$smterms$`BR(yindex)`[,4], col="red4", lwd=3, lty=2)
lines(rboot$smterms$`BR(yindex)`[,6], col="red4", lwd=3, lty=2)
opar <- par()
par(mar = c(5, 4, 0.3, 0.5))
plot(fit3, select=2, col="red4", lwd=3, rug=F, scale=0, ylab=expression(beta), xlab = "Spazio")
