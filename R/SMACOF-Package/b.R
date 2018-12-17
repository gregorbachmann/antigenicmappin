n <- c(10, 20, 30, 50, 100)
ndim <- c(1, 2, 3, 4, 5)
nrep <- 500

## --- ratio MDS
stressmat1 <- matrix(NA, ncol = length(ndim), nrow = length(n), dimnames = list(n, ndim))
for (i in 1:length(ndim)) {
  for (j in 1:length(n)) {
    cat(ndim[i], "/", n[j], "\n")
    stressmat1[j, i] <- mean(randomstress(n = n[j], ndim = ndim[i], nrep = nrep, type = "ratio"))    
  }
}
stressmat1


## --- interval MDS
stressmat2 <- matrix(NA, ncol = length(ndim), nrow = length(n), dimnames = list(n, ndim))
for (i in 1:length(ndim)) {
  for (j in 1:length(n)) {
    cat(ndim[i], "/", n[j], "\n")
    stressmat2[j, i] <- mean(randomstress(n = n[j], ndim = ndim[i], nrep = nrep, type = "interval"))    
  }
}
stressmat2

## --- ordinal MDS
stressmat3 <- matrix(NA, ncol = length(ndim), nrow = length(n), dimnames = list(n, ndim))
for (i in 1:length(ndim)) {
  for (j in 1:length(n)) {
    cat(ndim[i], "/", n[j], "\n")
    stressmat3[j, i] <- mean(randomstress(n = n[j], ndim = ndim[i], nrep = nrep, type = "ordinal")) 
  }
}
stressmat3

stressnorms <- list(stressmat1 = stressmat1, stressmat2 = stressmat2, stressmat3 = stressmat3)
## These stress norms can be subject to plotting (cf. Fig. 1)

## ---------------------------------- MDS Examples -----------------------------------
## --- Lawler dataset
Lawler                                                 ## similarity matrix (correlations) from Lawler (1967)
LawlerD <- sim2diss(Lawler)                            ## convert into dissimilarities
fitlaw1 <- mds(LawlerD, type = "ordinal")              ## ordinal MDS fit 
fitlaw2 <- mds(LawlerD, type = "interval")             ## interval MDS fit

## Figure 2 in paper
op <- par(mfrow = c(2,2))
plot(fitlaw1, main = "Configuration Plot Ordinal MDS")
plot(fitlaw1, plot.type = "Shepard", main = "Shepard Diagram Ordinal MDS")
plot(fitlaw2, main = "Configuration Plot Interval MDS")
plot(fitlaw2, plot.type = "Shepard", main = "Shepard Diagram Interval MDS")
par(op)
## The ordinal solution is degenerate

set.seed(123)                                         ## set a seed such to reproduce permutation results
fitlaw1 <- mds(LawlerD)                               ## ratio MDS 
fitlaw1
res.perm1 <- permtest(fitlaw1, nrep = 500)            ## ratio MDS fit on 500 permuted dissimilarity matrices
res.perm1                                             ## can't reject H0
fitlaw2 <- mds(LawlerD, type = "interval")            ## interval MDS 
fitlaw2
res.perm2 <- permtest(fitlaw2, nrep = 500)            ## interval MDS fit on 500 permuted dissimilarity matrices
res.perm2                                             ## reject H0

## --- Wenchuan dataset
data(Wenchuan)
Wdelta <- dist(t(Wenchuan))                           ## Euclidean distances                   
fit.wenchuan1 <- mds(Wdelta, type = "interval",init="random")        ## ordinal MDS fit  
plot(fit.wenchuan1)
fit.wenchuan1

## scree plot
stressvec <- NULL
for (i in 1:16) stressvec[i] <- mds(Wdelta, ndim = i, type = "ordinal")$stress
plot(1:16, stressvec, pch = 20, type = "b", xlab = "Number of Dimensions", ylab = "Stress-1", main = "MDS Scree Plot", xaxp = c(1, 16, 15))

## multiple random starts
set.seed(123)
fit.wenchuan <- NULL  
for(i in 1:100) fit.wenchuan[[i]] <- mds(Wdelta, type = "ordinal", init = "random") 
ind <- which.min(sapply(fit.wenchuan, function(x) x$stress))
fit.wenchuan2 <- fit.wenchuan[[ind]]
sppwen <- sort(fit.wenchuan2$spp, decreasing = TRUE)

fit.wenchuan2a <- mds(Wdelta, type = "interval", init = fit.wenchuan2$init)  ## interval MDS
fit.wenchuan2a
fit.wenchuan2b <- mds(Wdelta, type = "mspline", init = fit.wenchuan2$init)   ## spline MDS
fit.wenchuan2b

## Shepard diagrams
op <- par(mfrow = c(1,3))
plot(fit.wenchuan2, plot.type = "Shepard", main = "Shepard Diagram Ordinal", ylim = c(0, 2.2))
plot(fit.wenchuan2a, plot.type = "Shepard", main = "Shepard Diagram Interval", ylim = c(0, 2.2))
plot(fit.wenchuan2b, plot.type = "Shepard", main = "Shepard Diagram Spline", ylim = c(0, 2.2))
par(op)

## SPP
Wdelta2 <- as.matrix(Wdelta)
ind <- which(rownames(Wdelta2) %in% c("lossint", "future", "dreams"))  ## remove points with high SPP
Wdelta3 <- Wdelta2[-ind, -ind]
fit.wenchuan3 <- mds(Wdelta3, type = "ordinal")
fit.wenchuan3

## bubble plot
plot(fit.wenchuan2, plot.type = "bubbleplot", main = "Bubble Plot Wenchuan", asp = 1)

## permutation test (on the raw data matrix)
set.seed(123)
nrep <- 500                                        
res.perm <- permtest(fit.wenchuan2, data = Wenchuan, method.dat = "euclidean", nrep = nrep)
res.perm
mperm <- mean(res.perm$stressvec)                      ## permutation stress norm
mperm
perm5 <- quantile(res.perm$stressvec, probs = 0.05)    ## lower 5% quantile (critical value)
perm5

hist(res.perm$stressvec, xlim = c(0.10, 0.40), xlab = "Stress Values", main = "Wenchuan Permutation Histogram")
abline(v = perm5, lty = 2)             ## critical value (dashed)
abline(v = fit.wenchuan2$stress)       ## observed stress value

## jackknife
resjack <- jackknife(fit.wenchuan2)
plot(resjack, main = "Wenchuan MDS Jackknife")

## final configuration plot with DSM-IV facets
colpal <- c(rainbow_hcl(3, c = 100))
pal <- palette(colpal)
memb <- c(1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3) 
plot(fit.wenchuan2, main = "Wenchuan Configurations", label.conf = list(label = FALSE), col = memb)
legend("topright", legend = c("intrusive recollection", "avoidance/numbing", "arousal"), col = colpal, pch = 20)
abline(-0.05, 0.2, col = "lightgray")
lines(c(0.1, -1),c(-0.03, 1), col = "lightgray")
palette(pal)
colpal <- c(rainbow_hcl(3, c = 100, l = 30))
pal <- palette(colpal)
text(fit.wenchuan2$conf[-7,], labels = rownames(fit.wenchuan2$conf)[-7], col = memb[-7], pos = 3, cex = 0.8)
text(fit.wenchuan2$conf[7,1], fit.wenchuan2$conf[7,2], labels = 
       rownames(fit.wenchuan2$conf)[7], col = memb[7], pos = 1, cex = 0.8)
palette(pal)
