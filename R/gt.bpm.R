gt.bpm <- function(object) {

if(object$BivD!="N" || object$PL!="P") stop("This test's implementation is currently valid for bivariate normal errors only.")

corr <- 0

eta1 <- object$gam1$linear.predictors
eta2 <- object$VC$X2 %*% coef(object$gam2) 

p11 <- pmax(pbinorm(eta1, eta2, cov12 = corr), 1000 * .Machine$double.eps)

if(object$Model=="B" || object$Model=="BSS") p10 <- pmax(pnorm(eta1) - p11, 1000 * .Machine$double.eps)

if(object$Model=="B"){ p01 <- pmax(pnorm(eta2) - p11, 1000 * .Machine$double.eps)
                       p00 <- pmax(1 - p11 - p10 - p01, 1000 * .Machine$double.eps)
                      }

d.n1n2 <- dbinorm(eta1, eta2, cov12 = corr) 

if(object$Model=="B")   dl.drho <- object$weights*d.n1n2*(object$respvec$y1.y2/p11 - object$respvec$y1.cy2/p10 - object$respvec$cy1.y2/p01 + object$respvec$cy1.cy2/p00) 
if(object$Model=="BSS") dl.drho <- object$weights*d.n1n2*(object$respvec$y1.y2/p11 - object$respvec$y1.cy2/p10)
if(object$Model=="BPO") dl.drho <- object$weights*d.n1n2*(object$y1/p11 - (1-object$y1)/(1-p11))  

G <- as.numeric(sum(dl.drho)*object$rho)
return(pchisq(G, 1, lower.tail = FALSE))

}
