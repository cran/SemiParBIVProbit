gt.bpm <- function(object) {

if(object$BivD!="N") stop("This test's implementation is valid for bivariate normal models only.")

corr <- 0

eta1 <- object$X1%*%coef(object$gam1)
eta2 <- object$X2%*%coef(object$gam2)

y1 <- object$y1
y2 <- object$y2

y1.y2 <- y1 * y2
y1.cy2 <- y1 * (1 - y2)

p11 <- pmax(abs(pbinorm(eta1, eta2, cov12 = corr)), 1000 * .Machine$double.eps)
p10 <- pmax(pnorm(eta1) - p11, 1000 * .Machine$double.eps)

if(object$sel==FALSE){
cy1.y2 <- (1 - y1) * y2
cy1.cy2 <- (1 - y1) * (1 - y2)
p01 <- pmax(pnorm(eta2) - p11, 1000 * .Machine$double.eps)
p00 <- pmax(1 - p11 - p10 - p01, 1000 * .Machine$double.eps)
}

d.n1n2 <- dbinorm(eta1, eta2, cov12 = corr) 

if(object$sel==FALSE) dl.drho <- d.n1n2*(y1.y2/p11 - y1.cy2/p10 - cy1.y2/p01 + cy1.cy2/p00) else dl.drho <- d.n1n2*(y1.y2/p11 - y1.cy2/p10)

G <- as.numeric(sum(dl.drho)*object$rho)
return(pchisq(G, 1, lower.tail = FALSE))

}
