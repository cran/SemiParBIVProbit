post.check <- function(x, main = "Density Estimate",
                           xlab = "Quantile Residuals", ...){

pp <- distrHsAT(x$y2, x$eta2, x$sigma2, x$nu, x$margins[2]) 
p <- pp$p2
d <- pp$pdf2

par(mfrow = c(1, 2))
qqp <- qnorm(p)
plot(density(qqp, adjust = 2), lwd = 2, main = main, xlab = xlab, type = "l", ylab = "Density")
rug(qqp)

qqnorm(qqp)
abline(0, 1, col = "red")

}

