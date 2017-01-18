bcorrec <- function(VC, esnl, lp){

###############################
# fixed quantities
###############################

n  <- VC$n
rc <- VC$rc
rsim <- VC$rsim

int1 <- int2 <- int3 <- NA
H <- XXX <- XXXX <- BH <- list()
G <- matrix(NA, n, lp)

if(VC$margins[2] == "ZTP") rZTP <- function(n, mu) qpois(runif(n, dpois(0, mu), 1), mu)

if(is.null(VC$X2)) X2 <- matrix(1, n, 1) else X2 <- VC$X2 
# need to generalise to three parameter case (the whole code below as well) 

# need to check y1m for discrete margin as this would have to be changed I guess
###############################

for(i in 1:n){

set.seed(1) # yes or no? we would use this to make sure the results are reproducible
if(VC$margins[2] == "N")     y <- rNO(   rsim,    mu =     esnl$eta[i],     sigma = sqrt(esnl$sigma2[i])) 
if(VC$margins[2] == "N2")    y <- rNO(   rsim,    mu =     esnl$eta[i],     sigma =      esnl$sigma2[i]) 
if(VC$margins[2] == "GU")    y <- rGU(   rsim,    mu =     esnl$eta[i],     sigma = sqrt(esnl$sigma2[i])) 
if(VC$margins[2] == "rGU")   y <- rRG(   rsim,    mu =     esnl$eta[i],     sigma = sqrt(esnl$sigma2[i])) 
if(VC$margins[2] == "LO")    y <- rLO(   rsim,    mu =     esnl$eta[i],     sigma = sqrt(esnl$sigma2[i])) 
if(VC$margins[2] == "LN")    y <- rLOGNO(rsim,    mu =     esnl$eta[i],     sigma = sqrt(esnl$sigma2[i])) 
if(VC$margins[2] == "WEI")   y <- rWEI(  rsim,    mu = exp(esnl$eta[i]),    sigma = sqrt(esnl$sigma2[i])) 
if(VC$margins[2] == "iG")    y <- rIG(   rsim,    mu = exp(esnl$eta[i]),    sigma = sqrt(esnl$sigma2[i])) 
if(VC$margins[2] == "GA")    y <- rGA(   rsim,    mu = exp(esnl$eta[i]),    sigma = sqrt(esnl$sigma2[i])) 
if(VC$margins[2] == "GAi")   y <- rGA(   rsim,    mu =     esnl$eta[i],     sigma = sqrt(esnl$sigma2[i])) 
if(VC$margins[2] == "DAGUM") y <- rGB2(  rsim,    mu = exp(esnl$eta[i]),    sigma = sqrt(esnl$sigma2[i]), nu = esnl$nu[i], tau = 1 ) 
if(VC$margins[2] == "SM")    y <- rGB2(  rsim,    mu = exp(esnl$eta[i]),    sigma = sqrt(esnl$sigma2[i]), nu = 1 , tau = esnl$nu[i]) 
if(VC$margins[2] == "BE")    y <- rBE(   rsim, mu = plogis(esnl$eta[i]),    sigma = sqrt(esnl$sigma2[i])) 
if(VC$margins[2] == "FISK")  y <- rGB2(  rsim,    mu = exp(esnl$eta[i]),    sigma = sqrt(esnl$sigma2[i]), nu = 1 , tau = 1 )
if(VC$margins[2] == "NBI")   y <- rNBI(  rsim,    mu = exp(esnl$eta[i]),    sigma = sqrt(esnl$sigma2[i])) 
if(VC$margins[2] == "NBII")  y <- rNBII( rsim,    mu = exp(esnl$eta[i]),    sigma = sqrt(esnl$sigma2[i])) 
if(VC$margins[2] == "PIG")   y <- rPIG(  rsim,    mu = exp(esnl$eta[i]),    sigma = sqrt(esnl$sigma2[i])) 
if(VC$margins[2] == "PO")    y <- rPO(   rsim,    mu = exp(esnl$eta[i])) 
if(VC$margins[2] == "ZTP")   y <- rZTP(  rsim,    mu = exp(esnl$eta[i])) 
rm(list=".Random.seed", envir=globalenv()) 

#print(y)
#print(c(esnl$eta[i], sqrt(esnl$sigma2[i])) )

if(VC$margins[2] %in% VC$m2)            dHs  <-      distrHs(y, esnl$eta[i], esnl$sigma2[i], esnl$sigma2.st[i], nu = 1, nu.st = 1, margin2=VC$margins[2], naive = TRUE)
if(VC$margins[2] %in% c(VC$m1d,VC$m2d)) dHs  <- distrHsDiscr(y, esnl$eta[i], esnl$sigma2[i], esnl$sigma2.st[i], nu = 1, nu.st = 1, margin2=VC$margins[2], naive = TRUE, y2m = VC$y1m)

pdf2                         <- dHs$pdf2
derpdf2.dereta2              <- dHs$derpdf2.dereta2 
derpdf2.dersigma2.st         <- dHs$derpdf2.dersigma2.st  
der2pdf2.dereta2             <- dHs$der2pdf2.dereta2
der2pdf2.dersigma2.st2       <- dHs$der2pdf2.dersigma2.st2         
der2pdf2.dereta2dersigma2.st <- dHs$der2pdf2.dereta2dersigma2.st  

comp1   <- 1 + exp(log(pdf2) + rc) 
comp2   <- pdf2/comp1
comp3   <- pdf2/comp1^2

int1[i] <- mean( log(comp1)/ pdf2 )





#integ <- NA
#
#int1f <- function(y) log( 1 + exp( log( dnorm(y, mu, sqrt(sig2)) ) + rc ) )
#
#for(i in 1:n){
#
#mu   <- esnl$eta[i]
#sig2 <- esnl$sigma2[i]
#
#integ[i] <- integrate(int1f, -Inf, Inf)
#
#}
#
#n - exp(-rc)*sum(integ)




##########

dl.dbe       <- derpdf2.dereta2/pdf2
dl.dsigma.st <- derpdf2.dersigma2.st/pdf2

G[i, ] <- c(   colMeans( t(t(comp2*dl.dbe/pdf2))%*%VC$X1[i,] ), colMeans( t(t(comp2*dl.dsigma.st/pdf2))%*%X2[i,] )    )


#eta2 <- esnl$eta[i]
#sigma2 <- esnl$sigma2[i]
#sigma2.st <- esnl$sigma2.st[i]
#X1 <- VC$X1
#
#int2f <- function(y2){ 
#
#pdf2 <-  dnorm(y2, eta2, sqrt(sigma2)) 
#
#epsilon <- 0.0000001 
#pdf2 <- ifelse(pdf2 < epsilon, epsilon, pdf2 )
#
#dermu2.dereta2  <- 1
#derpdf2.dermu2  <- (1/(sqrt(2 * pi * sigma2))) * (exp(-0.5 * (y2 - eta2)^2/sigma2) * (0.5 * (2 * (y2 - eta2))/sigma2))  
#derpdf2.dereta2 <- derpdf2.dermu2*dermu2.dereta2  
#
#
#dersigma2.dersigma2.st <- exp(sigma2.st)  
#derpdf2.sigma2  <- (2 * ((y2 - eta2)^2/(2 * sigma2)^2) - 1/(2 * sigma2)) * exp(-((y2 - eta2)^2/(2 * sigma2)))/sqrt(2 * (pi * sigma2)) 
#derpdf2.dersigma2.st  <- derpdf2.sigma2 * dersigma2.dersigma2.st   
#
#
#dl.dbe       <- derpdf2.dereta2/pdf2
#dl.dsigma.st <- derpdf2.dersigma2.st/pdf2
#
#comp1   <- 1 + exp(log(pdf2) + rc) 
#comp2   <- pdf2/comp1
#
##comp2*dl.dbe*X1[i, 3]
#
#comp2*dl.dsigma.st
#
#
#}
#
#
#integrate(int2f, -Inf, Inf)
#



 
##########






#eta2 <- esnl$eta[i]
#sigma2 <- esnl$sigma2[i]
#sigma2.st <- esnl$sigma2.st[i]
#X1 <- VC$X1
#
#int3f <- function(y2){ 
#
#pdf2 <-  dnorm(y2, eta2, sqrt(sigma2)) 
#
#epsilon <- 0.0000001 
#pdf2 <- ifelse(pdf2 < epsilon, epsilon, pdf2 )
#
#dermu2.dereta2  <- 1
#der2mu2.dereta2eta2    <- 0 
#
#derpdf2.dermu2  <- (1/(sqrt(2 * pi * sigma2))) * (exp(-0.5 * (y2 - eta2)^2/sigma2) * (0.5 * (2 * (y2 - eta2))/sigma2))  
#derpdf2.dereta2 <- derpdf2.dermu2*dermu2.dereta2  
#
#
#dersigma2.dersigma2.st <- exp(sigma2.st)  
#derpdf2.sigma2  <- (2 * ((y2 - eta2)^2/(2 * sigma2)^2) - 1/(2 * sigma2)) * exp(-((y2 - eta2)^2/(2 * sigma2)))/sqrt(2 * (pi * sigma2)) 
#derpdf2.dersigma2.st  <- derpdf2.sigma2 * dersigma2.dersigma2.st   
#
#
#dl.dbe       <- derpdf2.dereta2/pdf2
#dl.dsigma.st <- derpdf2.dersigma2.st/pdf2
#
#comp1   <- 1 + exp(log(pdf2) + rc) 
#comp2   <- pdf2/comp1
#
##comp2*dl.dbe*X1[i, 3]
#
#der2pdf2.dermu2       <- ((y2 - eta2)^2/sigma2 - 1) * exp(-(0.5 * ((y2 - eta2)^2/sigma2)))/(sigma2 * sqrt(2 * (pi * sigma2)))
#
#
#der2pdf2.dereta2 <- der2pdf2.dermu2* dermu2.dereta2^2 + derpdf2.dermu2*der2mu2.dereta2eta2        
#
#
#d2l.be.be        <- (der2pdf2.dereta2*pdf2 - (derpdf2.dereta2)^2)/pdf2^2 
#
#d2l.be.be * comp2
#
#}
#
#
#integrate(int3f, -Inf, Inf)
#
#
#
#


d2l.be.be        <- (der2pdf2.dereta2*pdf2 - (derpdf2.dereta2)^2)/pdf2^2 
d2l.sigma.sigma  <- (der2pdf2.dersigma2.st2*pdf2-(derpdf2.dersigma2.st)^2)/pdf2^2 
d2l.be.sigma     <- (der2pdf2.dereta2dersigma2.st*pdf2 - derpdf2.dereta2*derpdf2.dersigma2.st)/pdf2^2 

for(j in 1:rsim){ # this operation can be vectorised # this part needs to be verified

XX <- c(VC$X1[i, ]*dl.dbe[j], X2[i, ]*dl.dsigma.st[j])
XXX[[j]] <- (XX%*%t(XX))*(comp3/pdf2)[j]

}





XXXX[[i]] <- Reduce("+", XXX)/length(XXX)


H1 <- apply(  outer(VC$X1[i, ]%*%t(VC$X1[i, ]), d2l.be.be * comp2/pdf2, FUN = "*"),        c(1, 2), FUN = mean)
H2 <- apply(  outer(   X2[i, ]%*%t(   X2[i, ]), d2l.sigma.sigma * comp2/pdf2, FUN = "*"),  c(1, 2), FUN = mean)
H3 <- apply(  outer(VC$X1[i, ]%*%t(   X2[i, ]), d2l.be.sigma * comp2/pdf2, FUN = "*"),     c(1, 2), FUN = mean)

H[[i]] <- rbind( cbind( H1      , H3  ), 
                 cbind( t(H3)   , H2  ) )  


BH[[i]] <- XXXX[[i]] + H[[i]]   


}


b <- n - exp(-rc)*sum(int1)
bp <- - colSums(G) 
bs <- - Reduce("+", BH) 

#print("end")

list(b = b, bp = bp, bs = bs)


}    