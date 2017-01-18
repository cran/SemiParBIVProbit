ass.ms <- function(BivD, nCa, theta){

if(BivD %in% c("J0","J180","J90","J270"))  theta <- ifelse(abs(theta) > 50, 50, abs(theta))
if(BivD %in% c("J90","J270"))              theta <- -theta 

if(BivD %in% c("F")){ signs <- sign(theta) 
                      theta <- ifelse(abs(theta) > 100, 100, abs(theta))
                      theta <- theta*signs }

if(BivD %in% c("C0","C180","G0","G180","C90","C270","G90","G270")) theta <- ifelse(abs(theta) > 100, 100, abs(theta))
if(BivD %in% c("C90","C270","G90","G270"))                         theta <- -theta

if(!(BivD %in% c("AMH","FGM"))) tau <- BiCopPar2Tau(family = nCa, par = theta)
if(BivD == "AMH")               tau <- 1 - (2/3)/theta^2*(theta + (1-theta)^2*log(1-theta))
if(BivD == "FGM")               tau <- 2/9*theta 

tau.a   <- mean(tau) 
theta.a <- mean(theta) 

list(theta.a = theta.a, theta = theta, tau = tau, tau.a = tau.a)

}
