ass.dp <- function(ass.s, BivD, scc, sccn, nCa){

if(BivD %in% scc)  ass.s <-  abs(ass.s)   
if(BivD %in% sccn) ass.s <- -abs(ass.s) 

if(!(BivD %in% c("AMH","FGM","PL","HO"))) i.rho <- BiCopTau2Par(family = nCa, tau = ass.s)
if(  BivD %in% c("AMH","FGM") )           i.rho <- BiCopTau2Par(family = 1,   tau = ass.s)
if(  BivD %in% c("PL") )                  i.rho <- BiCopTau2Par(family = 3,   tau = ass.s)
if(  BivD %in% c("HO") )                  i.rho <- ass.s - 1

if(BivD %in% c("N","AMH","FGM","T"))         i.rho <- atanh( i.rho ) 
if(BivD == "F")                              i.rho <- ifelse( abs(i.rho) < 0.0000001, 0.0000001, i.rho ) 
if(!(BivD %in% c("N","AMH","FGM","F","T")))  i.rho <- abs(i.rho)

if(BivD %in% c("C0","C180","C90","C270",
               "C0C90","C0C270","C180C90","C180C270"))             i.rho <-  log(i.rho)   
if(BivD %in% c("J0","J180","G0","G180","J90","J270","G90","G270",
               "J0J90","J0J270","J180J90","J180J270",
               "G0G90","G0G270","G180G90","G180G270"))             i.rho <-  log(i.rho - 1)   

names(i.rho) <- "theta.star"   

i.rho

}

