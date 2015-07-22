\name{hiv}
\alias{hiv}
\alias{hiv.polys}
\docType{data}


\title{HIV Zambian data}

\description{HIV Zambian data by region, together with polygons describing the regions' shapes.}

\usage{
data(hiv)
data(hiv.polys)
}



\format{ \code{hiv} is a 5000 row data frame with the following columns 



\describe{


\item{hivconsent}{binary variable indicating consent to test for HIV.}
\item{hiv}{binary variable indicating whether an individual is HIV positive.}
\item{age}{age in years.}
\item{education}{years of education.}
\item{region}{code identifying region, and matching \code{names(hiv.polys)}. It can take nine possible values: 1 central, 2 copperbelt, 3 eastern, 
             4 luapula, 5 lusaka, 6 northwestern, 7 northern, 8 southern, 9 western.}
\item{marital}{never married, currently married, formerly married.}
\item{std}{had a sexually transmitted disease.}
\item{highhiv}{had high risk sex.}
\item{condom}{used condom during last intercourse.}
\item{aidscare}{equal to 1 if would care for an HIV-infected relative.}
\item{knowsdiedofaids}{equal to 1 if know someone who died of HIV.}
\item{evertestedHIV}{equal to 1 if previously tested for HIV.}           
\item{smoke}{smoker.}         
\item{et}{bemba, lunda (luapula), lala, ushi, lamba, tonga, luvale,
lunda (northwestern), mbunda, kaonde, lozi, chewa, nsenga,
ngoni, mambwe, namwanga, tumbuka, other.}              
\item{language}{English, Bemba, Lozi, Nyanja, Tonga, other.}                     
\item{interviewerID}{interviewer identifier.}                        
\item{sw}{survey weights.}  

}



\code{hiv.polys} contains the polygons defining the areas in the format described below.
} 


\details{The data frame \code{hiv} relates to the regions whose boundaries are coded in \code{hiv.polys}.
\code{hiv.polys[[i]]} is a 2 column matrix, containing the vertices of the polygons defining the boundary of the ith 
region. \code{names(hiv.polys)} matches \code{hiv$region} (order unimportant).

}

\source{

The data have been produced using the procedure described in: 

McGovern M.E., Barnighausen T., Marra G. and Radice R. (2015), On the Assumption of Joint Normality in Selection Models: A Copula Approach Applied to Estimating HIV Prevalence. \emph{Epidemiology}, 26(2), 229-237.


}

\examples{

\dontrun{ 

#########################################################
#########################################################

library("SemiParBIVProbit")

data("hiv", package = "SemiParBIVProbit")        
data("hiv.polys", package = "SemiParBIVProbit")  

xt <- list(polys = hiv.polys) 
# neighbourhood structure info for MRF  
            
#########################################################
# Bivariate brobit model with non-random sample selection
#########################################################          
            
sel.eq <- hivconsent ~ s(age) + s(education) + 
                       s(region, bs = "mrf", xt = xt) + marital + std + 
                       highhiv + condom + aidscare + knowsdiedofaids + 
                       evertestedHIV + smoke + et + language + 
                       s(interviewerID, bs = "re") 
 
out.eq <- hiv ~ s(age) + s(education) + s(region, bs = "mrf", xt = xt) +  
                marital + std + highhiv + condom + aidscare + 
                knowsdiedofaids + evertestedHIV + smoke  + et +  
                language     
                       
theta.eq <- ~ s(region, bs = "mrf", xt = xt)                        
  

bss <- SemiParBIVProbit(list(sel.eq, out.eq, theta.eq), 
                        data = hiv, BivD = "J90", Model = "BSS")

conv.check(bss)

set.seed(1)
summary(bss)

par(mfrow = c(2, 3), mar=c(4, 4.7, 3.1, 2), 
    cex.lab = 1.6, cex.axis = 1.6) 
for (i in 1:3) plot(bss, eq = 1, seWithMean = TRUE, scheme = 1,   
                    scale = 0, select = i, jit = TRUE)
for (i in 1:3) plot(bss, eq = 2, seWithMean = TRUE, scheme = 1,
                    scale = 0, select = i, jit = TRUE)

#dev.copy(postscript, "smoothHIV.eps", height=7)
#dev.off()  

prev(bss, sw = hiv$sw, type = "naive") 

set.seed(1)
prev(bss, sw = hiv$sw, type = "univariate") 

prev(bss, sw = hiv$sw) 


lr <- length(hiv.polys) 
prevBYreg  <- matrix(NA, lr, 2)
thetaBYreg <- NA

for(i in 1:lr) {
prevBYreg[i,1] <- prev(bss, sw = hiv$sw, ind = hiv$region==i, 
                       type = "univariate")$res[2]
prevBYreg[i,2] <- prev(bss, sw = hiv$sw, ind = hiv$region==i)$res[2]
thetaBYreg[i]  <- bss$theta[hiv$region==i][1]
}


zlim <- range(prevBYreg)  # to establish a common prevalence range

par(mfrow = c(1, 3), cex.axis = 1.3)
polys.map(hiv.polys, prevBYreg[,1], zlim = zlim, lab = "",  
          cex.lab = 1.5, cex.main = 1.5, 
          main = "HIV - Imputation Model")
polys.map(hiv.polys, prevBYreg[,2], zlim = zlim, cex.main = 1.5, 
          main = "HIV - Selection Model")
polys.map(hiv.polys, thetaBYreg, rev.col = FALSE, cex.main = 1.7, 
          main = expression(paste("Copula parameter (",hat(theta),")")))

#dev.copy(postscript, "PrevMaps.eps", height=3.5)
#dev.off() 

set.seed(1)
CItheta <- summary(bss)$CItheta
CItheta[1,]

}

#

}





