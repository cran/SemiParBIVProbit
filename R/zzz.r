#print.SemiParBIVProbit.version <- function()
#{ library(help=SemiParBIVProbit)$info[[1]] -> version
#  version <- version[pmatch("Version",version)]
#  um <- strsplit(version," ")[[1]]
#  version <- um[nchar(um)>0][2]
#  hello <- paste("\nThis is SemiParBIVProbit ",version,".\nFor overview type 'help(\"SemiParBIVProbit-package\")'.\n",sep="")
#  packageStartupMessage(hello)
#}

.onAttach <- function(...) { 

#require(MASS, quietly = TRUE, warn.conflicts = FALSE)
#require(VGAM, quietly = TRUE, warn.conflicts = FALSE)
#require(magic, quietly = TRUE, warn.conflicts = FALSE)
#require(statmod, quietly = TRUE, warn.conflicts = FALSE)
#require(mgcv, quietly = TRUE, warn.conflicts = FALSE)
#require(trust, quietly = TRUE, warn.conflicts = FALSE)
#require(mvtnorm, quietly = TRUE, warn.conflicts = FALSE)

  library(help=SemiParBIVProbit)$info[[1]] -> version
  version <- version[pmatch("Version",version)]
  um <- strsplit(version," ")[[1]]
  version <- um[nchar(um)>0][2]
  hello <- paste("\nThis is SemiParBIVProbit ",version,".\nFor overview type 'help(\"SemiParBIVProbit-package\")'.\n",sep="")
  packageStartupMessage(hello)
  
}






