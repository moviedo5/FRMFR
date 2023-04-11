globalVariables('icnt')

#' @keywords internal
 

#' @import fda 
#' @import splines
#' @import MASS
#' @import mgcv
#' @import nlme
#' @import methods
#' @import utils
#' @import grDevices
# @import graphics
#' @import stats
  
  # @import rpart deleted
  
  #import(foreach,"getDoParWorkers","getDoParRegistered","getDoParName","getDoParVersion","foreach","registerDoSEQ","setDoSeq","setDoPar")
  #importFrom(parallel, "makeCluster", "stopCluster", "detectCores", "clusterExport", "clusterEvalQ")
#' @importFrom doParallel registerDoParallel
#' @importFrom iterators  icount
#' @importFrom kSamples ad.test

#' @importFrom foreach %dopar% getDoParWorkers getDoParRegistered getDoParName getDoParVersion foreach registerDoSEQ setDoSeq setDoPar
#' @import parallel
# @importFrom grDevices adjustcolor colorRampPalette dev.cur dev.interactive dev.new devAskNewPage gray heat.colors palette
#' @importFrom graphics abline boxplot contour curve filled.contour image legend lines matlines pairs par persp plot points polygon rect rug stars text title
# @importFrom methods callGeneric
# @importFrom utils combn installed.packages modifyList setTxtProgressBar txtProgressBar globalVariables
  
  

  .onAttach <- function(lib, pkg,...){
    pkg.info <- drop(read.dcf(file=system.file("DESCRIPTION", package="fda.usc.devel"),
                              fields=c("Title","Version","Date")))
    foo <- suppressWarnings(foreach::"%dopar%"(foreach::foreach(i=1), {}))
     packageStartupMessage(
      #paste("----------------------------------------------------------------------------------\n",pkg.info["Title"]),"\n",
      #	 "Functional Data Analysis in R \n",
      #paste(" fda.usc version ", pkg.info["Version"]," (built on ", pkg.info["Date"], ") is now loaded\n", sep=""),
paste("fda.usc.devel is running sequentially usign foreach package\n"),
paste("Please, execute ops.fda.usc() once to run in local parallel mode\n"),
paste("Deprecated functions (v-2.0.0): min.basis, min.np, anova.hetero,\nanova.onefactor and anova.RPm are renamed by: optim.basis, optim.np, \nfanova.hetero, fanova.onefactor, fanova.RPm\n"), 
paste("New functions (v-2.1.1): fregre.mlm.fr, fregre.lm.fr, fregre.sam.fr\nand fregre.kam.fr\n"),      
            #" ops.fda.usc() changes the parameters of the package\n",
      "----------------------------------------------------------------------------------\n"
    )
  }
  
# fda.usc 2.0.02019-11-14