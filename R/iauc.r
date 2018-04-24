# iauc.r

## 1. The evaluation times need to be modified!
## 2. The name of data&status need to be changed! (D_PFS, D_PFS_FLAG)
calculate.integratedAUC <- function(predicted, D_PFS, D_PFS_FLAG, times = 30.5 * c(14, 16, 18, 20, 22)) {
  suppressPackageStartupMessages(library("timeROC"))
  suppressPackageStartupMessages(library("risksetROC"))
  tempAUC <- timeROC(T=D_PFS, delta=D_PFS_FLAG, marker=predicted,cause=1, times=times)
  iaucs <- IntegrateAUC(tempAUC$AUC, tempAUC$times, tempAUC$survProb,tmax = max(tempAUC$times))
  return(iaucs)
}
