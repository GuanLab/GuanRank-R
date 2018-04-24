# cindex.r

## 1. The name of data&status need to be changed! (D_PFS, D_PFS_FLAG)
calculate.concordanceIndex <- function(predicted, D_PFS, D_PFS_FLAG) {
  suppressPackageStartupMessages(library("survcomp"))
  cIndex <- survcomp::concordance.index(x = predicted, surv.time = D_PFS, surv.event  = D_PFS_FLAG)
}

