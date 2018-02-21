# GuanRank

This is the R version of GuanRank.

Example R code:

source("guanrank.R")
tbl=data.frame(time=c(430,257,185,298,506),status=c(0,1,0,1,0))
rownames(tbl)=paste0("patient",1:5)
tbl

|          | time | status |
| -------- | ---- | ------ |
| patient1 |  430 |      0 |
| patient2 |  257 |      1 |
| patient3 |  185 |      0 |
| patient4 |  298 |      1 |
| patient5 |  506 |      0 |

guanrank(tbl)

|          | time | status |      rank | 
| -------- | ---- | ------ | --------- |
| patient3 |  185 |      0 | 0.6000000 |
| patient2 |  257 |      1 | 1.0000000 |
| patient4 |  298 |      1 | 0.6666667 |
| patient1 |  430 |      0 | 0.2000000 |
| patient5 |  506 |      0 | 0.2000000 |

guanrank_elf(tbl) # new version

|          | time | status |      rank |
| -------- | ---- | ------ | --------- |
| patient3 |  185 |      0 | 0.4588353 |
| patient2 |  257 |      1 | 1.0000000 |
| patient4 |  298 |      1 | 0.8353414 |
| patient1 |  430 |      0 | 0.0000000 |
| patient5 |  506 |      0 | 0.0000000 |


To do:
add document about these functions


