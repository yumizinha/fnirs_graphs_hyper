library(glasso)
library(tidyverse)

# definição do total de duplas utilizadas no experimento
todos = 25
j = 1

prof = read.table(paste("./Data/", "prof", j, "_001_oxyhb_bloco", bloco, ".txt", sep=""))
aluno = read.table(paste("./Data/", "aluno", j, "_001_oxyhb_bloco", bloco, ".txt", sep=""))
GSR = rowMeans(prof)
for(i in 1:ncol(prof)){prof[, i] = lm(prof[, i]~GSR)$resid}
GSR = rowMeans(aluno)
for(i in 1:ncol(aluno)){aluno[, i ] = lm(aluno[, i]~GSR)$resid}
cov_chan = cov(cbind(prof, aluno))
for(ii in 1:32)
{
  for(jj in 1:32)
  {
    aux_ii = ii <= 16
    aux_jj = jj <= 16
    if(ii != jj) cov_chan[ii, jj] = cov_chan[ii, jj]*abs(aux_ii - aux_jj)
  }
}

d = dim(cov_chan)[2]
rho = 0.00001
heatmap(glasso(cov_chan, rho)$w)

