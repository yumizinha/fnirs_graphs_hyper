##### ANALISE DE Modularidade
library(igraph)
library(tidyverse)
library(fdrtool)
source("./hyper_plot.R")
Z = NULL


###############################################
# Processamento dos dados, selecao de canais
# Cálculo do sinal com short distance
##############################################

# definição do total de duplas utilizadas no experimento
todos = 25


# Separar apenas o experimento dos demais dados.
sep_blocks <- function(subj_triggers, type)
{
  j = subj_triggers$Subject
  fname = paste("./Data/", type, j, "_001_oxyhb.rds", sep="")
  data = read_rds(fname)
  aux = data[subj_triggers$Presenting_ini:nrow(subj_triggers),]
  fsave = paste("./Data/", type, j, "_001_oxyhb_aula", sep = "")
  
  
  # Codigo estava no for de bloco, que nao existe mais
  fname = paste(fsave,".txt", sep="")
  write.table(aux, fname)
    
}


# Processamento dos triggers e divisão dos dados em quatro arquivos
triggers = read_csv2("./Data/triggers_frames_aula.csv")
for(ii in 1:nrow(triggers))
{
  subj_triggers = triggers[ii, ]
  sep_blocks(subj_triggers, "prof")
  sep_blocks(subj_triggers, "aluno")
}




