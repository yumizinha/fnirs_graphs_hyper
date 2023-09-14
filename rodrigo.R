##### ANALISE DE Modularidade
library(igraph)
library(tidyverse)
Z = NULL

###############################################
# Processamento dos dados, selecao de canais
# Cálculo do sinal com short distance
##############################################

# definição do total de duplas utilizadas no experimento
todos = 25

# Processamento de short distance
# A montagem utilizada possui os seguintes canais short distance:
short = c(3, 6, 9, 12, 16, 19, 22, 24)
short_n = as.character(short)
short_n = paste(rep("V", length(short)), short_n, sep="")

# os demais canais (long) são os que possuem sinais a serem utilizados na análise
long = (1:24)[-short]
long_n = as.character(long)
long_n = paste(rep("V", length(long)), long_n, sep="")

# Processamento sinal integral. Resíduos da série inteira considerando short distance
remove_short <- function(data)
{
  for(ii in long)
  {
    aux = paste(short_n, collapse = "+")
    formula = as.formula(paste("V", ii, "~", aux, sep=""))
    data[,ii] <- lm(formula, data)$residuals
  }
  data
}

# separa os dados nos blocos do experimento, conforme triggers pre-definidos e salvo em arquivo
# sao quatro blocos: baseline (inicio do código), atividade independente, introdução e atividade
sep_blocks <- function(subj_triggers, type)
{
  j = subj_triggers$Subject
  fname = paste("./Data/", type, j, "_001_oxyhb.rds", sep="")
  data = read_rds(fname)
  aux = list(
    data[subj_triggers$Baseline_ini:subj_triggers$Baseline_fim,],
    data[subj_triggers$Independent_ini:subj_triggers$Independent_fim,],
    data[subj_triggers$Presenting_ini:subj_triggers$Presenting_fim,],
    data[subj_triggers$Activities_ini:subj_triggers$Activities_fim,])
  fsave = paste("./Data/", type, j, "_001_oxyhb_bloco", sep = "")
  for(bloco in 1:4) {
    fname =  paste(fsave, bloco, ".txt", sep="")
    write.table(aux[[bloco]], fname)
  }
}

for(j in 1:todos) { #for da dupla
  #Leitura de Dados
  fname = paste("./Data/", "prof", j, "_001_oxyhb.rds", sep="")
  paste("./Data/", "prof", j, "_001_oxyhb.txt", sep="") %>% 
    read.table() %>% 
    as_tibble() %>% 
    remove_short() %>%
    select(long_n) %>% 
    write_rds(fname)
  
  fname = paste("./Data/", "aluno", j, "_001_oxyhb.rds", sep="")
  paste("./Data/", "aluno", j, "_001_oxyhb.txt", sep="") %>% 
    read.table() %>% 
    as_tibble() %>% 
    remove_short() %>%
    select(long_n) %>% 
    write_rds(fname)
}

# Processamento dos triggers e divisão dos dados em quatro arquivos
triggers = read_csv2("./Data/triggers_frames.csv")
for(ii in 1:nrow(triggers))
{
  subj_triggers = triggers[ii, ]
  sep_blocks(subj_triggers, "prof")
  sep_blocks(subj_triggers, "aluno")
}



########################################################################
# Leitura dos dados conforme bloco a ser analisado
# A variavel bloco determina qual arquivo de dados me interessa, sendo:
# bloco 1 - baseline. Coleta de 1min de repouso
# bloco 2 - Atividades independentes da dupla
# bloco 3 - Primeira parte da aula: predominância em expositiva
# bloco 4 - Continuação da aula: predominância em interativa
########################################################################


bloco = 3
todos = c(1:12,14:25) #duplas que serao analisadas. Nesse caso, a dupla 13 nao atende requisitos.


# Leitura de Dados
for(j in todos){ #for da dupla
  prof = read.table(paste("./Data/", "prof", j, "_001_oxyhb_bloco", bloco, ".txt", sep=""))
  aluno = read.table(paste("./Data/", "aluno", j, "_001_oxyhb_bloco", bloco, ".txt", sep=""))
  
  #GSR - pega os resíduos do sinal médio de cada canal dos professores e alunos
  GSR = rowMeans(prof)
  for(i in 1:ncol(prof)){prof[, i] = lm(prof[, i]~GSR)$resid}
  
  GSR = rowMeans(aluno)
  for(i in 1:ncol(aluno)){aluno[, i ] = lm(aluno[, i]~GSR)$resid}
  
  
  # Preprocessamento
  # Matriz de conectividade funcional - encontra a correlação dos canais
  # conexão entre cada canal
  # matriz 32x32: correlação dos 16 canais de professor com 16 canais de alunos.
  # todos os canais contra todos os canais e zerando corr<0.1
  
  CONNECT = cor(cbind(prof, aluno))   #,method="spearman")
  #para imprimir cada matriz write.csv(CONNECT, 'MATRIZ5.csv')
  
  CONNECT[which(CONNECT < 0.2)] = 0
  SUBJECT = c(rep(1, ncol(prof)), rep(2, ncol(aluno)))
  
  
  #Analise de Grafos a partir da modularidade (zorig)
  g = graph.adjacency(CONNECT, weighted=TRUE, mode="undirected", diag=FALSE)
  Zorig = modularity(g, SUBJECT)
  
  
  #################################################
  # ANALISE CENTRALIDADE DIVISAO EM DOIS CEREBROS #
  # PLOT DE GRAFOS COM TEMPLATE DAS CABEÇAS       #
  #################################################
  
  require(igraph)
  require(caret)
  require(tidyverse)
  require(png)
  require(lmtest)
  require(xts)
  
  # molde de duas cabecas
  molde <- readPNG("molde_hypper.png")
  
  
  # Calculando ajuste das coordenadas calculadas no eeg_positions
  # deslocando cérebro do aluno para separar as cabeças corretamente
  # Alterei o deslocamento para coincidir com o template do cérebro carregado
  coords_fnirs_prof = read.csv("coords_fNIRS.csv") %>%
    mutate(x= x-1.52) %>%
    mutate(y=y-0.028)
  
  coords_fnirs_aluno = read.csv("coords_fNIRS.csv") %>% 
    mutate(x = x+1.53) %>%
    mutate(y=y-0.028)
  
  # unindo as novas coordenadas para plot do grafo
  coords_fnirs = rbind(coords_fnirs_prof, coords_fnirs_aluno) %>% 
    select(x, y) %>% 
    as.matrix()
  
  
  # Leitura dos dados de cada uma das duplas. 
  # Os BDs indicados por _P são professores
  # BDs indicados sem _P são os alunos. Somente oxyhb foi indicado.
  
  # se tiver c(1, 1) é porque gerei da última vez um grafo de cada vez.
  # Ver dentro do for se está fixando uma única dupla e deletar se for gerar todos
  
  par(mfrow =c(1, 1))
  bloco =3
  
  #for(j in 1:todos){ #for da dupla
  
    j = 5
      #Leitura de Dados
    
    prof = read.table(paste("./Data/", "prof", j, "_001_oxyhb_bloco", bloco, ".txt", sep=""))
    aluno = read.table(paste("./Data/", "aluno", j, "_001_oxyhb_bloco", bloco, ".txt", sep=""))
    
    #GSR - pega os resíduos do sinal médio de cada canal dos professores e alunos
    GSR = rowMeans(prof)
    for(i in 1:ncol(prof)){prof[, i] = lm(prof[, i]~GSR)$resid}
    
    GSR = rowMeans(aluno)
    for(i in 1:ncol(aluno)){aluno[, i ] = lm(aluno[, i]~GSR)$resid}
    
    # Preprocessamento
    # Matriz de conectividade funcional - encontra a correlação dos canais
    # conexão entre cada canal
    # matriz 36x36: correlação dos 16 canais de professor com 16 canais de alunos.
    # todos os canais contra todos os canais e zerando corr<0.1
    
    CONNECT = cor(cbind(prof, aluno))   #,method="spearman")
    CONNECT[which(CONNECT < 0.2)] = 0
    SUBJECT = c(rep(1, ncol(prof)), rep(2, ncol(aluno)))
    
    # zerando conexoes intracerebrais
    CONNECT[1:16, 1:16] = 0
    CONNECT[17:32, 17:32] = 0
    
    # grafo nao direcionado a partir da matriz de adjacencia e calculo centralidade
    g = graph.adjacency(CONNECT, weighted = TRUE, mode = "undirected", diag = FALSE)
    centrality = eigen_centrality(g, directed = FALSE, scale = TRUE,
                                  weights = NULL, options = arpack_defaults)
    
    
    
    # configurando paleta de cor por centralidade
    fine = 500 # this will adjust the resolving power.
    pal = colorRampPalette(c('red','green'))
    #this gives you the colors you want for every point
    graphCol = pal(fine)[as.numeric(cut(centrality$vector,breaks = fine))]
    
    
    # colocar centrality$vector como variavel em vertex.color
    plot(g, vertex.color=graphCol, layout = coords_fnirs, vertex.size = 10, rescale = FALSE)
    
    
    # Tentando colocar o pano de fundo do hypper
    lim <- par()
    rasterImage(molde,
                xleft=-2.6, xright=2.6, 
                ybottom=-1.1, ytop=1.1)
    
    
  #}
  