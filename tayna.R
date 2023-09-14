##### ANALISE DE Modularidade
library(igraph)
library(tidyverse)
Z = NULL


########################################################################
# Leitura dos dados conforme bloco a ser analisado
# A variavel bloco determina qual arquivo de dados me interessa, sendo:
# bloco 1 - baseline. Coleta de 1min de repouso - NAO ENVIADO
# bloco 2 - Atividades independentes da dupla - NAO ENVIADO
# bloco 3 - Primeira parte da aula: predominância em expositiva
# bloco 4 - Continuação da aula: predominância em interativa
########################################################################

bloco = 4
#bloco = 3
todos = c(1:12,14:25) #duplas que serao analisadas. Nesse caso, a dupla 13 nao atende requisitos.


# Leitura de Dados do txt conforme o bloco selecionado
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
  # todos os canais contra todos os canais e zerando a partir de um limiar (cuttoff)
  # correlação negativa não é considerada por questões de ativação neuronal:
  # um aluno em repouso e professor falando pode ter correlação altamente negativa
  # mas não é uma sincronização do tipo que estamos interessados.
  
  CONNECT = cor(cbind(prof, aluno))   #,method="spearman")
  cutoff = 0.25
  # para imprimir cada matriz:
  # write.csv(CONNECT, 'MATRIZ5.csv')
  
  # Simplifiquei zerando as conexões mais fracas. 
  # Porém, não necessariamente é importante fazer isso
  CONNECT[which(CONNECT < cutoff)] = 0
  SUBJECT = c(rep(1, ncol(prof)), rep(2, ncol(aluno)))
  
  
  #Analise de Grafos a partir da modularidade (zorig)
  g = graph.adjacency(CONNECT, weighted=TRUE, mode="undirected", diag=FALSE)
  Zorig = modularity(g, SUBJECT)
  
 
}

  
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
  
# molde de duas cabecas - desenho em que os nós dos grafos ficam
molde <- readPNG("molde_hypper.png")
  
# 1) Essa parte do código é somente para posicionar os grafos conforme as coordenadas espaciais
# 1a) utilizei uma biblioteca em python para achar a posição no plano cartesiano
# (calculando ajuste das coordenadas calculadas no eeg_positions)
# 1b) fiz os deslocamentos nos eixos x e y para que as posições coincidam com a máscara
# (deslocando cérebro do aluno para separar as cabeças corretamente)

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

par(mfrow =c(25, 2))
todos = c(1:12,14:25)
# cria um pdf em que poderá visusalizar os grafos gerados
pdf("centrality.pdf", width = 17, height = 7)


# Alterar bloco 3 e 4 para as funções abaixo gerarem a matriz de centralidades

# bloco =3
bloco = 4

# Caso vá estudar as centralidades dos grafos, criei essas matrizes.
# Spoiler alert: Não deu em nada ainda.
#matriz_centralidade_3 = matrix(NA, nrow = max(todos), ncol = 32)
matriz_centralidade_4 = matrix(NA, nrow = max(todos), ncol = 32)


for(j in todos){ #for da dupla
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
  CONNECT[which(CONNECT < cutoff)] = 0
  SUBJECT = c(rep(1, ncol(prof)), rep(2, ncol(aluno)))
  
  # zerando conexoes intracerebrais
  CONNECT[1:16, 1:16] = 0
  CONNECT[17:32, 17:32] = 0
  
  # grafo nao direcionado a partir da matriz de adjacencia e calculo centralidade
  g = graph.adjacency(CONNECT, weighted = TRUE, mode = "undirected", diag = FALSE)
  centrality = eigen_centrality(g, directed = FALSE, scale = TRUE,
                                weights = NULL, options = arpack_defaults)

  # Centralidades para outra análise
  #  matriz_centralidade_2[j,]= centrality$vector
  #  matriz_centralidade_3[j,]= centrality$vector
  matriz_centralidade_4[j,]= centrality$vector
  
  
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
  
  
}
dev.off()


# Codigo para escrever as centralidades dos grafos
#write_csv(as.data.frame(matriz_centralidade_3), "introducao.csv")
#write_csv(as.data.frame(matriz_centralidade_4), "interativo.csv")
#write_csv(as.data.frame(matriz_centralidade_2), "independent.csv")
