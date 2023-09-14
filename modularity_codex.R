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



# pedaço de código que ocultei por não precisar mais pegar a média do sinal. O intercepto ao fazer
#lm e resíduos de short distance já resolve o alinhamento dos sinais.
#for(j in 1:todos){ #for da dupla
  #Leitura de Dados
 # prof = read.table(paste("./Data/", "prof", j, "_001_oxyhb.txt", sep=""))
 # aluno = read.table(paste("./Data/", "aluno", j, "_001_oxyhb.txt", sep=""))
  
  #GSR - pega os resíduos do sinal médio de cada canal dos professores e alunos
  #GSR = rowMeans(prof)
  #for(i in 1:ncol(prof)){prof[, i] = lm(prof[, i]~GSR)$resid}
  
  #GSR = rowMeans(aluno)
  #for(i in 1:ncol(aluno)){aluno[, i ] = lm(aluno[, i]~GSR)$resid}




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
  
  CONNECT[which(CONNECT < 0.25)] = 0
  SUBJECT = c(rep(1, ncol(prof)), rep(2, ncol(aluno)))
  
  
  #Analise de Grafos a partir da modularidade (zorig)
  g = graph.adjacency(CONNECT, weighted=TRUE, mode="undirected", diag=FALSE)
  Zorig = modularity(g, SUBJECT)
  
  
  #Bootstrap
  Zboot = NULL
  NBOOT = 1000
  for(boot in 1:NBOOT){
    aux = sample(nrow(prof)-1, 1)
    auxsign = sample(c(-1, 1), 1)
    BOOTprof = prof[c((aux+1):nrow(prof), 1:aux),]*auxsign
    
    #Matriz de conectividade funcional
    CONNECT = cor(cbind(BOOTprof, aluno))  #,method="spearman")
    CONNECT[which(CONNECT < 0.25)] = 0
    SUBJECT = c(rep(1, ncol(prof)), rep(2, ncol(aluno)))
    g = graph.adjacency(CONNECT, weighted = TRUE, mode = "undirected", diag = FALSE)
    Zboot = c(Zboot, modularity(g, SUBJECT))
  }
  pvalue = length(which(Zboot<=Zorig))/NBOOT
  print(c(j, Zorig, pvalue))
}


################################################
# Descritiva
###############################################



##### Analise de hubs intra e intercerebro (bridge)
require(igraph)
DegreePROF = matrix(0, length(todos), 16)
DegreeALUNO = matrix(0, length(todos), 16)
INTERCEREBRO = matrix(0, length(todos), 32)
par(mfrow = c(2, 2))

for(j in 1:todos){
  print(j)
  #PROCESSAMENTO IDENTICO AO ANTERIOR
  #Leitura de Dados
  prof = read.table(paste("./Data/", "prof", j, "_001_oxyhb_bloco", bloco, ".txt", sep=""))
  aluno = read.table(paste("./Data/", "aluno", j, "_001_oxyhb_bloco", bloco, ".txt", sep=""))
  
  #GSR
  GSR = rowMeans(prof)
  for(i in 1:ncol(prof)){prof[, i] = lm(prof[, i]~GSR)$resid}
  GSR = rowMeans(aluno)
  for(i in 1:ncol(aluno)){aluno[, i] = lm(aluno[, i]~GSR)$resid}
  
  CONNECT = cor(cbind(prof, aluno))   #,method="spearman")
  CONNECT[which(CONNECT < 0.1)] = 0
  
  #grafos intracerebro
  matrizPROF = CONNECT[1:16, 1:16]
  diag(matrizPROF) = 0 #tira a diagonal principal
  matrizALUNO = CONNECT[17:32, 17:32]
  diag(matrizALUNO) = 0
  
  #calculo do degree intracerebro
  # VERIFICAR PQ colMeans ao invés de colSums - corrigi de acordo com link
  DegreePROF[j,] = colSums(matrizPROF)
  DegreeALUNO[j,] = colSums(matrizALUNO)
  
  #Analise do intercerebro
  #Zera todas as conexoes intracerebro
  CONNECT[1:16, 1:16] = 0
  CONNECT[17:32, 17:32] = 0
  SUBJECT = c(rep(1, ncol(prof)), rep(2, ncol(aluno)))
  g = graph.adjacency(CONNECT, weighted=TRUE, mode="undirected", diag=FALSE)
  plot(g, vertex.color = SUBJECT)
  Zorig = modularity(g, SUBJECT)
  INTERCEREBRO[j,] = colSums(CONNECT)
}


par(mfrow =c(2, 2))
#Descritivas Degree PROFESSOR
Z = DegreePROF
MEDIA = array(0, ncol(Z))
SD = array(0, ncol(Z))

for(i in 1:ncol(Z)){
  MEDIA[i] = mean(Z[, i])
  SD[i] = sd(Z[, i])
}

COHENd = MEDIA/SD
names(COHENd) = colnames(prof)
barplot(COHENd, main = "Hub Teacher", xlab = "Channel", ylab = "Cohen-D")

#Descritivas Degree ALUNO
Z = DegreeALUNO
MEDIA = array(0, ncol(Z))
SD = array(0, ncol(Z))

for(i in 1:ncol(Z)){
  MEDIA[i] = mean(Z[, i])
  SD[i] = sd(Z[, i])
}
COHENd = MEDIA/SD
names(COHENd) = colnames(aluno)
barplot(COHENd, main = "Hub Student", xlab = "Channel", ylab = "Cohen-D")


#Descritivas Degre INTERCEREBRO - Bridges
Z = INTERCEREBRO
MEDIA = array(0, ncol(Z))
SD = array(0, ncol(Z))
for(i in 1:ncol(Z)){
  MEDIA[i] = mean(Z[, i])
  SD[i] = sd(Z[, i])
}
COHENd = MEDIA/SD
names(COHENd) = colnames(CONNECT)
barplot(COHENd[1:16], main="Bridges Teacher", xlab="Channel", ylab="Cohen-D")
barplot(COHENd[17:32], main="Bridges Student", xlab="Channel", ylab="Cohen-D")



###############################################################
# Análise de clustering Espectral
# tentativa de construir os grafos com o clustering espectral 
# usando a matriz de correlacao como similaridade
###############################################################


require(stats)

# cria matriz de afinidade para clustering espectral
make.affinity <- function(S, n.neighboors=2) {
  N <- length(S[,1])
  
  if (n.neighboors >= N) {  # fully connected
    A <- S
  } else {
    A <- matrix(rep(0,N^2), ncol=N)
    for(i in 1:N) { # for each line
      # only connect to those points with larger similarity 
      best.similarities <- sort(S[i,], decreasing=TRUE)[1:n.neighboors]
      for (s in best.similarities) {
        j <- which(S[i,] == s)
        A[i,j] <- S[i,j]
        A[j,i] <- S[i,j] # to make an undirected graph, ie, the matrix becomes symmetric
      }
    }
  }
  A  
}

par(mfrow = c(2, 2))
todos = c(1:12,14:25)
pdf("spectral_clustering.pdf")

for(j in 1:todos) { 

  print(j)
  #PROCESSAMENTO IDENTICO AO ANTERIOR
  #Leitura de Dados
  prof = read.table(paste("./Data/", "prof", j, "_001_oxyhb_bloco", bloco, ".txt", sep=""))
  aluno = read.table(paste("./Data/", "aluno", j, "_001_oxyhb_bloco", bloco, ".txt", sep=""))
  
  #GSR
  GSR = rowMeans(prof)
  for(i in 1:ncol(prof)){prof[, i] = lm(prof[, i]~GSR)$resid}
  GSR = rowMeans(aluno)
  for(i in 1:ncol(aluno)){aluno[, i] = lm(aluno[, i]~GSR)$resid}
  
  data = cbind(prof, aluno)
  S_orig = cor(data)
  S = S_orig
  #S[1:18, 1:18] = 0
  #S[19:36, 19:36] = 0
  S = S + diag(1, 32)
  A <- make.affinity(S, length(todos))
  D <- diag(apply(A, 1, sum))
  U <- D - A
  k   <- 3
  evL <- eigen(U, symmetric=TRUE)
  Z   <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]
  cc = 3
  km <- kmeans(Z, centers = cc, nstart=5)
  cl = km$cluster
  g = graph.adjacency(S, weighted=TRUE, mode="undirected", diag=FALSE)
  modularity(g, cl)
  
  MM = matrix(NA, nrow = 32, ncol = 32)
  for(aa in 1:32)
  {
    for(bb in 1:32)
    {
      MM[aa, bb] <- (cl[aa] == cl[bb])
    }
  }
  colnames(MM) = rep(1:16, 2)
  rownames(MM) = rep(1:16, 2)
  g = graph.adjacency(MM, weighted=TRUE, mode="undirected", diag=FALSE)
  SUBJECT = c(rep(1, ncol(prof)), rep(2, ncol(aluno)))
  plot(g, vertex.color = SUBJECT)
}
dev.off()


#######################################################################################
# BOOTSTRAP RETIRANDO 1 CANAL DE CADA VEZ
# Leitura dos dados de cada uma das duplas. Os BDs indicados por _P são professores
# BDs indicados sem _P são os alunos. Somente oxyhb foi indicado (VER JOAO).
#######################################################################################

data = as.list(rep(NA, 16))
for(y in 1:16) {
  print(y)
  aux_da_yuyu = matrix(NA, nrow = length(todos), ncol = 3)
  for(j in todos){ #for da dupla
    print(j)
    #Leitura de Dados
    prof = read.table(paste("./Data/", "prof", j, "_001_oxyhb_bloco", bloco, ".txt", sep=""))
    prof = prof[,-y]
    aluno = read.table(paste("./Data/", "aluno", j, "_001_oxyhb_bloco", bloco, ".txt", sep=""))
    aluno = aluno[,-y]
    
    #GSR - pega os resíduos do sinal médio de cada canal dos professores e alunos
    GSR = rowMeans(prof)
    for(i in 1:ncol(prof)){prof[, i] = lm(prof[, i]~GSR)$resid}
    
    GSR = rowMeans(aluno)
    for(i in 1:ncol(aluno)){aluno[, i ] = lm(aluno[, i]~GSR)$resid}
    
    # Preprocessamento
    # Matriz de conectividade funcional - encontra a correlação dos canais
    # conexão entre cada canal
    # matriz 36x36: 18 canais de professor com 18canais de alunos.
    # todos os canais contra todos os canais e zerando corr<0.1
    
    CONNECT = cor(cbind(prof, aluno))   #,method="spearman")
    CONNECT[which(CONNECT < 0.1)] = 0
    SUBJECT = c(rep(1, ncol(prof)), rep(2, ncol(aluno)))
    
    
    #Analise de Grafos
    g = graph.adjacency(CONNECT, weighted=TRUE, mode="undirected", diag=FALSE)
    Zorig = modularity(g, SUBJECT)
    
    
    #Bootstrap
    Zboot = NULL
    NBOOT = 10000
    for(boot in 1:NBOOT){
      aux = sample(nrow(prof)-1, 1)
      auxsign = sample(c(-1, 1), 1)
      BOOTprof = prof[c((aux+1):nrow(prof), 1:aux),]*auxsign
      
      #Matriz de conectividade funcional
      CONNECT = cor(cbind(BOOTprof, aluno))  #,method="spearman")
      CONNECT[which(CONNECT < 0.1)] = 0
      SUBJECT = c(rep(1, ncol(prof)), rep(2, ncol(aluno)))
      g = graph.adjacency(CONNECT, weighted = TRUE, mode = "undirected", diag = FALSE)
      Zboot = c(Zboot, modularity(g, SUBJECT))
    }
    pvalue = length(which(Zboot<=Zorig))/NBOOT
    aux_da_yuyu[j, ] = c(j, Zorig, pvalue)
    print(aux_da_yuyu)
  }
  data[[y]] = aux_da_yuyu
  print(paste("Terminei o canal!", y))
  print(data[[y]])
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

# molde de duas cabecas desenhado para plot do grafo
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

# se estiver par(mfrow =c(1, 1)) é porque gerei da última vez um grafo de cada vez.
# Ver dentro do for se está fixando uma única dupla e deletar se for gerar todos
#par(mfrow =c(1, 1))

par(mfrow =c(25, 2))
todos = c(1:12,14:25)
pdf("centrality.pdf", width = 17, height = 7)


# Alterar bloco 3 e 4 para as funções abaixo gerarem a matriz de centralidades
#bloco = 4
 bloco =3
#bloco = 2
#matriz_centralidade_2 = matrix(NA, nrow = max(todos), ncol = 32)
matriz_centralidade_3 = matrix(NA, nrow = max(todos), ncol = 32)
#matriz_centralidade_4 = matrix(NA, nrow = max(todos), ncol = 32)


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
  CONNECT[which(CONNECT < 0.0)] = 0
  SUBJECT = c(rep(1, ncol(prof)), rep(2, ncol(aluno)))
  
  # zerando conexoes intracerebrais
  CONNECT[1:16, 1:16] = 0
  CONNECT[17:32, 17:32] = 0
  
  # grafo nao direcionado a partir da matriz de adjacencia e calculo centralidade
  g = graph.adjacency(CONNECT, weighted = TRUE, mode = "undirected", diag = FALSE)
  centrality = eigen_centrality(g, directed = FALSE, scale = TRUE,
                                weights = NULL, options = arpack_defaults)
#  matriz_centralidade_2[j,]= centrality$vector
  matriz_centralidade_3[j,]= centrality$vector
#  matriz_centralidade_4[j,]= centrality$vector
  
  
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

write_csv(as.data.frame(matriz_centralidade_3), "introducao.csv")
write_csv(as.data.frame(matriz_centralidade_4), "interativo.csv")
write_csv(as.data.frame(matriz_centralidade_2), "independent.csv")


####################################################################################################
# Teste T e Mann Whitney Wilcoxon para comparação de medidas de centralidade por autovetor 
# Por definição, a centralidade dos grafos não segue distribuição normal,
# mas sendo rebelde aqui para verificação e comparar ambos os testes
####################################################################################################

introducao = read_csv("introducao.csv")
interativo = read_csv("interativo.csv")


# Fazendo teste T mesmo assim:

# medidas repetidas em mesmo indivíduo (comparando duas condições):
# considerando as diferenças entre as condições:
canais = paste(c(rep("P_", 16), rep("A_", 16)), colnames(CONNECT), sep = "")

aux = tibble(na.omit(interativo - introducao))
colnames(aux) = canais

# matriz para leitura dos canais e identificação, já que o csv da centralidade foi na ordem
# e não considera as localizações identificadas nos grafos:

pvalores <- matrix(NA, nrow = length(canais), ncol = 4)
rownames(pvalores) = canais
estatisticas <- matrix(NA, nrow = length(canais), ncol = 2)
rownames(estatisticas) = canais

#  teste dos canais, NÃO CONSIDERANDO BONFERRONI
for(ii in 1:length(canais))
{
  pvalores[ii, 1] <- t.test(aux[[canais[ii]]])$p.value
  pvalores[ii, 2] <- wilcox.test(aux[[canais[ii]]])$p.value
  pvalores[ii, 3] <- wilcox.test(aux[[canais[ii]]], alternative = "less")$p.value
  pvalores[ii, 4] <- wilcox.test(aux[[canais[ii]]], alternative = "greater")$p.value
  estatisticas[ii, 1] <- t.test(aux[[canais[ii]]])$statistic
  estatisticas[ii, 2] <- wilcox.test(aux[[canais[ii]]])$statistic
}
fdrtool(pvalores[,3], statistic = "pvalue")
# pvalores[,1] # teste t
round(pvalores[,2], 3) # Mann-Whitney
which(round(pvalores[,2], 2) <= 0.05)

hyper_plot(as.numeric(pvalores[,2] < 0.05/32))
hyper_plot(as.numeric(pvalores[,3] < 0.05/32)) # introducao > interativo
hyper_plot(as.numeric(pvalores[,4] < 0.05/32)) # introducao < interativo


# 1) teste feito com bonferrroni e só resultou em um canal. Testei sem cutoff nas arestas tb
# Resp: Não funcionou
# 1b) Ao invés de Bonferroni fazer FDR
# Nem implementei. FDR é recomendado para dados independentes e não é isso que ocorre
# 2) fazer centralidade de grau e não eigenvector nos grafos totais (sem cuttoff)
# troquei e também não teve resultado significativo
# 2b) Agrupar em áreas de interesse: máximas/média de centralidade (PF esq, PF dir, TPJ) de cada um.
# 3) Focar resultados na modularidade das coisas.






###########################################################
# USANDO CAUSALIDADE DE GRANGER COMO MATRIZ DE ADJACENCIA #
###########################################################


# Ajustamos dois modelos para cada uma das "vias": prever série do professor em função
# do aluno + prever professor só com lag da série do professor 
# e prever o aluno em função do professor + prever aluno em função apenas do lag do aluno.
# Com o resultado, pegaria a razão da variância dos resídios e inverteria, para ter uma matriz
# de adjacência com os índices. Ainda, zeraria os casos em que o teste de Granger mostrariam
# que o modelo não é significativo.
# Resultado: Não deu certo para ajuste dos modelos.



for(j in 1:todos){ #for da dupla
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
  # matriz 32x32: correlação dos 16 canais de professor com 16 canais de alunos.
  # todos os canais contra todos os canais e zerando corr<0.1
  
  # Resíduos dos modelos com as séries temporais
  CONNECT_g = matrix(0, nrow = 32, ncol = 32)
  
  
  dt = 10
  # criar variável t-1 até t-10
  # pra cada equação vou ter 10 do passado x e 10 do passado y
  for (k in 1:16) {
    for (l in 1:16) {
      y = xts(aluno[,l], as.Date(1:length(aluno[,l])))
      x = xts(prof[, k], as.Date(1:length(prof[, k])))
      
      y1 = lm(y[-(1:dt),] ~ lag(x, k=1:dt)[-(1:dt),] + lag(y, k=1:dt)[-(1:dt),])$resid
      y2 = lm(y[-(1:dt),] ~ lag(y, k=1:dt)[-(1:dt),])$resid
      x1 = lm(x[-(1:dt),] ~ lag(x, k=1:dt)[-(1:dt),] + lag(y, k=1:dt)[-(1:dt),])$resid
      x2 = lm(x[-(1:dt),] ~ lag(x, k=1:dt)[-(1:dt),])$resid
      
      CONNECT_g[k,l+16] = var(y1)/var(y2)
      CONNECT_g[l+16,k] = var(x1)/var(x2)
    }
  }
  
  
  
  
  # matriz pra pegar o p-valor das séries do teste Granger
  CONNECT_g = matrix(0, nrow = 32, ncol = 32)
  # pra cada coluna de prof, calcular a causalidade de Granger pra mesma coluna de aluno:
  for (k in 1:16) {
    for (l in 1:16) {
      A = grangertest(prof[,k], aluno[,l], order = 1)
      B = grangertest(aluno[,l], prof[,k], order = 1)
      CONNECT_g[k,l+16] = A$"Pr(>F)"[2]
      CONNECT_g[l+16,k] = B$"Pr(>F)"[2]
      
    }
  }
  
  CONNECT_g < 0.005
  
}




####################################################################################################
# VERIFICACAO DE PARES TROCADOS. PROFESSORES COM OUTROS ALUNOS E PROFESSORES COM OUTROS PROFESSORES.#
# ALUNOS COM ALUNOS
####################################################################################################


##### ANALISE DE Modularidade
require(igraph)
Z = NULL



#### PROFESSOR COM CADA ALUNO, MESMO SEM INTERACAO COM ELES NO EXPERIMENTO
for(j in 1:todos){ #for da dupla
  #Leitura de Dados
  for (k in 1:todos) {
    prof = read.table(paste("./Data/", "prof", j, "_001_oxyhb_bloco", bloco, ".txt", sep=""))
    aluno = read.table(paste("./Data/", "aluno", k, "_001_oxyhb_bloco", bloco, ".txt", sep=""))
    
    m = min(nrow(prof), nrow(aluno))
    prof = prof[1:m,]
    aluno = aluno[1:m,]
    
    
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
    
    CONNECT[which(CONNECT < 0.15)] = 0
    SUBJECT = c(rep(1, ncol(prof)), rep(2, ncol(aluno)))
    
    
    #Analise de Grafos a partir da modularidade (zorig)
    g = graph.adjacency(CONNECT, weighted=TRUE, mode="undirected", diag=FALSE)
    Zorig = modularity(g, SUBJECT)
    
    
    #Bootstrap
    Zboot = NULL
    NBOOT = 1000
    for(boot in 1:NBOOT){
      aux = sample(nrow(prof)-1, 1)
      auxsign = sample(c(-1, 1), 1)
      BOOTprof = prof[c((aux+1):nrow(prof), 1:aux),]*auxsign
      
      #Matriz de conectividade funcional
      CONNECT = cor(cbind(BOOTprof, aluno))  #,method="spearman")
      CONNECT[which(CONNECT < 0.15)] = 0
      SUBJECT = c(rep(1, ncol(prof)), rep(2, ncol(aluno)))
      g = graph.adjacency(CONNECT, weighted = TRUE, mode = "undirected", diag = FALSE)
      Zboot = c(Zboot, modularity(g, SUBJECT))
    }
    pvalue = length(which(Zboot<=Zorig))/NBOOT
    print(c(j,k, Zorig, pvalue))
  }
}


##### PROFESSOR COM OUTROS PROFESSORES

#### PROFESSOR COM CADA PROFESSOR, MESMO SEM INTERACAO COM ELES NO EXPERIMENTO

for(j in 1:todos){ #for da dupla
  #Leitura de Dados
  for (k in 1:todos) {

    prof = read.table(paste("./Data/", "prof", j, "_001_oxyhb_bloco", bloco, ".txt", sep=""))
    prof2 = read.table(paste("./Data/", "prof", k, "_001_oxyhb_bloco", bloco, ".txt", sep=""))
    
    
        
    m = min(nrow(prof), nrow(prof2))
    prof = prof[1:m,]
    prof2 = prof2[1:m,]
    
    
    #GSR - pega os resíduos do sinal médio de cada canal dos professores e alunos
    GSR = rowMeans(prof)
    for(i in 1:ncol(prof)){prof[, i] = lm(prof[, i]~GSR)$resid}
    
    GSR = rowMeans(prof2)
    for(i in 1:ncol(prof2)){prof2[, i ] = lm(prof2[, i]~GSR)$resid}
    
    
    # Preprocessamento
    # Matriz de conectividade funcional - encontra a correlação dos canais
    # conexão entre cada canal
    # matriz 36x36: correlação dos 18 canais de professor com 18 canais de alunos.
    # todos os canais contra todos os canais e zerando corr<0.1
    
    CONNECT = cor(cbind(prof, prof2))   #,method="spearman")
    #para imprimir cada matriz write.csv(CONNECT, 'MATRIZ5.csv')
    
    CONNECT[which(CONNECT < 0.15)] = 0
    SUBJECT = c(rep(1, ncol(prof)), rep(2, ncol(prof2)))
    
    
    #Analise de Grafos a partir da modularidade (zorig)
    g = graph.adjacency(CONNECT, weighted=TRUE, mode="undirected", diag=FALSE)
    Zorig = modularity(g, SUBJECT)
    
    
    #Bootstrap
    Zboot = NULL
    NBOOT = 1000
    for(boot in 1:NBOOT){
      aux = sample(nrow(prof)-1, 1)
      auxsign = sample(c(-1, 1), 1)
      BOOTprof = prof[c((aux+1):nrow(prof), 1:aux),]*auxsign
      
      #Matriz de conectividade funcional
      CONNECT = cor(cbind(BOOTprof, prof2))  #,method="spearman")
      CONNECT[which(CONNECT < 0.15)] = 0
      SUBJECT = c(rep(1, ncol(prof)), rep(2, ncol(prof2)))
      g = graph.adjacency(CONNECT, weighted = TRUE, mode = "undirected", diag = FALSE)
      Zboot = c(Zboot, modularity(g, SUBJECT))
    }
    pvalue = length(which(Zboot<=Zorig))/NBOOT
    print(c(j,k, Zorig, pvalue))
  }
}


##### ALUNO COM OUTROS ALUNOS
for(j in 1:todos){ #for da dupla
  #Leitura de Dados
  for (k in 1:todos) {

    aluno2 = read.table(paste("./Data/", "aluno", j, "_001_oxyhb_bloco", bloco, ".txt", sep=""))
    aluno = read.table(paste("./Data/", "aluno", k, "_001_oxyhb_bloco", bloco, ".txt", sep=""))
    
    
        
    m = min(nrow(aluno2), nrow(aluno))
    aluno2 = aluno2[1:m,]
    aluno = aluno[1:m,]
    
    
    #GSR - pega os resíduos do sinal médio de cada canal dos professores e alunos
    GSR = rowMeans(aluno2)
    for(i in 1:ncol(aluno2)){aluno2[, i] = lm(aluno2[, i]~GSR)$resid}
    
    GSR = rowMeans(aluno)
    for(i in 1:ncol(aluno)){aluno[, i ] = lm(aluno[, i]~GSR)$resid}
    
    
    # Preprocessamento
    # Matriz de conectividade funcional - encontra a correlação dos canais
    # conexão entre cada canal
    # matriz 36x36: correlação dos 18 canais de professor com 18 canais de alunos.
    # todos os canais contra todos os canais e zerando corr<0.1
    
    CONNECT = cor(cbind(aluno2, aluno))   #,method="spearman")
    #para imprimir cada matriz write.csv(CONNECT, 'MATRIZ5.csv')
    
    CONNECT[which(CONNECT < 0.15)] = 0
    SUBJECT = c(rep(1, ncol(aluno2)), rep(2, ncol(aluno)))
    
    
    #Analise de Grafos a partir da modularidade (zorig)
    g = graph.adjacency(CONNECT, weighted=TRUE, mode="undirected", diag=FALSE)
    Zorig = modularity(g, SUBJECT)
    
    
    #Bootstrap
    Zboot = NULL
    NBOOT = 1000
    for(boot in 1:NBOOT){
      aux = sample(nrow(aluno2)-1, 1)
      auxsign = sample(c(-1, 1), 1)
      BOOTprof = aluno2[c((aux+1):nrow(aluno2), 1:aux),]*auxsign
      
      #Matriz de conectividade funcional
      CONNECT = cor(cbind(BOOTprof, aluno))  #,method="spearman")
      CONNECT[which(CONNECT < 0.15)] = 0
      SUBJECT = c(rep(1, ncol(aluno2)), rep(2, ncol(aluno)))
      g = graph.adjacency(CONNECT, weighted = TRUE, mode = "undirected", diag = FALSE)
      Zboot = c(Zboot, modularity(g, SUBJECT))
    }
    pvalue = length(which(Zboot<=Zorig))/NBOOT
    print(c(j,k, Zorig, pvalue))
  }
}

