##### ANALISE DE Modularidade
require(caret)
library(fdrtool)
library(igraph)
require(lmtest)
require(png)
library(statGraph)
library(tidyverse)
require(xts)

analisa_bloco <- function(bloco)
{
  grafos = list()
  CONNECT_LIST = list()
  todos = c(1:12,14:25)
  matriz_centralidade = matrix(NA, nrow = max(todos), ncol = 32)
  num_arestas = NULL
  peso_arestas = NULL
  mods = NULL
  for(j in todos){
    prof = read.table(paste("./Data/", "prof", j, "_001_oxyhb_bloco", bloco, ".txt", sep=""))
    aluno = read.table(paste("./Data/", "aluno", j, "_001_oxyhb_bloco", bloco, ".txt", sep=""))
    GSR = rowMeans(prof)
    for(i in 1:ncol(prof)){prof[, i] = lm(prof[, i]~GSR)$resid}
    GSR = rowMeans(aluno)
    for(i in 1:ncol(aluno)){aluno[, i ] = lm(aluno[, i]~GSR)$resid}
    
    CONNECT = cor(cbind(prof, aluno))   #,method="spearman")
    #aux = c(rep("P-", 16), rep("A-", 16))
    #colnames(CONNECT) = paste0(aux, colnames(CONNECT))
    #rownames(CONNECT) = paste0(aux, rownames(CONNECT))
    
    CONNECT[which(CONNECT < 0.15)] = 0
    SUBJECT = c(rep(1, ncol(prof)), rep(2, ncol(aluno)))
    
    g = graph.adjacency(CONNECT, weighted=TRUE, mode="undirected", diag=FALSE)
    zorig = modularity(g, SUBJECT)
    mods = c(mods, zorig)
    
    CONNECT[1:16, 1:16] = 0
    CONNECT[17:32, 17:32] = 0
    num_arestas = c(num_arestas, sum(CONNECT != 0))
    peso_arestas = c(peso_arestas, sum(CONNECT))
    
    g = graph.adjacency(CONNECT, weighted=TRUE, mode="undirected", diag=FALSE)
    grafos = append(grafos, list(g))
    
    CONNECT_LIST = append(CONNECT_LIST, list(CONNECT != 0))
  }
  
  # Takahashi clusters
  # Takahashi, D. Y., Sato, J. R., Ferreira, C. E. and Fujita A. (2012) Discriminating Different Classes
  # of Biological Networks by Analyzing the Graph Spectra Distribution. _PLoS ONE_, *7*, e49949.
  # doi:10.1371/journal.pone.0049949.
  clust = graph.kmeans(grafos, 2) %>% as.numeric()
  aux = tibble(Dupla = todos, clust)
  
  # Medidas resumo do grafo
  data_grafo <- tibble(Dupla = todos, 
                       num_arestas = num_arestas,
                       peso_arestas = peso_arestas,
                       mods = mods)
  
  # Variáveis explicativas das duplas
  data_quali <- read_csv2("./Data/qualitativo_duplas.csv")
  
  # Takahashi-cluster (kmeans = 2) + explicativa + medidas resumo
  data = inner_join(data_grafo, data_quali, by = "Dupla") %>% 
    inner_join(aux, by = "Dupla")
  data = data %>% 
    mutate(escola_publica = (Aluno == "publico"),
           formacao_exatas = (Professor == "exatas"),
           gen_aluno_f = (Gen_aluno == "F"),
           gen_prof_f = (Gen_prof == "F")) %>% 
    select(num_arestas, peso_arestas, mods, dupla = Dupla, clust,
           escola_publica, formacao_exatas, gen_aluno_f, gen_prof_f)
   write_csv(data, paste0("./Data/resumo_duplas_", bloco, ".csv"))
  
  # Resumo dos clusters
  data %>% 
    group_by(clust) %>% 
    summarise(num_arestas = mean(num_arestas),
              peso_arestas = mean(peso_arestas), 
              mods = mean(mods),
              escola_publica = mean(escola_publica),
              formacao_exatas = mean(formacao_exatas),
              gen_aluno_f = mean(gen_aluno_f),
              gen_prof_f = mean(gen_prof_f)) %>% 
    print(width=Inf)
  
  # GLM: num_arestas por covariáveis
  data %>% 
    glm(num_arestas ~ escola_publica +  gen_aluno_f + gen_prof_f,
        family = poisson(link = "log"), data = .) %>% 
    summary() %>% 
    print()
  
  # LM: peso_arestas por covariáveis
  data %>% 
    lm(peso_arestas ~ escola_publica +  gen_aluno_f + gen_prof_f,
       data = .) %>% 
    summary() %>% 
    print()
  
  # LM: modularidade por covariáveis
  data %>% 
    lm(mods ~ escola_publica +  gen_aluno_f + gen_prof_f,
       data = .) %>% 
    summary() %>% 
    print()
  
  # Lista de CONNECTS por dupla
  CONNECT_LIST
}

# Função que gera o "grafo médio", mas não exatamente
gera_grafo_medio <- function(CONNECT_LIST, quantil = 0.9)
{
  MED_CONNECT = CONNECT_LIST[[1]]
  dims = dim(MED_CONNECT)
  for(ii in 1:dims[1])
  {
    for(jj in 1:dims[2])
    {
      con_i_j = NULL
      for(kk in 1:length(CONNECT_LIST))
      {
        con_i_j = c(con_i_j, CONNECT_LIST[[kk]][ii, jj])
      }
      MED_CONNECT[ii, jj] = mean(con_i_j)
    }
  }
  quant = quantile(MED_CONNECT, probs = quantil)
  g <- MED_CONNECT > quant
}
  
plot_connect <- function(CONNECT,
                         fname = "centrality_medio.pdf")
{
  molde <- readPNG("molde_hypper.png")
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
  
  par(mfrow =c(25, 2))
  todos = c(1:12,14:25)
  pdf(fname, width = 17, height = 7)
  
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

  dev.off()
}

# ANÁLISE BLOCO 3
CONNECT_LIST = analisa_bloco(3)
CONNECT_MEDIO = gera_grafo_medio(CONNECT_LIST, quantil = 0.95)
plot_connect(CONNECT_MEDIO, "centrality_medio_3.pdf")

# ANÁLISE BLOCO 4
CONNECT_LIST = analisa_bloco(4)
CONNECT_MEDIO = gera_grafo_medio(CONNECT_LIST, quantil = 0.95)
plot_connect(CONNECT_MEDIO, "centrality_medio_4.pdf")
