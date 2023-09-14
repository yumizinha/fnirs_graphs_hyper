require(igraph)
require(caret)
require(tidyverse)
require(png)
require(lmtest)
require(xts)

hyper_plot <- function(valores)
{
  molde <- readPNG("molde_hypper.png")
  
  # Calculando ajuste das coordenadas calculadas no eeg_positions
  # deslocando cérebro do aluno para separar as cabeças corretamente
  # Alterei o deslocamento para coincidir com o template do cérebro carregado
  coords_fnirs_prof = read.csv("coords_fNIRS.csv") %>%
    mutate(x = x - 1.52) %>%
    mutate(y = y - 0.028)
  
  coords_fnirs_aluno = read.csv("coords_fNIRS.csv") %>% 
    mutate(x = x + 1.53) %>%
    mutate(y = y - 0.028)
  
  # unindo as novas coordenadas para plot do grafo
  coords_fnirs = rbind(coords_fnirs_prof, coords_fnirs_aluno) %>% 
    select(x, y) %>% 
    as.matrix()
  
  dd <- length(valores)
  CONNECT <- matrix(0, dd, dd)
  colnames(CONNECT) <- names(valores)
  rownames(CONNECT) <- names(valores)
  g = graph.adjacency(CONNECT, weighted = TRUE, mode = "undirected", diag = FALSE)
  
  # configurando paleta de cor por centralidade
  fine = 500 # this will adjust the resolving power.
  pal = colorRampPalette(c('red','green'))
  #this gives you the colors you want for every point
  graphCol = pal(fine)[as.numeric(cut(valores, breaks = fine))]
  
  # colocar centrality$vector como variavel em vertex.color
  plot(g, vertex.color=graphCol, layout = coords_fnirs, vertex.size = 10, rescale = FALSE)
  
  # Tentando colocar o pano de fundo do hypper
  lim <- par()
  rasterImage(molde,
              xleft=-2.6, xright=2.6, 
              ybottom=-1.1, ytop=1.1)
}
