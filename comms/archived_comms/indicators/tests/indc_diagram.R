library(threed)
library(dplyr)
library(ggplot2)

rm(list = ls())

mesh3dobj$cube %>%
  as.data.frame()

eye <- look_at_matrix(eye = c(6, 2, 9), at = c(0, 0, 0))

obj1 <- mesh3dobj$cube %>%
  transform_by(invert_matrix(eye)) %>%
  perspective_projection() %>%
  as.data.frame()

ggplot(obj1) +
  geom_polygon(aes(x = x, y = y, group = zorder, fill = 0.5 * fnx + fny), colour = 'black', size = 0.2) +
  theme_minimal() +
  theme(
    legend.position = 'none') + #,
    #axis.text = element_blank()) +
  coord_equal()


cc <- mesh3dobj$cube
tm <- matrix(0, nrow = 4, ncol = 8)
tm[2,] <- 2

cc$vb <- cc$vb + tm

obj2 <- cc %>%
  transform_by(invert_matrix(eye)) %>%
  perspective_projection() %>%
  as.data.frame() %>%
  mutate(element_id = element_id + 6,
         zorder = factor(as.numeric(as.character(zorder)) + 6))

ggplot(obj2) +
  geom_polygon(aes(x = x, y = y, group = zorder, fill = 0.5 * fnx + fny), colour = 'black', size = 0.2) +
  theme_minimal() +
  theme(
    legend.position = 'none') + #,
  #axis.text = element_blank()) +
  coord_equal()

rbind(obj1, obj2) %>%
  ggplot() +
  geom_polygon(aes(x = x, y = y, group = zorder, fill = 0.5 * fnx + fny), colour = 'black', size = 0.2) +
  theme_minimal() +
  theme(
    legend.position = 'none') + #,
  #axis.text = element_blank()) +
  coord_equal()


obj1 <- mesh3dobj$cube %>%
  as.data.frame()
cc <- mesh3dobj$cube
tm <- matrix(0, nrow = 4, ncol = 8)
tm[2,] <- 0.5

cc$vb <- cc$vb + tm

obj2 <- cc %>%
  transform_by(invert_matrix(eye)) %>%
  perspective_projection() %>%
  as.data.frame()

ggplot(obj1) +
  geom_polygon(aes(x = x, y = y, group = zorder, fill = 0.5 * fnx + fny), colour = 'black', size = 0.2) +
  theme_minimal() +
  theme(
    legend.position = 'none') + #,
  #axis.text = element_blank()) +
  coord_equal()



# DiagrammeR experiment ---------------------------------------------------


library(DiagrammeR)

a_graph <- create_graph() %>%
    add_node(label = "species 1",
             node_aes = node_aes(shape = "rectangle")) %>%
    add_node(label = "species 1",
             node_aes = node_aes(shape = "rectangle")) %>%
    add_node(label = "species 2") %>%
    add_edge(from = 1, to = 2)

a_graph <- create_graph() %>%
    add_n_nodes(label = paste0("species", 1:2),
                n = 2,
                node_aes = node_aes(shape = "rectangle")) %>%
    add_n_nodes(label = paste0("year", 1:2),
                n = 2,
                node_aes = node_aes(shape = "rectangle")) %>%
    add_n_nodes(label = paste0("year", 1:2),
                n = 2,
                node_aes = node_aes(shape = "rectangle")) %>%
    add_edge(from = 1, to = 3) %>%
    add_edge(from = 1, to = 4) %>%
    add_edge(from = 2, to = 5) %>%
    add_edge(from = 2, to = 6)

render_graph(a_graph, layout = "tree")

grViz("
      digraph Random{
      graph [
      layout = circo,
      overlap =T,
      outputorder = edgesfirst,
      bgcolor='white',
      splines=line]

      # controls l type setup
      edge[labelfontname='Arial',fontSize=13,color='red',fontcolor='navy']
      node [shape = box,style='filled',
      fillcolor='indianred4',width=2.5,
      fontSize=20,fontcolor='snow',
      fontname='Arial']

      # node shape
      a [label = 'A']
      b [label = 'B']
      c [label='D']
      a->b[color='red']
      b->c[color='dodgerblue']
      }")


grViz("
      digraph levels{
      graph [
      layout = circo]

      # controls l type setup
      edge[labelfontname='Arial',fontSize=13,color='red',fontcolor='navy']
      node [shape = box,style='filled',
      fillcolor='indianred4',width=2.5,
      fontSize=20,fontcolor='snow',
      fontname='Arial']

      # node shape
      a [label = 'species 1']
      b [label = 'species 2']
      c [label = 'year 1']
      a->c[color='red']
      b->c[color='dodgerblue']
      }")

