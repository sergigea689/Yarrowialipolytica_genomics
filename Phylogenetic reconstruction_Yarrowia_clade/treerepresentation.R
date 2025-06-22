#Representing phylogenetic tree and annotating it

#Install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
install.packages("ape")

#Libraries
library(ggtree)
library(ape)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(readxl)
library(reshape2)
library(openxlsx)

#Directories
#Cambiar el directorio
getwd()
setwd("C:/Users/USUARIO/Documents/GitHub/Yarrowialipolytica_genomics/Phylogenetic reconstruction_Yarrowia_clade")

#Load the tree in newick format
tree <- read.tree("final_aligned.fasta.treefile")  # Carga un archivo en formato Newick

#Fast visualization
ggtree(tree) + geom_tiplab()
ggtree(tree, layout="circular")
#Load metadata to annotate
popgen <- read_excel("update_manfinetreeposition.xlsx")


#store tip names
tips <- tree$tip.label
#Change to meaningful names

names <- c("CBS 6659", "CBS 10144", "MUCL 42905", "H222", "888", "DBVPG 4400",
           "Y. galli", "Y. brassicae", "Y. divulgata", "Y. deformans", "C. hispaniensis", "N. commutata",
           "Y. alimentaria", "Y. phangngaensis", "Y. bubula", "Y. hollandica", "Y. keelungensis", "Y. porcina",
           "Y. osloensis", "Y. yakushimensis", "NCYC 3071", "CLIB 122", "DBVPG 4557", "CLIB 84")

tree$tip.label <- names
ggtree(tree) + geom_tiplab()

#Rooting
"N. commutata" %in% tree$tip.label
tree_rooted <- root(tree, outgroup = "N. commutata", resolve.root = TRUE)
ggtree(tree_rooted) +
  geom_tiplab() +
  ggtitle("Árbol enraizado en N. commutata")

#Genetic distances
library(tibble)
library(dplyr)


# Calcular matriz de distancias genéticas entre tips
dist_matrix <- cophenetic.phylo(tree_rooted)

# Extraer distancia desde N. commutata a todos los demás tips
dist_to_ncommutata <- dist_matrix["N. commutata", ]

# Convertir a data.frame para unir al árbol
dist_df <- data.frame(
  label = names(dist_to_ncommutata),
  dist_to_ncommutata = as.numeric(dist_to_ncommutata)
)

#Colour branches according to genetic distances

ggtree(tree_rooted) %<+% dist_df +
  geom_tree(aes(color = dist_to_ncommutata), size = 0.9) +
  geom_tiplab(size = 4, fontface = "italic") +
  scale_color_viridis_c(
    option = "inferno",    # más contraste en rangos estrechos
    direction = -1,
    end = 0.98,
    limits = c(2.15, 2.42),
    breaks = c(2.15, 2.25, 2.35, 2.42),
    name = "Genetic\ndistance to\nN. commutata"
  ) +
  theme_tree2(base_size = 14) +
  ggtitle("Genetic distance to *N. commutata* across species")

#Heatmap
phylo <- as.phylo(tree_rooted)
rownames(dist_df) <- dist_df$label

ggtree(phylo) %<+% dist_df -> p



dist_df <- dist_df %>%
  mutate(
    recent_label = case_when(
      dist_to_ncommutata == 0 ~ "Outgroup",
      dist_to_ncommutata < 2.234605 ~ "Closest species",
      dist_to_ncommutata > 2.234605 ~ "Most distant"
    ),
    species_group = if_else(str_detect(label, fixed(".")), "Other species", "Y. lipolytica")
  )
# Colores pastel para recent_label
colors_recent <- c(
  "Closest species" = "#a6cee3",  # azul pastel
  "Most distant"  = "#fb9a99",  # rojo pastel
  "Outgroup"       = "#b2df8a"   # verde pastel
)

# Colores pastel para species_group
colors_species <- c(
  "Y. lipolytica" = "#c78734",  # naranja pastel
  "Other species" = "#85549e"   # lila pastel
)

# Representación con gheatmap usando los bloques
# Asumiendo que 'tree_rooted' es el árbol y 'dist_df' está ordenado por el orden de tree$tip.label

p <- ggtree(tree_rooted) %<+% dist_df +
  geom_tiplab(aes(color = species_group), size = 6, fontface = "italic") +
  scale_color_manual(values = colors_species, name = "Species group") +
  theme_tree2(base_size = 14) +
  ggtitle("Phylogenetic tree colored by species group")

p_heatmaps <- gheatmap(p, 
                       dist_df[, "recent_label", drop = FALSE],
                       offset = 0.2,  # más a la derecha
                       width = 0.05, 
                       colnames = FALSE) +
  scale_fill_manual(values = colors_recent, name = "Distance to outgroup") +
  theme_tree2(base_size = 14) +
  ggtitle("Phylogenetic tree with species and genetic group annotation")

print(p_heatmaps)
#Subset metadata to strains in tree
anno <- subset(popgen, popgen$treesamp %in% tree$tip.label)
anno$morph <- as.factor(anno$morph)

#Check coinciddences among sample labels

setdiff(tree$tip.label, anno$treesamp)


#Represent tree plus annotations
p <- ggtree(tree, layout="circular") %<+% anno
# p <- p + geom_tiplab()

color_palette_habitat <- c("food" = "#62abe2", "host" = "#ff3131", "industry" = "#505050", "marine"= "#004aad", "soil" = "#c85700")

# Definir una paleta de formas
color_palette_continent <- c("Europe" = "gold2", "Africa" = "gold2", "North America" = "purple4", "South America" = "red4", "Asia" = "green4")  # 16: círculo, 17: triángulo

colormorph <- c("1" = "gold", "2" = "green", "3" = "cyan")
# p + geom_tippoint(aes(colour = habitat, shape = continent), size = 3) +
#   scale_color_manual(values = color_palette) +  # Asignar colores específicos
#   scale_shape_manual(values = shape_palette) +
#   theme_tree()  # Mejora la visualización del árbol

p + 
  geom_tippoint(aes(colour = continent), size = 3, position = position_nudge(x = 0.005)) +  # Primera bolita desplazada
  geom_tippoint(aes(colour = habitat), size = 3, position = position_nudge(x = -0.005)) +  # Segunda bolita desplazada
  scale_color_manual(values = c(color_palette_continent, color_palette_habitat)) +  # Colores para ambas bolitas
  theme_tree()  # Mejorar la visualización del árbol

treeim <- p + geom_tippoint(aes(colour = morph), size = 3, position = position_nudge(x = 0.005)) +
  scale_color_manual(values = c(colormorph))
ggsave("morphtree.png", plot = treeim, width = 10, height = 10, dpi = 300)

p + 
  geom_tippoint(aes(colour = habitat), size = 3, position = position_nudge(x = 0.001)) +  # Primera bolita desplazada
  geom_tippoint(aes(shape = continent), size = 2, position = position_nudge(x = 0.015)) +  # Segunda bolita desplazada
  scale_color_manual(values = c(color_palette_habitat)) +  # Colores para ambas bolitas
  scale_shape_manual(values = color_palette_continent) +
  theme_tree()  # Mejorar la visualización del árbol


#Let's try and heatmap the metadata

# Ejemplo: Datos de anotación para el heatmap
data_heatmap <- as.data.frame(anno[, c(1, 11, 12)])

colnames(data_heatmap) <- c("label", "continent", "habitat")
rownames(data_heatmap) <- data_heatmap$label

#Create a sub dataframe with habitat features

hab <- data_heatmap[, c(1, 3)] 
annohab <- as.data.frame(hab[, -1])
colnames(annohab) <- "habitat"
rownames(annohab) <- hab$label

#Create a sub dataframe with continent features

cont <- data_heatmap[, c(1, 2)] 
annocont <- as.data.frame(cont[, -1])
colnames(annocont) <- "continent"
rownames(annocont) <- cont$label

#Cretae the tree graph

p <- ggtree(tree, layout = "circular") +geom_tree(color = "black", size = 0.3) +
  geom_treescale(x = 0.5, y = -1, width = 0.05, fontsize = 3)

class(p)
# Add the heatmap of one feature
p <- gheatmap(p, annohab, offset = 0.0001, width = 0.04, 
              colnames_angle = 0, colnames_offset_y = 0.1)+ 
  scale_fill_manual(values = color_palette_habitat, na.value = "grey")

p



p <-  gheatmap(p, annocont, offset = 0.03, width = 0.04, 
              colnames_angle = 0, colnames_offset_y = 0.1)+ 
  scale_fill_manual(values = c(color_palette_continent, color_palette_habitat), na.value = "grey")

p

colnames(data_heatmap)

ggsave("high_res_tree.png", plot = p, width = 10, height = 10, dpi = 300)


#Rotate to have the tree in the same position as other analysis
##Identify nodes
tree <- as.phylo(tree)
g <- ggtree(tree, layout = "rectangular") + 
  geom_nodelab(aes(subset=!isTip, label=node), hjust = -.1, color = "red")  # Etiqueta los nodos internos
g
#Rotate nodes

#We need functions because there is a bug
nodeid.tbl_tree <- utils::getFromNamespace("nodeid.tbl_tree", "tidytree")
rootnode.tbl_tree <- utils::getFromNamespace("rootnode.tbl_tree", "tidytree")
offspring.tbl_tree <- utils::getFromNamespace("offspring.tbl_tree", "tidytree")
offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item", "tidytree")
child.tbl_tree <- utils::getFromNamespace("child.tbl_tree", "tidytree")
parent.tbl_tree <- utils::getFromNamespace("parent.tbl_tree", "tidytree")

rotree <- rotate(tree, node = 129)
rotree <- rotate(rotree, node = 143)


g <- ggtree(rotree, layout = "circular")   # Etiqueta los nodos internos
g
flip(g, 130, 144)
g
rotate
class(tree)
# flip
g <- ggtree(tree)
g
g1 <- flip(ggtree(tree), node1 = 130, node2 = 144)
tree_rotated <- flip(tree, 130, 144)


# library(msa)
# library(phangorn)
# library(ggtree)
# library(phytools)
# library(Biostrings)
# install.packages("phangorn")

#We want to get the tree order for the labels
tip_labels <- g$data$label

la <- ggtree(tree) +
  geom_tiplab(size = 1)
la




#Cálculo de la distancia genética

library(ape)
  # Cargar tu árbol filogenético
dist <- cophenetic.phylo(tree)       # Matriz de distancias de cophenetic
dist["cepa1", "cepa2"]     

write.csv(as.matrix(dist), "distancias_geneticas.csv", row.names = TRUE)
write.table(as.matrix(dist), "distancias_geneticas.tsv", sep = "\t", quote = FALSE, row.names = TRUE)

#Visualize the matrix

library(pheatmap)
pheatmap(as.matrix(dist))

#Obtaining genetic distances in suitable format
df_distancias <- as.data.frame(as.table(dist))
colnames(df_distancias) <- c("Strain1", "Strain2", "Geneticdistance")


#Calculate divergence times

# Mutation rate in Yarrowia lipolytica
mu <- 4e-10  

# Estimate divergence time in terms of generations that have passed
df_distancias$divergence_time <- df_distancias$Geneticdistance / (2 * mu)

# Show results
head(df_distancias)

#Translate these results into years

# Two cases: min (1 division/day) y max (8 divisions/day)
generaciones_por_anio_min <- 365      # 1 vez al día
generaciones_por_anio_max <- 2920     # 8 veces al día

# Calcular el tiempo en años para ambos casos
df_distancias$Tiempo_Años_min <- df_distancias$divergence_time / generaciones_por_anio_max
df_distancias$Tiempo_Años_max <- df_distancias$divergence_time / generaciones_por_anio_min


#Gráfico de dispersión
# Suponiendo que ya tienes el dataframe df_distancias con las columnas 'Distancia_Genetica', 'Tiempo_Anios_min' y 'Tiempo_Anios_max'

ggplot(df_distancias, aes(x = Geneticdistance, y = Tiempo_Años_min)) +
  geom_point(color = "blue", alpha = 0.5) +   # Puntos para el tiempo de divergencia mínimo
  geom_point(aes(y = Tiempo_Años_max), color = "red", alpha = 0.5) +  # Puntos para el tiempo de divergencia máximo
  labs(x = "Distancia Genética", y = "Tiempo de Divergencia (Años)") +
  theme_minimal() +
  ggtitle("Visualización de Tiempos de Divergencia vs Distancia Genética")


#RelTime time calibrated tree

# Cargar las bibliotecas necesarias
library(ape)

# Paso 1: Cargar el árbol con tiempos relativos (asumiendo que el archivo es en formato Newick)
# Reemplaza el archivo de entrada con la ubicación de tu archivo
tree <- read.tree("cladetimecalibratedrelativetree")

# Paso 2: Especificar la tasa de mutación en sustituciones por base por división celular
mutation_rate <- 4e-10  # Tasa de mutación en sustituciones por base por ciclo de división

# Paso 3: Especificar el rango de generaciones por día (1 a 8 generaciones por día)
generations_per_day_min <- 1
generations_per_day_max <- 8

# Paso 4: Calcular las longitudes de las ramas en generaciones
# Usamos la tasa de mutación para convertir los tiempos relativos a generaciones
# El tiempo relativo en el árbol necesita ser multiplicado por un factor basado en la tasa de mutación

tree$edge.length_generations <- tree$edge.length / mutation_rate

# Paso 5: Convertir las generaciones a años usando el rango de generaciones por día
# Convertir generaciones a días: 1 generación = 1 día / (generaciones por día)
# Convertir días a años: 1 año = 365 días

tree$edge.length_years_min <- tree$edge.length_generations / generations_per_day_max / 365
tree$edge.length_years_max <- tree$edge.length_generations / generations_per_day_min / 365

# Paso 6: Verificar los tiempos calculados en años
# Mostramos las longitudes de las ramas convertidas a años (rango mínimo y máximo)
print(tree$edge.length_years_min)
print(tree$edge.length_years_max)
# Paso 7: Guardar el árbol con los tiempos absolutos en años (rango mínimo y máximo)
write.tree(tree, file = "arbol_absoluto_anios_min_max.nwk")

# Opcional: Graficar el árbol para visualizar los tiempos
plot(tree, main = "Árbol con tiempos en años (rango mínimo y máximo)")


# Paso 1: Crear un dataframe con las longitudes de las ramas y las etiquetas correspondientes
edge_data <- data.frame(
  node1 = tree$edge[, 1], 
  node2 = tree$edge[, 2], 
  time_years_min = tree$edge.length_years_min, 
  time_years_max = tree$edge.length_years_max
)

# Paso 2: Inicializar los vectores para los tiempos de los nodos
# El número de nodos es la cantidad de hojas + la cantidad de nodos internos
node_times_min <- rep(0, length(tree$tip.label) + length(tree$node.label))  
node_times_max <- node_times_min

# Paso 3: Sumar los tiempos de las ramas a los nodos correspondientes
for (i in 1:nrow(tree$edge)) {
  node_times_min[tree$edge[i, 2]] <- node_times_min[tree$edge[i, 2]] + tree$edge.length_years_min[i]
  node_times_max[tree$edge[i, 2]] <- node_times_max[tree$edge[i, 2]] + tree$edge.length_years_max[i]
}

# Paso 4: Crear las matrices de distancias entre los nodos
# Usamos las diferencias de tiempos entre los nodos
distance_matrix_min <- outer(node_times_min, node_times_min, FUN = function(x, y) abs(x - y))
distance_matrix_max <- outer(node_times_max, node_times_max, FUN = function(x, y) abs(x - y))

# Paso 5: Convertir las matrices de distancias a un dataframe
df_distance_min <- as.data.frame(distance_matrix_min)
df_distance_max <- as.data.frame(distance_matrix_max)

# Paso 6: Etiquetar las filas y columnas con los nombres correctos
# Las primeras filas corresponden a las hojas (tipos), y las siguientes a los nodos internos
tip_labels <- tree$tip.label
node_labels <- paste("Node", (length(tree$tip.label) + 1):(length(tree$tip.label) + length(tree$node.label)), sep = "_")

# Paso 7: Asegurarse de que las matrices de distancias estén correctamente dimensionadas (24x24 en tu caso)
df_distance_max <- df_distance_max[1:24, 1:24]
df_distance_min <- df_distance_min[1:24, 1:24]
# Unir las etiquetas de las hojas y nodos
colnames(df_distance_min) <- c(tip_labels, node_labels)
rownames(df_distance_min) <- c(tip_labels, node_labels)
colnames(df_distance_max) <- c(tip_labels, node_labels)
rownames(df_distance_max) <- c(tip_labels, node_labels)



# Asignar nombres a las filas y columnas (nombres de los nodos y tips)
rownames(df_distance_min) <- colnames(df_distance_min) 
rownames(df_distance_max) <- colnames(df_distance_max) 


# Paso 10: Crear el heatmap usando pheatmap
# Heatmap de tiempos mínimos (en años)
library(pheatmap)
pheatmap(df_distance_min, main = "Heatmap de Distancias - Tiempos Mínimos (Años)")

# Heatmap de tiempos máximos (en años)
pheatmap(df_distance_max, main = "Heatmap de Distancias - Tiempos Máximos (Años)")


#Antoher way for the RELTime calibrated tree
# Cargar las librerías necesarias
library(ape)
library(pheatmap)

# Paso 1: Cargar el árbol con tiempos relativos (formato Newick)
tree <- read.tree("cladetimecalibratedrelativetree")

# Paso 2: Especificar la tasa de mutación (sustituciones por base por división celular)
mutation_rate <- 4e-10  

# Paso 3: Especificar el rango de generaciones por día (1 a 8)
gen_day_min <- 1
gen_day_max <- 8

# Paso 4: Convertir tiempos relativos a generaciones
# Suponemos que los valores en tree$edge.length son proporcionales al número de sustituciones.
# Dividimos por la tasa de mutación para obtener generaciones.
tree$edge.length_generations <- tree$edge.length / mutation_rate

# Paso 5: Convertir generaciones a años
# Usamos el límite más restrictivo (mayor número de generaciones por día) para obtener el tiempo mínimo en años.
tree$edge.length_years_min <- tree$edge.length_generations / gen_day_max / 365
# Para obtener el tiempo máximo (menos generaciones por día):
tree$edge.length_years_max <- tree$edge.length_generations / gen_day_min / 365

# Para el ejemplo, vamos a trabajar con los tiempos mínimos.
tree_abs <- tree
tree_abs$edge.length <- tree_abs$edge.length_years_min  # Asignamos los tiempos en años (mínimos)
tree_abs$edge.length <- tree_abs$edge.length_years_max#do when finished

# Tree is already in absolute times
# Ahora calculamos los tiempos de divergencia entre los nodos (en años)
branching_times <- branching.times(tree_abs)  # Obtiene los tiempos de los nodos

# Imprimir los tiempos de los nodos para asegurarnos de que son absolutos
print(branching_times)

# Paso 7: Obtener la matriz de MRCA para los tips
mrca_matrix <- mrca(tree_abs)  # Encontramos los nodos más recientes comunes (MRCA)

# Paso 8: Calcular la matriz de divergencia entre los tips
tip_labels <- tree_abs$tip.label
n_tips <- length(tip_labels)
divergence_matrix <- matrix(0, n_tips, n_tips)
rownames(divergence_matrix) <- tip_labels
colnames(divergence_matrix) <- tip_labels

# Usar los tiempos absolutos para calcular la divergencia
for (i in 1:n_tips) {
  for (j in 1:n_tips) {
    if (i == j) {
      divergence_matrix[i, j] <- 0  # No hay divergencia entre el mismo tip
    } else {
      mrca_node <- mrca_matrix[i, j]  # Encontramos el nodo MRCA para este par de tips
      divergence_matrix[i, j] <- branching_times[as.character(mrca_node)]  # Tiempo de divergencia en años
    }
  }
}

# Convertir la matriz de divergencia a un dataframe
df_divergence <- as.data.frame(divergence_matrix)

# Paso 9: Visualizar el heatmap del tiempo de divergencia
pheatmap(as.matrix(df_divergence), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         main = "Tiempo de Divergencia (en años)")

# Visualizar el árbol con los tiempos de los nodos
library(ggtree)
plot(tree_abs, main = "Árbol con tiempos de nodos calibrados")
nodelabels(branching_times)  # Etiquetar los nodos con sus tiempos de divergencia

p <- ggtree(tree_abs, layout = "rectangular", branch.length = "time", 
            main = "Árbol con tiempos de divergencia (en años)") +
  geom_tiplab(size = 3, hjust = -0.1) +  # Etiquetas de las tips (pequeñas y visibles)
  geom_nodelab(aes(label = branching_times[as.character(node)]), size = 3, color = "blue") +  # Etiquetas de los nodos
  theme_tree2() +  # Mejorar el tema del árbol
  theme(plot.title = element_text(size = 16, face = "bold"))  # Título grande para publicación

# Mostrar el gráfico
print(p)

#Representing KYA
library(ggtree)
library(ggplot2)



# Asegúrate de que las coordenadas de los nodos se calculen correctamente
p <- ggtree(tree_abs, layout = "rectangular", branch.length = "time") +
  geom_tiplab(size = 3, hjust = -0.1) +  # Etiquetas de las tips (pequeñas y visibles)
  geom_nodelab(aes(label = paste0(round(branching_times[as.character(node)] / 1000))), 
               size = 4, color = "blue", 
               nudge_y = 0.2, nudge_x = 0.5) +  # Etiquetas de los nodos en KYA, desplazadas
  ggtitle("Árbol con tiempos de divergencia (en KYA)") +  # Título del gráfico
  theme_tree2() +  # Mejorar el tema del árbol
  theme(plot.title = element_text(size = 16, face = "bold"))  # Título grande para publicación

# Mostrar el gráfico
print(p)

