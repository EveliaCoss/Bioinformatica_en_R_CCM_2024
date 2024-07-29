# Mini tutorial de phyloseq
# Author: Daniel Vaulot
# Modificado por: Evelia Coss
# 29 Julio 2024
# Tutorial basado en el tutorial https://vaulot.github.io/tutorials/Phyloseq_tutorial.html

# --- Cargar paquetes -----
library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("tidyverse")       # Needed for converting column to row names,  filter and reformat data frames

# --- Paso 1. Descarga y descomprir datos ---- 
#  Descarga y descomprime la siguiente carpeta https://github.com/vaulot/R_tutorials/archive/master.zip, guardala en `data/`
# Dentro de tu R project.

getwd() # Verificar ubicacion
# Cambiar de ubicacion si es necesario en Session/Set Working Directory/Choose directory/

# --- Paso 2. Importar datasets ----
otu_mat<- read_excel("../data/CARBOM data.xlsx", sheet = "OTU matrix")
tax_mat<- read_excel("../data/CARBOM data.xlsx", sheet = "Taxonomy table")
samples_df <- read_excel("../data/CARBOM data.xlsx", sheet = "Samples")

# --- Paso 3. Renombrar objeto phyloseq ----

# define the row names from the otu column
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("otu") 

# Idem for the two other matrixes
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("otu")

samples_df <- samples_df %>% 
  tibble::column_to_rownames("sample") 

# Verificar formato
class(otu_mat)
class(tax_mat)

# --- Paso 4.Cambiar formato de data.frame a matrix ------------
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

# Verificar formato
class(otu_mat)
class(tax_mat)

# --- Paso 5. Transformar a un objeto phyloseq ------------

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

# Verificar archivos generados
head(OTU,3) # Muestras
head(TAX,3) # filogenia
head(samples,3)

carbom <- phyloseq(OTU, TAX, samples)
carbom

# Verificar nombre de las muestras
sample_names(carbom)

# Tipo de informacion contenida 
rank_names(carbom)

# Nombre de las columnas (variables)
sample_variables(carbom)

# --- Paso 6. Seleccionar solo las muestras que se van a analizar ------------
carbom <- subset_samples(carbom, Select_18S_nifH =="Yes")
carbom

# Seleccionar algunas filogenias de interes
carbom <- subset_taxa(carbom, Division %in% c("Chlorophyta", "Dinophyta", "Cryptophyta", 
                                              "Haptophyta", "Ochrophyta", "Cercozoa"))
carbom <- subset_taxa(carbom, !(Class %in% c("Syndiniales", "Sarcomonadea")))
carbom

# --- Paso 7. Normalizar el número de lecturas en cada muestra utilizando la media de la profundidad de secuenciación -----

total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)


# --- Paso 8. Visualizacion grafica --------

# Grafica sencilla
plot_bar(carbom, fill = "Division")

# Grafica acoplada a ggplot2
plot_bar(carbom, fill = "Division") +  # base
  geom_bar(aes(color=Division, fill=Division), stat="identity", position="stack") +
  scale_colour_brewer(palette = "Set1")

# Reagrupar las fracciones de las muestras, obteniendo 2 grupos (Pico vs Nano)
carbom_fraction <- merge_samples(carbom, "fraction")
plot_bar(carbom_fraction, fill = "Division") + 
  geom_bar(aes(color=Division, fill=Division), stat="identity", position="stack")

#  Seleccionar solo Chlorophyta
carbom_chloro <- subset_taxa(carbom, Division %in% c("Chlorophyta"))
sample_variables(carbom_chloro)
plot_bar(carbom_chloro, x="Genus", # x = Tax
         fill = "Genus", # colorear por Genero 
         facet_grid = level~fraction) + # dividir por level y fraction
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") +
  theme_bw() 

# Mas temas aqui https://ggplot2.tidyverse.org/reference/ggtheme.html

# Heatmap
plot_heatmap(carbom, method = "NMDS", distance = "bray")



