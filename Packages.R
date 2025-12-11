# Imports and Packages

#General packages
install.packages("ggplot2")
install.packages("dplyr")
install.packages("readr")
install.packages("tidyr")
install.packages("corrplot")
install.packages("purrr")
install.packages("psych")
install.packages("GGally")
install.packages('writexl')
install.packages('stringr')
install.packages('rlang')
install.packages('farff')
install.packages('readxl')
install.packages('magrittr')

#dimRed
install.packages("dimRed")
install.packages("coRanking")


#DR methods
install.packages("diffusionMap")
install.packages("igraph")
install.packages("MASS")
install.packages("RSpectra")
install.packages("Matrix")
install.packages("RANN")
install.packages("pcaL1")
install.packages("DRR") 
install.packages("mds")
install.packages("fastICA")
install.packages("Rtsne")
install.packages("vegan")
install.packages("kernlab")

#DRQuality
install.packages("DRquality")
install.packages("remotes")
#remotes::install_github("Mthrun/ProjectionBasedClustering")
install.packages("FastKNN")
install.packages("FCPS")
install.packages("pcaPP")
install.packages("cccd")
install.packages("spdep")
install.packages("pracma")

# Main packages
library(dimRed)
library(ggplot2)
library(readr)
library(dplyr)
library(corrplot)
library(coRanking)
library(psych)
library(GGally)
library(purrr)
library(tidyr)
library(writexl)
library(stringr)
library(rlang)
library(farff)
library(readxl)
library(magrittr)

#DRQuality
library(DRquality)
library(ProjectionBasedClustering)
library(FastKNN)
library(FCPS)
library (pcaPP)
library(cccd)
library(spdep)
library(pracma)
#pkgbuild::has_build_tools(debug = TRUE)

# DR Methods
library(fastICA)     #ICA
library(Rtsne)       #t-SNE
library(vegan)       #Isomap, nMDS
library(diffusionMap)#Diffusion Maps
library(pcaL1)       #PCA L1
library(DRR)         #DRR
library(mds)         #MDS
library(igraph)      #KamadaKawai, DRL and FruchtermanReingold
library(kernlab)     #Kernel PCA
library(RANN)        #HLLE
library(Matrix)
library(RSpectra)
library(MASS)


