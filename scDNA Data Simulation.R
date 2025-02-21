#####################################################################
#####################################################################
####################### scDNA Data Simulation #######################
############################ Xiangyu Xu #############################
#####################################################################
#####################################################################

###### Set Workspace ######
rm(list = ls())
this.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(this.dir)

###### Packages ######


###### Dataset ######
AlleleFreq_all = read.csv("AlleleFreq_all.csv", header = TRUE)
AlleleFreq_all[is.na(AlleleFreq_all)] = 0
AlleleFreq_AfAm = read.csv("AlleleFreq_AfAm.csv", header = TRUE)
AlleleFreq_AfAm[is.na(AlleleFreq_AfAm)] = 0
AlleleFreq_Cauc = read.csv("AlleleFreq_Cauc.csv", header = TRUE)
AlleleFreq_Cauc[is.na(AlleleFreq_Cauc)] = 0
AlleleFreq_Hispanic = read.csv("AlleleFreq_Hispanic.csv", header = TRUE)
AlleleFreq_Hispanic[is.na(AlleleFreq_Hispanic)] = 0
AlleleFreq_Asian = read.csv("AlleleFreq_Asian.csv", header = TRUE)
AlleleFreq_Asian[is.na(AlleleFreq_Asian)] = 0

#####################################################################
############################# Version 1 #############################
#####################################################################

# Note 1: This version simulates diploid cells only.
# Note 2: This version does not consider the co-ancestry coefficient.
# Note 3: This version does not consider mutation.
# Note 4: A heterozygous locus with ADI has to have one ADO involved.
# Note 5: Heterozygous and homozygous are considered the same way in the case of ADO and ADI.

scDNASimulation_v1 = function(){
  
  AlleleFreq_all = read.csv("AlleleFreq_all.csv", header = TRUE)
  AlleleFreq_all[is.na(AlleleFreq_all)] = 0
  AlleleFreq_AfAm = read.csv("AlleleFreq_AfAm.csv", header = TRUE)
  AlleleFreq_AfAm[is.na(AlleleFreq_AfAm)] = 0
  AlleleFreq_Cauc = read.csv("AlleleFreq_Cauc.csv", header = TRUE)
  AlleleFreq_Cauc[is.na(AlleleFreq_Cauc)] = 0
  AlleleFreq_Hispanic = read.csv("AlleleFreq_Hispanic.csv", header = TRUE)
  AlleleFreq_Hispanic[is.na(AlleleFreq_Hispanic)] = 0
  AlleleFreq_Asian = read.csv("AlleleFreq_Asian.csv", header = TRUE)
  AlleleFreq_Asian[is.na(AlleleFreq_Asian)] = 0
  
  NOC =  as.numeric(readline(prompt = "Number of contributors: "))
  if (is.na(NOC)){
    NOC = 2 # If the user did not enter, NOC is by default to be 2.
  } else if (NOC <= 0 | NOC != as.integer(NOC)){
    stop("NOC has to be a positive integer.")
  }
  
  D = as.numeric(readline(prompt = "ADO rate D: "))
  if (D < 0 | D > 1){
    stop("D must be between 0 and 1")
  }
  
  e = as.numeric(readline(prompt = "ADI rate e: "))
  if (e < 0 | D > 1){
    stop("e must be between 0 and 1")
  }
  
  scDNASimulation = vector("list", length = NOC)
  
  for (i in 1:NOC){
    
    ind_pop = as.numeric(readline(prompt = paste("Subpopulation of contributor ", i, " (1: African American; 2: Caucasian; 3: Hispanic; 4: Asian; 5: All): ", sep = "")))
    if (!(ind_pop %in% c(1, 2, 3, 4))){ind_pop = 5} # If not specified, the genotype will be simulated from the entire population.
    population = switch(ind_pop, "African American", "Caucasian", "Hispanic", "Asian", "all")
    AlleleFreq = switch(ind_pop, AlleleFreq_AfAm, AlleleFreq_Cauc, AlleleFreq_Hispanic, AlleleFreq_Asian, AlleleFreq_all)
    
    num_cells = as.numeric(readline(prompt = paste("The number of cells from contributor ", i, ": ", sep = "")))
    if (is.na(num_cells)){
      num_cells = 100 # If the user did not enter, number of cells is by default to be 100.
    } else if (num_cells <= 0 | num_cells != as.integer(num_cells)){
      stop("The number of cells has to be a positive integer.")
    }
    
    genotype = array(NA, dim = c(num_cells, 20, 2))
    dimnames(genotype) = list(paste0("Cell_", 1:num_cells), colnames(AlleleFreq[, -1]), c("Allele_1", "Allele_2"))
    genotype_true = array(NA, dim = c(20, 2))
    dimnames(genotype_true) = list(colnames(AlleleFreq[, -1]), c("Allele_1", "Allele_2"))
    
    for (k in 1:20){
      AlleleFreq_k = AlleleFreq[, c(1, k+1)]
      AlleleFreq_k = AlleleFreq_k[AlleleFreq_k[, 2] != 0, ]
      alleles = sample(AlleleFreq_k[, 1], size = 2, replace = TRUE, prob = AlleleFreq_k[, 2])
      alleles_true = alleles
      genotype_true[k, ] = alleles_true
      for (j in 1:num_cells){
        alleles = alleles_true
        for (m in 1:2){
          if (D > 0){
            if (runif(1, 0, 1) <= D){
              alleles[m] = 0
              if (e > 0){
                if(runif(1, 0, 1) <= e){
                  AlleleFreq_temp = AlleleFreq_k[AlleleFreq_k[, 1] != alleles_true[m], ]
                  alleles[m] = sample(AlleleFreq_temp[, 1], size = 1)
                }
              }
            }
          }
        }
        genotype[j, k, ] = alleles
      }
    }
    
    scDNASimulation[[i]] = list("Population" = population, "NumCells" = num_cells, "Genotype" = genotype, "GenotypeTrue" = genotype_true)
    
  }
  
  return(scDNASimulation)

}

scDNAMixture_v1 = function(simulation){
  
  NOC = length(simulation)
  
  num_cells_simulated = array(NA, dim = NOC)
  num_cells = array(NA, dim = NOC)
  
  for (i in 1:NOC){
    
    num_cells_simulated[i] = simulation[[i]]$NumCells
    
    num_cells[i] = as.numeric(readline(prompt = paste("The number of cells in the mixture from contributor ", i, ": ", sep = "")))
    if (is.na(num_cells[i]) | num_cells[i] > num_cells_simulated[i]){
      num_cells[i] = num_cells_simulated[i] # If the user did not enter or the number is greater than the number of cells simulated, number of cells is by default to be the number of cells simulated.
    } else if (num_cells[i] < 0 | num_cells[i] != as.integer(num_cells[i])){
      stop("The number of cells has to be a non-negative integer.")
    }
  }
  
  num_cells_mixture = sum(num_cells)
  num_cells_mixed = 0
  
  scDNAMixture = list("NOCSimulated" = NOC, "NOC" = sum(num_cells != 0), "NumCells" = num_cells)
  
  genotype_mixture = array(NA, dim = c(num_cells_mixture, 20, 2))
  dimnames(genotype_mixture) = list(paste0("Cell_", 1:num_cells_mixture), rownames(simulation[[1]]$GenotypeTrue), c("Allele_1", "Allele_2"))
  
  for (i in 1:NOC){
    if (num_cells[i] > 0){
      ind_cells_mixture = sample(num_cells_simulated[i], size = num_cells[i], replace = FALSE)
      genotype = simulation[[i]]$Genotype[ind_cells_mixture, , ]
      contributor_name = paste("Contributor_", i, sep = "")
      scDNAMixture[[contributor_name]] = genotype
      genotype_mixture[(num_cells_mixed + 1):(num_cells_mixed + num_cells[i]), , ] = genotype
      num_cells_mixed = num_cells_mixed + num_cells[i]
    }
  }
  
  scDNAMixture[["Mixture"]] = genotype_mixture
  
  return(scDNAMixture)
  
}

simulation = scDNASimulation_v1()
# The result will be a list of contributors, each with population, The number of cells simulated, and the corresponding genotype (2 alleles) at the 20 core CODIS loci, and the true genotype.
# simulation[[n]] gives the list for the n-th contributor.
# simulation[[n]]$Population, simulation[[n]]$NumCells, and simulation[[n]]$Genotype gives the corresponding results.
# Genotype is a three-dimentional array, with the first dimension being the cells (NumCells in total), the second dimension being the loci (20 in total), and the third dimension being the pair of alleles.
# simulation[[2]]$Genotype[3, 5, ] gives the genotype of the fourth locus (D13S317) of the third cell simulated of the second contributor.

mixture = scDNAMixture_v1(simulation)
# The result will be a list containing information such as NOC, how many cells in the mixture belongs to each person, their selected cells, and the mixture.
# In the mixture, the rows are ordered from contributor 1 to the last contributor.










########################################
########################################
########################################

library(ggplot2)
library(dplyr)
install.packages("patchwork")
library(patchwork)
library(ggplotify)
install.packages("igraph")
library(igraph)

# ggplot(AlleleFreq_AfAm, aes(x = Alleles, y = CSF1PO)) + geom_bar(stat = "identity")
# barplot(AlleleFreq_AfAm[,2])
# plot(AlleleFreq_AfAm[,1],AlleleFreq_AfAm[,2],type="b")
# plot(AlleleFreq_AfAm[,1],AlleleFreq_AfAm[,2],type="l")
# lines(AlleleFreq_Cauc[,1],AlleleFreq_Cauc[,2],col=2)
# lines(AlleleFreq_Hispanic[,1],AlleleFreq_Hispanic[,2],col=3)
# lines(AlleleFreq_Asian[,1],AlleleFreq_Asian[,2],col=4)

# Plot.allele = function(i){
#   par(mfrow=c(1,2))
#   plot(AlleleFreq_AfAm[,1],AlleleFreq_AfAm[,i+1],type="l", xlab = "Alleles", ylab = "Probability", main = colnames(AlleleFreq_AfAm)[i+1])
#   lines(AlleleFreq_Cauc[,1],AlleleFreq_Cauc[,i+1],col=2)
#   lines(AlleleFreq_Hispanic[,1],AlleleFreq_Hispanic[,i+1],col=3)
#   lines(AlleleFreq_Asian[,1],AlleleFreq_Asian[,i+1],col=4)
# }

# Plot.allele(2)

PlotAlleleFreq = function(i){
  
  xmin = floor(min(AlleleFreq_all[AlleleFreq_all[, i + 1] != 0, "Alleles"])) - 1
  xmax = ceiling(max(AlleleFreq_all[AlleleFreq_all[, i + 1] != 0, "Alleles"])) + 1
  
  plot1 = ggplot() +
    xlim(xmin, xmax) +
    geom_line(data = AlleleFreq_AfAm, aes(x = Alleles, y = AlleleFreq_AfAm[, i + 1], color = "African American")) +
    geom_point(data = AlleleFreq_AfAm, aes(x = Alleles, y = AlleleFreq_AfAm[, i + 1], color = "African American")) +
    geom_line(data = AlleleFreq_Cauc, aes(x = Alleles, y = AlleleFreq_Cauc[, i + 1], color = "Caucasian")) +
    geom_point(data = AlleleFreq_Cauc, aes(x = Alleles, y = AlleleFreq_Cauc[, i + 1], color = "Caucasian")) +
    geom_line(data = AlleleFreq_Hispanic, aes(x = Alleles, y = AlleleFreq_Hispanic[, i + 1], color = "Hispanic")) +
    geom_point(data = AlleleFreq_Hispanic, aes(x = Alleles, y = AlleleFreq_Hispanic[, i + 1], color = "Hispanic")) +
    geom_line(data = AlleleFreq_Asian, aes(x = Alleles, y = AlleleFreq_Asian[, i + 1], color = "Asian")) +
    geom_point(data = AlleleFreq_Asian, aes(x = Alleles, y = AlleleFreq_Asian[, i + 1], color = "Asian")) +
    scale_color_manual(values = c("African American" = "black", 
                                  "Caucasian" = "blue", 
                                  "Hispanic" = "red", 
                                  "Asian" = "orange")) +
    labs(title = colnames(AlleleFreq_all)[i+1], x = "Alleles", y = "Probability", color = "Population") +
    theme_minimal() +
    theme(legend.position = c(0.9, 0.8))
  
  plot2 = function(){
    distance_matrix = as.matrix(dist(scale(rbind(AlleleFreq_AfAm[,i+1], AlleleFreq_Cauc[,i+1], AlleleFreq_Hispanic[,i+1], AlleleFreq_Asian[,i+1])), method = "euclidean"))
    distance_matrix = 100/distance_matrix
    diag(distance_matrix) = 0
    rownames(distance_matrix) = c("AfAm", "Cauc", "Hispanic", "Asian")
    colnames(distance_matrix) = c("AfAm", "Cauc", "Hispanic", "Asian")
    g = graph.adjacency(distance_matrix, mode = "undirected", weighted = TRUE)
  
    V(g)$color = c("yellow", "yellow", "yellow", "yellow")
    E(g)$color = "gray"
  
    plot(g, edge.width = E(g)$weight, vertex.size = 30, vertex.label.cex = 1.5)
  }
  
  plot2 = as.ggplot(plot2)
  
  print(plot1+plot2)

}

dev.new()
par(mfrow = c(20, 1))
for (i in 1:20){
  PlotAlleleFreq(i)
}

dev.new()
par(mfrow=c(4,5),mar=c(.5,.5,.5,.5))
for (i in 1:20){
  distance_matrix = as.matrix(dist(scale(rbind(AlleleFreq_AfAm[,i+1], AlleleFreq_Cauc[,i+1], AlleleFreq_Hispanic[,i+1], AlleleFreq_Asian[,i+1])), method = "euclidean"))
  distance_matrix = 100/distance_matrix
  diag(distance_matrix) = 0
  rownames(distance_matrix) = c("AfAm", "Cauc", "Hispanic", "Asian")
  colnames(distance_matrix) = c("AfAm", "Cauc", "Hispanic", "Asian")
  g = graph.adjacency(distance_matrix, mode = "undirected", weighted = TRUE)
  
  V(g)$color = c("yellow", "yellow", "yellow", "yellow")
  E(g)$color = "gray"
  
  plot(g, edge.width = 5 * (E(g)$weight/max(E(g)$weight))^2, vertex.size = 30, vertex.label.cex = 1.5, main = colnames(AlleleFreq_all)[i+1])
}
