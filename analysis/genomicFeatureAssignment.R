#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")

source("scripts/helperFunctions.R")
source("analysis/incidence.R")
library(GenomicRanges)
library(ggplot2) 
library(reshape) 
library(reshape2)
library(made4)
library(cluster)
library(spatstat) # "im" function 

retrieveCores <- function(dir){
  return(read.table(dir, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
}

retrieveOrganoidSlices <- function(dir){
  return(read.table(dir, header = TRUE, sep = "\t", stringsAsFactors = FALSE))
}

# TODO: NEED TO UPDATE RETRIEVE TRAINING SET - SPECIFICALLY HOW IT DEFINES THE SCORES
retrieveTrainingSet <- function(loaded_samples, Acores, Dcores, ADcores, organoidSlicesFile){
  melted_training_set <- data.frame(stringsAsFactors = FALSE)  
  matrix_training_set <- data.frame(stringsAsFactors = FALSE)
  
  organoidSlices <- retrieveOrganoidSlices(organoidSlicesFile)
  for(sample in loaded_samples){
    #
    # Retrieve necessary data for feature value calculation (bins with cnlr and COREs)
    #
    final_incidence_table <- data.frame(sample= sample)
    if(!missing(ADcores)){
      print("Adding ADcores")
      incidence.input.core <- ADcores[,c(6,7)] # TODO: CORE should have headers so that I can call based on colname instead of index
      colnames(incidence.input.core) <- c("start", "end")
      incidence.input.events <- organoidSlices[organoidSlices$profID == sample, c("gstart", "gend")] 
      colnames(incidence.input.events) <- c("start", "end")
      incidence.output.table <- incidence(incidence.input.core, incidence.input.events, dropevents="Greedy",assoc="I")
      final_incidence_table <- cbind(final_incidence_table, t(as.data.frame(incidence.output.table)))
    } else {
      if(!missing(Acores)){
        print("Adding Acores")
        incidence.input.core <- Acores[,c(6,7)] # TODO: CORE should have headers so that I can call based on colname instead of index
        colnames(incidence.input.core) <- c("start", "end")
        incidence.input.events <- organoidSlices[organoidSlices$profID == sample & organoidSlices$lesion.type == 1, c("gstart", "gend")] 
        colnames(incidence.input.events) <- c("start", "end")
        incidence.output.table <- incidence(incidence.input.core, incidence.input.events, dropevents="Greedy",assoc="I")
        final_incidence_table <- cbind(final_incidence_table, t(as.data.frame(incidence.output.table)))
      }
      if(!missing(Dcores)){
        print("Adding Dcores")
        incidence.input.core <- Dcores[,c(6,7)] # TODO: CORE should have headers so that I can call based on colname instead of index
        colnames(incidence.input.core) <- c("start", "end")
        incidence.input.events <- organoidSlices[organoidSlices$profID == sample & organoidSlices$lesion.type == 0, c("gstart", "gend")] 
        colnames(incidence.input.events) <- c("start", "end")
        incidence.output.table <- incidence(incidence.input.core, incidence.input.events, dropevents="Greedy",assoc="I")
        final_incidence_table <- cbind(final_incidence_table, t(as.data.frame(incidence.output.table)))
      } 
    }
    rownames(final_incidence_table) <- sample
    final_incidence_table <- final_incidence_table[, -c(1)]
    matrix_training_set <- rbind(matrix_training_set, final_incidence_table)
  }
  
  #
  # Convert melt the matrix_training_set
  #
  melted_training_set <- do.call(rbind, lapply(seq(nrow(matrix_training_set)), function(index){
    return(do.call(rbind, lapply(colnames(matrix_training_set[index, ]), function(coreId, index){
      coreEntry <- data.frame(score = matrix_training_set[index, coreId], coreId = coreId, sampleId = rownames(matrix_training_set[index, ])[1])
      return(coreEntry)
    }, index)))
  }))
  
  return(list(melted=melted_training_set, matrix=matrix_training_set))
}

attachLabelsToSet <- function(matrix_training_set, labelData){
  sampleList <- rownames(matrix_training_set)
  labelLists <- lapply(names(labelData), function(label){
    aucList <- unlist(sapply(sampleList, function(sample, label){
      labelMatrix <- labelData[[label]]  
      auc <- c(labelMatrix[labelMatrix$SampleId == sample, ]$AUC, NA)[1]
      return(auc)
    }, label))
    return(aucList)
  })
  names(labelLists) <- names(labelData)
  labelDataframe <- do.call(cbind.data.frame, labelLists)
  labeled_matrix_training_set <- cbind(labelDataframe, matrix_training_set)
  return(labeled_matrix_training_set)
}

visualizeUnclusteredHeatmap <- function(training_set){
  ggplot(data = training_set, aes(x = coreId, y = sampleId)) + 
    geom_tile(aes(fill = score), color = "white", size = 1) + 
    scale_fill_gradient2(low = "blue", mid="white", high = "tomato") + 
    xlab("core ID") + 
    theme_grey(base_size = 10) + 
    ggtitle("Heatmap (ggplot)") + 
    theme(axis.ticks = element_blank(), 
          panel.background = element_blank(), 
          plot.title = element_text(size = 12, colour = "gray50")) 
}

clusterTrainingSet <- function(training_set, visualize = FALSE){
  # Unmelt training set for correlation analysis
  training_set_matrix <- dcast(data = training_set,formula = sampleId~coreId,fun.aggregate = sum,value.var = "score")
  sampleIds <- training_set_matrix$sampleId
  training_set_matrix <- training_set_matrix[,-c(1)]
  training_set_matrix <- t(training_set_matrix)
  colnames(training_set_matrix) <- sampleIds
  training_set_matrix <- as.data.frame(training_set_matrix)
  
  # Remove samples with 0 variance
  nonzero_variance_samples <- unlist(lapply(colnames(training_set_matrix), function(sample){
    if(var(training_set_matrix[, sample]) != 0){
      return(sample)
    }
  }))
  training_set_matrix <- training_set_matrix[,c(nonzero_variance_samples)]
  
  # Calculate distance matrix
  corRaw <- cor(training_set_matrix)
  dissimilarity <- 1 - corRaw
  distance.sample <- as.dist(dissimilarity)
  distance.core <- as.dist(t(dissimilarity))
  
  # Run hierarchical clustering
  hc.sample <- hclust(distance.sample)
  hc.core <- hclust(distance.core)
  
  if(visualize == TRUE){
    color.palette  <- colorRampPalette(c("blue", "white", "tomato"))(n=600)
    heatmap.2(t(training_set_matrix),  trace="none", dendrogram="row", density.info = 'none', scale='none', col = color.palette)
  }
  
  return(hc.sample)  
}