#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")

setwd("~/Git-Projects/Git-Research-Projects/CNprep-Slicing-CORE-Analysis/")
source("analysis/genomicFeatureAssignment.R")

#
# Load sample to retrieve feature set for
# TODO: There seems to be a scope conflict - samples is getting overwritten
#
samples <- load_samples(classes = c("T","F", "M"), sampleList = "./resources/sampleList.csv")

#
# Retrieve CORE features
#
Acores <- retrieveCores("./output/coresResults/prev_run_8_2_2018_2/selectedCores/AselectedCoresBP.bed") # BED file of amplification recurrent regions
Dcores <- retrieveCores("./output/coresResults/prev_run_8_2_2018_2/selectedCores/DselectedCoresBP.bed") # BED file of deletion recurrent regions
ADcores <- retrieveCores("./output/coresResults/prev_run_8_2_2018_2/selectedCores/ADselectedCoresBP.bed") # BED file of both recurrent regions

aucData <- readRDS("./resources/listSampleTESAUC.RDS")



#
# TODO: NEED TO UPDATE RETRIEVE TRAINING SET - SPECIFICALLY HOW IT DEFINES THE SCORES
#
training_set <- retrieveTrainingSet(samples, Acores, Dcores, organoidSlicesFile = "./resources/slicingOutput/table/prev_run_8_2_2018_3/organoidSlices.txt")
training_set$matrix <- attachLabelsToSet(matrix_training_set = training_set$matrix, labelData = aucData)

visualizeUnclusteredHeatmap(training_set$melted)
# TODO: Deal with 0 variance issue
hc <- clusterTrainingSet(training_set$melted, visualize = TRUE)
plot(hc)

cd_local("mlOutput")
write.csv(training_set$matrix, file ="coreTrainingSet_8_3_2018_1.csv")
