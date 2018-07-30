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
Acores <- retrieveCores("./output/selectedCores/AcoresBP.bed") # BED file of amplification recurrent regions
Dcores <- retrieveCores("./output/selectedCores/DcoresBP.bed") # BED file of deletion recurrent regions
ADcores <- retrieveCores("./output/selectedCores/ADcoresBP.bed") # BED file of both recurrent regions

aucData <- readRDS("./resources/listSampleTESAUC.RDS")


#
# Retrieve training set
# TODO: HERE
setwd("~/Git-Projects/Git-Research-Projects/FACETS_write_files/")

reference <- "hN31"
sample_dir <- "./output/FACETS_Reference_hN31_7_28_18_2/"

#
# TODO: NEED TO UPDATE RETRIEVE TRAINING SET - SPECIFICALLY HOW IT DEFINES THE SCORES
#
training_set <- retrieveTrainingSet(samples, Acores, Dcores, sample_subdir = "/",  reference = reference, binDir = sample_dir)
training_set$matrix <- attachLabelsToSet(matrix_training_set = training_set$matrix, labelData = aucData)

visualizeUnclusteredHeatmap(training_set$melted)
hc <- clusterTrainingSet(training_set$melted, visualize = TRUE)
plot(hc)

cd_local("mlOutput")
write.csv(training_set$matrix, file ="coreTrainingSet_7_26_2018_1.csv")
