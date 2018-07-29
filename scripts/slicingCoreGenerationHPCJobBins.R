#
# This script generates CORES from input files in BP units (instead of SNP/bin units).
# This script is not interactive and is meant to be ran on a HPC job from a shell script.
#

args <- commandArgs(trailingOnly = TRUE)

#
# Set script arguments
#
event <- "A" 
outputCsv <- "coreTable.csv"
outputObj <- "newCOREobj.rds"
if (length(args) == 1){
	event <- args[1]
} else if (length(args) == 2){
	event <- args[1]
	outputCsv <- args[2]	
} else if (length(args) == 3){
	event <- args[1]
	outputCsv <- args[2]
	outputObj <- args[3]
}

setwd("~/code/slicing_cores/scripts")
source("coreGenerationLibrary.R")
source("helperFunctions.R")

events <- c(event)
slicing_dir <- "./resources/slicingOutput/prev_run_7_28_18_5/"
#
# Get CORE input
#
cd_slicing_core()
samples <- load_samples(classes = c("T", "M", "F"), sampleList = "./resources/sampleList.csv")
chromosomeSizes <- readRDS("./resources/chromosomeSizes.rds")
inputCORESegments <- loadSlicingRegions(dir =  slicing_dir, samples = samples, events = events, probes = FALSE, silent = TRUE)
print("Prepared all inputs - now running CORE")

#
# Run CORE
#
outputCOREobj <- runCORE(inputCORESegments, distrib="Grid", maxmark=500, nshuffle=500, seedme = 1234, njobs=4)
print("CORE run complete")

#
# Save CORE object
# WARNING: The coreTable in the CORE object may not be in chromosomal location units. See CORE table output
#
saveRDS(outputCOREobj, outputObj)
print(paste("Saved CORE obj to", outputObj, sep = ""))


#
# Save CORE table
#
COREtable <- retrieveCORETable(outputCOREobj, rescaleOutput = FALSE)
write.csv(COREtable, outputCsv)
print(paste("Saved coreTable as", outputCsv, sep = ""))

print("script complete")
