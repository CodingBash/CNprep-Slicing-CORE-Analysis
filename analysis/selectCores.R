#
# Quick interactive script to sub-select cores based on p-value
# TODO: Make into callable script with arguments
#

#
# Load source libraries
#
setwd("~/Git-Projects/Git-Research-Projects/CNprep-Slicing-CORE-Analysis")
source("scripts/helperFunctions.R")

#
# Set the core p_threshold
#
p_threshold <- 0.002

rdsFile <- "./output/coresResults/ADnewCOREobjBP_slicing2.rds"
tableFile <- "./output/coresResults/ADcoreTableBP_slicing2.csv"
outputBed <- "./output/selectedCores/ADcoresBP.bed"

#
# Get core information
#
obj <- readRDS(rdsFile) # Retrieve CORE object
coreTable <- read.table(tableFile, header = TRUE, sep = ",") # Retrieve CORE table CSV (since the table's scale is chromosome-based instead of absolute)

#
# Subset core information using p_threshold then write to file
#
coreTable <- coreTable[which(obj$p<p_threshold),] # Filter cores based on p-value threshold
coreTable <- coreTable[,c(2,3,4)] # Now convert to BED format
write.table(coreTable, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE, file = outputBed) # Write to bed file
