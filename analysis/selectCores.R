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
p_threshold <- 0.05

events <- c("A","D","AD")
output <- "output/coresResults"
dir <- "prev_run_8_2_2018_2"
#dir <- "prev_run_7_28_2018_x_1"
for(event in events){
  rdsFile <- paste0("./", output, "/", dir, "/", event, "newCOREobjBP_slicing2.rds")
  tableFile <- paste0("./", output, "/", dir, "/", event, "coreTableBP_slicing2.csv")
  outputBed <- paste0("./", output, "/", dir, "/selectedCores/", event, "selectedCoresBP.bed")
  
  #
  # Get core information
  #
  obj <- readRDS(rdsFile) # Retrieve CORE object
  coreTable <- read.table(tableFile, header = TRUE, sep = ",") # Retrieve CORE table CSV (since the table's scale is chromosome-based instead of absolute)
  
  #
  # Subset core information using p_threshold then write to file
  #
  coreTable <- coreTable[which(obj$p<p_threshold),] # Filter cores based on p-value threshold
  coreTable <- coreTable[,c(2,3,4,5)] # Now convert to BED format
  coreTable$p <- obj$p[unlist(obj$p < p_threshold)]
  
  # TODO: reposition this before score and p-value
  coreTable$gstart <- obj$coreTable[,"start"][which(obj$p<p_threshold)]
  coreTable$gend <- obj$coreTable[,"end"][which(obj$p<p_threshold)]
  
  # TODO: Keep headers in core table - need to udpate downstream CORE analysis then
  write.table(coreTable, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE, file = outputBed) # Write to bed file
}