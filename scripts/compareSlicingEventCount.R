setwd("~/Git-Projects/Git-Research-Projects/CNprep-Slicing-CORE-Analysis") # Set working directory to where the scripts are located at
source("scripts/coreGenerationLibrary.R")
source("scripts/helperFunctions.R")

slicing_dir_hN31 <- "./resources/slicingOutput/prev_run_7_28_18_5/"
slicing_dir_NA12878 <- "./resources/slicingOutput/prev_run1/"

samples <- load_samples(classes = c("T"), sampleList = "./resources/sampleList.csv")
events <- c("A", "D")

NA12878_X <- list()
hN31_Y <- list()
list.index <- 1
for (sample in samples){
  print(sample)
  try({
    inputCORESegments_X <- loadSlicingRegions(dir =  slicing_dir_NA12878, samples = sample, events = events, probes = FALSE, silent = TRUE)  
    inputCORESegments_Y <- loadSlicingRegions(dir =  slicing_dir_hN31, samples = sample, events = events, probes = FALSE, silent = TRUE)  
    NA12878_X[list.index] <- nrow(inputCORESegments_X)
    hN31_Y[list.index] <- nrow(inputCORESegments_Y)
    list.index <- list.index + 1
  }, silent = TRUE)
}
NA12878_X <- unlist(NA12878_X)
hN31_Y <- unlist(hN31_Y)

plot(NA12878_X, hN31_Y, xlim = c(5, 50), ylim = c(5, 50))
abline(lm(I(hN31_Y) ~ 0 + NA12878_X), col = "red")
abline(1,1, col = "blue")
dev.off()


