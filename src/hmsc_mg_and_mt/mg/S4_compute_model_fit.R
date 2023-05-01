# Computing model fit - updated 22/03/2023

# Packages
library(Hmsc)


# 1- Open fitted model ---------------------
# Directory
localDir = "."
ModelDir = file.path(localDir, "src/hmsc_mg_and_mt/mg/models")
MFDir = file.path(localDir, "src/hmsc_mg_and_mt/mg/model_fit")

thin = 10
samples = 250
nChains = 4
filename = paste("model","_thin_", as.character(thin),
                 "_samples_", as.character(samples),
                 "_chains_",as.character(nChains),".Rdata",sep = "")
load(file = file.path(ModelDir,filename))

# Create partitions
set.seed(12)
partition_1 <- createPartition(hM = m, nfolds = 2)
partition_2 <- c(rep(1,sum(m$XData$trial == "CA")),
                 rep(2,sum(m$XData$trial == "CB")))


# 2- Model fit ---------------------
set.seed(1)
predY.MF <- computePredictedValues(m, expected = TRUE)
MF <- evaluateModelFit(hM = m, predY = predY.MF)
filename <- file.path(MFDir, paste("MF_",
                                  "chains_",as.character(nChains),
                                  "_thin_", as.character(thin),
                                  "_samples_", as.character(samples),
                                  ".rds", sep = ""))
saveRDS(MF,file = filename)
MF <- readRDS(file = filename)
mean(MF$R2)


# Model fit using 2-fold CV: samples randomly assigned to folds
set.seed(1)
predY.CV_1 <- computePredictedValues(m,
                                     expected = TRUE,
                                     partition = partition_1,
                                     nChains = 1,
                                     nParallel = 1)
MF.CV_1 <- evaluateModelFit(hM = m, predY = predY.CV_1)
# create file name and save
filename <- file.path(MFDir, paste("MF_CV1_",
                                  "chains_",as.character(nChains),
                                  "_thin_", as.character(thin),
                                  "_samples_", as.character(samples),
                                  ".rds",sep = ""))
saveRDS(MF.CV_1,file = filename)
MF_CV_1 <- readRDS(file = filename)
mean(MF_CV_1$R2)

# Model fit using 2-fold CV: each fold contains the samples from one trial
set.seed(1)
predY.CV_2 <- computePredictedValues(m,
                                     expected = TRUE,
                                     partition = partition_2,
                                     nChains = 1,
                                     nParallel = 1)
MF.CV_2 <- evaluateModelFit(hM = m, predY = predY.CV_2)
# create file name and save
filename <- file.path(MFDir, paste("MF_CV2_",
                                  "chains_",as.character(nChains),
                                  "_thin_", as.character(thin),
                                  "_samples_", as.character(samples),
                                  ".rds", sep = ""))
saveRDS(MF.CV_2, file = filename)
MF_CV_2 <- readRDS(file = filename)
mean(MF_CV_2$R2)
