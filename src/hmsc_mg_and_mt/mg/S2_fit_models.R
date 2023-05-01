# Fitting the model - updated 22/03/2023

# Packages
library(Hmsc)

# 1- Load unfitted model  ---------------------
# Directory
localDir = "."
ModelDir = file.path(localDir, "src/hmsc_mg_and_mt/mg/models")

# Chain characteristics
samples_list = c(5, 250, 250)
thin_list = c(1,1,10)
nChains = 4
load(file.path(ModelDir,"Unfitted_model.Rdata"))

unf_model <- m


# 2- Fit model ---------------------
set.seed(1)
for (Lst in 1:length(samples_list)) {
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  # fit the model
  m = sampleMcmc(unf_model,
                 samples = samples,
                 thin = thin,
                 adaptNf = rep(ceiling(0.4*samples*thin), unf_model$nr),
                 transient = ceiling(0.5*samples*thin),
                 nChains = nChains,
                 nParallel = nChains)
  # create file name
  filename = paste("model",
                   "_thin_", as.character(thin),
                   "_samples_", as.character(samples),
                   "_chains_", as.character(nChains),
                   ".Rdata", sep = "")
# save file
  save(m, file = file.path(ModelDir,filename))
}
