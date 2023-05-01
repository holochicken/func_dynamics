# Examining MCMC convergence - updated 22/03/2023

# Packages
library(Hmsc)
library(colorspace)
library(vioplot)


# 1- Load fitted model ---------------------
# Directory
localDir = "."
ModelDir = file.path(localDir, "src/hmsc_mg_and_mt/mg/models")
PanelDir = file.path(localDir, "src/hmsc_mg_and_mt/mg/panels")

# Chain characteristics
samples_list = c(5,250,250)
thin_list = c(1,1,10)
nst = length(thin_list)
nChains = 4

# 2- Evaluate convergence ---------------------
set.seed(1)
for (i in 1) {
  ma = NULL
  na = NULL
  ma_rho = NULL
  na_rho = NULL
  for (Lst in 1:nst) {
    thin = thin_list[Lst]
    samples = samples_list[Lst]
    filename = paste("model","_thin_", as.character(thin),
                     "_samples_", as.character(samples),
                     "_chains_", as.character(nChains),
                     ".Rdata", sep = "")

    load(file = file.path(ModelDir, filename))
    mpost = convertToCodaObject(m,
                                spNamesNumbers = c(T,F),
                                covNamesNumbers = c(T,F))
    ## beta - Species niches - response of MAGs to the fixed effects
    psrf.beta = gelman.diag(mpost$Beta,multivariate = FALSE)$psrf
    ## rho - Influence of phylogeny on species niches - phylogenetic signal
    psrf.rho = gelman.diag(mpost$Rho,multivariate = FALSE)$psrf
    ## beta
    if (is.null(ma)) {
      ma = psrf.beta[,1]
      na = paste0("model_", as.character(thin),
                  ",", as.character(samples))
    } else {
      ma = cbind(ma,psrf.beta[,1])
      na = c(na,paste0("model_", as.character(thin),
                       ",", as.character(samples)))
    }
    ## rho
    if (is.null(ma_rho)) {
      ma_rho = psrf.rho[,1]
      na_rho = paste0("model_", as.character(thin),
                      ",", as.character(samples))
    } else {
      ma_rho = cbind(ma_rho,psrf.rho[,1])
      na_rho = c(na_rho,paste0("model_", as.character(thin),
                               ",", as.character(samples)))
    }
  }
  # create file name
  panel.name = paste("MCMC_convergence",
                     "_model_", as.character(i),
                     ".pdf", sep = "")
  # Plot
  pdf(file = file.path(PanelDir,panel.name))
  ## beta
  par(mfrow = c(2,1))
  vioplot(ma,
          col = rainbow_hcl(1),
          names = na,
          ylim = c(0, max(ma)),
          main = "psrf(beta)")
  vioplot(ma,
          col = rainbow_hcl(1),
          names = na,
          ylim = c(0.9,1.1),
          main = "psrf(beta)")

  ## rho
  par(mfrow = c(2,1))
  vioplot(ma_rho,
          col = rainbow_hcl(1),
          names = na_rho,
          ylim = c(0, max(ma_rho)),
          main = "psrf(rho)")
  vioplot(ma_rho,
          col = rainbow_hcl(1),
          names = na_rho,
          ylim = c(0.9,1.1),
          main = "psrf(rho)")
  dev.off()
}

