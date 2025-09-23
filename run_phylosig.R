library(phytools)
library(ape)

# load data
setwd("~/Dropbox/Work/Research/Projects/trainee_projects/Adeyemi/comparative_epi_aging/")

tree <- read.newick("skin_species.nwk")
slope_dat <- readRDS("skin_rates_50.rds")
se_dat <- readRDS("skin_std_err_50.rds")

# define function that applies phylosig to a single row of data
run_phylosig <- function(i, tree, slope_dat, se_dat) {
  # convert data frame row to numeric vector
  slopes <- as.numeric(slope_dat[i,])
  std_errs <- as.numeric(se_dat[i,])
  # add names back
  names(slopes) <- colnames(slope_dat)
  names(std_errs) <- colnames(se_dat)
  # re-order vectors
  slopes <- slopes[tree$tip.label]
  std_errs <- std_errs[tree$tip.label]
  # run phylosig
  x <- phylosig(tree, x = slopes, se = std_errs, method="lambda")
  output <- as.data.frame(matrix(nrow=1, ncol=2))
  output[1,] <- c(x$lambda, x$sig2)
  colnames(output) <- c("lambda", "sig2")
  return(output)
}

# use sapply() to run the function for each row of data
results <- sapply(1:nrow(slope_dat), run_phylosig,
                  tree = tree,
                  slope_dat = slope_dat,
                  se_dat = se_dat)

results <- as.data.frame(t(results))
rownames(results) <- rownames(slope_dat)

# look into parSapply() from the R package 'parallel' to parallelize above operation
