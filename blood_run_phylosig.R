library(phytools)
library(ape)

# load data
# setwd("~/Dropbox/Work/Research/Projects/trainee_projects/Adeyemi/comparative_epi_aging/")
library(phytools)
library(ape)
blood_slope_dat <- readRDS("~/Downloads/blood_regression_slope_wide.rds")
blood_se_dat <- readRDS("~/Downloads/blood_regression_std_error_wide.rds")
blood_tree <- read.tree("~/Downloads/blood_species.nwk")

# define function that applies phylosig to a single row of data
run_phylosig <- function(i, tree, slope_dat, se_dat) {
  if (i %% 500 == 0) message("Processing CpG site ", i, " of ", nrow(slope_dat))
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
  x <- phylosig(tree, x = slopes, se = std_errs, method="lambda", test = TRUE)
  output <- as.data.frame(matrix(nrow=1, ncol=3))
  output[1,] <- c(x$lambda, x$sig2, x$P)
  colnames(output) <- c("lambda", "sig2", "p_val")
  return(output)
}

# use sapply() to run the function for each row of data
blood_phylo_results <- sapply(1:nrow(slope_dat), run_phylosig,
                  tree = blood_tree,
                  slope_dat = blood_slope_dat,
                  se_dat = blood_se_dat)

blood_phylo_results <- as.data.frame(t(blood_phylo_results))
rownames(results) <- rownames(blood_slope_dat)

# look into parSapply() from the R package 'parallel' to parallelize above operation




run_phylosig <- function(i, tree, slope_dat, se_dat) {
  if (i %% 500 == 0) cat("Processing CpG site", i, "of", nrow(slope_dat), "\n")
  
  slopes <- as.numeric(slope_dat[i, ])
  std_errs <- as.numeric(se_dat[i, ])
  
  names(slopes) <- colnames(slope_dat)
  names(std_errs) <- colnames(se_dat)
  
  slopes <- slopes[tree$tip.label]
  std_errs <- std_errs[tree$tip.label]
  
  # Try-catch block
  out <- tryCatch({
    x <- phylosig(tree, x = slopes, se = std_errs, method = "lambda", test = TRUE)
    c(lambda = x$lambda,
      sig2   = x$sig2,
      p_val  = if (!is.null(x$P)) x$P else NA_real_)
  }, error = function(e) {
    c(lambda = NA_real_, sig2 = NA_real_, p_val = NA_real_)
  })
  
  return(out)
}

# use sapply() to run the function for each row of data
blood_phylo_results <- sapply(1:nrow(slope_dat), run_phylosig,
                              tree = blood_tree,
                              slope_dat = blood_slope_dat,
                              se_dat = blood_se_dat)

blood_phylo_results <- as.data.frame(t(blood_phylo_results))
rownames(blood_phylo_results) <- rownames(blood_slope_dat)


saveRDS(blood_phylo_results, file = "~/Downloads/blood_phylo_result.rds")


