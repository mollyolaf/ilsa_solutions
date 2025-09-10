library(lsasim)
library(MASS)
library(matrixcalc)
library(fastDummies)
library(performance)

###############################
# 1. Setting up constants and parameters
# (These remain the same across all replications)
###############################

# generating questionnaires and thetas
# proportions for ordinal variables (continuous are "1")
# 3 theta, 1 continuous for school, 3 continuous, 1 binary, 4 ordinal
props <- list(1, 1, 1, 
              1, #q3 is school
              1, 1, 1, #q4 to q6 are continuous
              c(.50, 1), # q7 binary variable will be used as direct covariate. 
              # correlations vary with other covariates and thetas.
              c(.20, .50, .80, 1), 
              c(.10, .40, .60, 1), 
              c(.20, .70, .90, 1), 
              c(.20, .50, .80, 1)
)

# generating covariance matrix between theta and bq and among bq
# first 3 variables will be used as thetas for 3 arbitrary domains
yw_cov <- matrix(c(
  1.00, 0.80, 0.85, 0.40, 0.25, 0.30, 0.10, 0.35, 0.20, 0.45, 0.35, 0.15,
  0.80, 1.00, 0.79, 0.40, 0.20, 0.25, 0.15, 0.25, 0.20, 0.30, 0.40, 0.05,
  0.85, 0.79, 1.00, 0.40, 0.30, 0.35, 0.20, 0.30, 0.25, 0.15, 0.35, 0.10,
  0.40, 0.40, 0.40, 1.00, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20,
  0.25, 0.20, 0.30, 0.20, 1.00, 0.45, 0.30, 0.40, 0.35, 0.05, 0.20, 0.30,
  0.30, 0.25, 0.35, 0.20, 0.45, 1.00, 0.25, 0.20, 0.15, 0.40, 0.25, 0.25,
  0.10, 0.15, 0.20, 0.20, 0.30, 0.25, 1.00, 0.35, 0.30, 0.10, 0.15, 0.20,
  0.35, 0.25, 0.30, 0.20, 0.40, 0.20, 0.35, 1.00, 0.40, 0.35, 0.20, 0.15,
  0.20, 0.20, 0.25, 0.20, 0.35, 0.15, 0.30, 0.40, 1.00, 0.20, 0.10, 0.25,
  0.45, 0.30, 0.15, 0.20, 0.05, 0.40, 0.10, 0.35, 0.20, 1.00, 0.40, 0.30,
  0.35, 0.40, 0.35, 0.20, 0.20, 0.25, 0.15, 0.20, 0.10, 0.40, 1.00, 0.25,
  0.15, 0.05, 0.10, 0.20, 0.30, 0.25, 0.20, 0.15, 0.25, 0.30, 0.25, 1.00
), nrow = 12, byrow = TRUE)

# check if positive definite
if(is.positive.definite(yw_cov)) {
  cat("The generated matrix is positive definite.\n")
} else {
  cat("An error occurred. The matrix is not positive definite.\n")
}

# number of  items (PISA 2022 technical report, chapter 2, p. 3)
Im <- 234 # math
Ir <- 197 # reading
Is <- 115 # science

n <- 154*28 # 154 schools with 28 students each
n_reps <- 100 # number of replications

# Define a function to calculate means and variances thetas
calc_m_v <- function(data, variables) {
  # Use sapply to calculate both the mean and variance for each variable passed in the list
  stats <- sapply(variables, function(var) {
    mean_val <- mean(data[[var]], na.rm = TRUE)
    var_val <- var(data[[var]], na.rm = TRUE)
    return(c(mean = mean_val, variance = var_val))
  })
  # Return the statistics as a named matrix
  return(stats)
}

# setting path
path <- "<your path here>"

# Create directory for replications if it doesn't exist
if (!dir.exists(paste0(path, "data/replications/"))) {
  dir.create(paste0(path, "data/replications/"), recursive = TRUE)
}

###############################
# 2. Main simulation loop
###############################

for (rep in 1:n_reps) {
  cat(sprintf("Starting replication %d of %d\n", rep, n_reps))
  
  # Set seed for reproducibility (different seed for each replication)
  set.seed(178150 + rep)
  
  ###############################
  # 3. Simulating multivariate normal thetas for 3 arbitrary domains
  ###############################
  
  bq <- questionnaire_gen(n, cov_matrix = yw_cov, cat_prop = props, theta = TRUE,
                          family = "gaussian")
  colnames(bq)[2:4] <- c("thm", "thr", "ths")
  
  # checking that theta distributions are close to intended
  th_stats <- calc_m_v(bq, c("thm", "thr", "ths"))
  if (rep == 1) {
    cat("Theta statistics for replication 1:\n")
    print(th_stats)
  }
  
  # checking that correlations are close to intended
  overall_correlation <- cor(bq[,2:4]) 
  if (rep == 1) {
    cat("Correlation matrix for replication 1:\n")
    print(overall_correlation)
  }
  
  #####################################
  # 4. Simulating test for 3 domains
  #####################################
  
  # simulating items
  # math
  item_pool_m <- lsasim::item_gen(n_2pl = Im, b_bounds = c(-2.000, 2.000), 
                                  a_bounds = c(0.305, 2.097))
  
  # reading
  item_pool_r <- lsasim::item_gen(n_2pl = Ir, b_bounds = c(-2.000, 2.000), 
                                  a_bounds = c(0.377, 1.8277))
  
  # science
  item_pool_s <- lsasim::item_gen(n_2pl = Is, b_bounds = c(-2.000, 2.000), 
                                  a_bounds = c(0.478, 2.188))
  
  # simulating blocks
  # math
  K_m <- 18
  blocks_m <- lsasim::block_design(n_blocks = K_m, item_parameters = item_pool_m, 
                                   item_block_matrix = NULL)
  
  # reading
  K_r <- 10
  blocks_r <- lsasim::block_design(n_blocks = K_r, item_parameters = item_pool_r, 
                                   item_block_matrix = NULL)
  
  # science
  K_s <- 6
  blocks_s <- lsasim::block_design(n_blocks = K_s, item_parameters = item_pool_s, 
                                   item_block_matrix = NULL)
  
  # simulating booklets
  # math
  book_math <- matrix(c(
    1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 
  ), nrow=18, byrow=TRUE)
  books_m <- lsasim::booklet_design(item_block_assignment=blocks_m$block_assignment, 
                                    book_design = book_math)
  
  # reading
  book_reading <- matrix(c(
    1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=10, byrow=TRUE)
  
  books_r <- lsasim::booklet_design(item_block_assignment=blocks_r$block_assignment, 
                                    book_design = book_reading)
  
  # science
  book_science <- matrix(c(
    1, 1, 0, 0, 0, 0, 
    0, 1, 1, 0, 0, 0, 
    0, 0, 1, 1, 0, 0, 
    0, 0, 0, 1, 1, 0, 
    0, 0, 0, 0, 1, 1, 
    1, 0, 0, 0, 0, 1), nrow=6, byrow=TRUE)
  books_s <- lsasim::booklet_design(item_block_assignment=blocks_s$block_assignment,
                                    book_design = book_science)
  
  # assigning books to examinees
  bookassign_m <- lsasim::booklet_sample(n_subj=n, book_item_design=books_m, book_prob = NULL, 
                                         resample = FALSE, e = 0.1, iter = 20)
  
  bookassign_r <- lsasim::booklet_sample(n_subj=n, book_item_design=books_r, book_prob = NULL, 
                                         resample = FALSE, e = 0.1, iter = 20)
  
  bookassign_s <- lsasim::booklet_sample(n_subj=n, book_item_design=books_s, book_prob = NULL, 
                                         resample = FALSE, e = 0.1, iter = 20)
  
  #####################################
  # 5. Simulating responses to each test
  #####################################
  
  # math
  sim_math <- response_gen(subject = bookassign_m$subject, item=bookassign_m$item, 
                           theta = bq$thm,
                           a_par = item_pool_m$a, b_par = item_pool_m$b)
  
  # reading
  sim_read <- response_gen(subject = bookassign_r$subject, item=bookassign_r$item, 
                           theta = bq$thr,
                           a_par = item_pool_r$a, b_par = item_pool_r$b)
  
  # science
  sim_sci <- response_gen(subject = bookassign_s$subject, item=bookassign_s$item, 
                          theta = bq$ths,
                          a_par = item_pool_s$a, b_par = item_pool_s$b)
  
  # changing the name of the variables to "m"/"r"/"s" then item number
  names(sim_math) <- gsub("^i", "m", names(sim_math))
  names(sim_read) <- gsub("^i", "r", names(sim_read))
  names(sim_sci)  <- gsub("^i", "s", names(sim_sci))
  
  # converting q3 into pseudo-school variable with 154 levels
  breaks <- quantile(bq$q3, probs = seq(0, 1, length.out = 155), na.rm = TRUE)
  bq$school <- as.factor(cut(bq$q3, breaks = breaks, labels = FALSE, include.lowest = TRUE))
  
  # merging all of the test data together into a single dataframe
  simdat <- cbind(sim_math[, c(235, 1:234)], sim_read[,-ncol(sim_read)], sim_sci[,-ncol(sim_sci)])
  
  ###################################
  # 6. Processing background questionnaire data
  ###################################
  
  # converting ordinal variables to dummy coded numeric variables
  dums <- fastDummies::dummy_cols(bq[,9:14], remove_first_dummy = TRUE)
  
  # combining all bq data
  bqdum <- cbind(bq, dums[,7:172])
  
  # creating PCs
  pcbq <- princomp(bqdum[,c(6:8,15:27)]) # schools excluded
  
  # selecting first 11 PCs to attach to BQ data
  pcs <- pcbq$scores[,1:11]
  bqdum <- cbind(bqdum, pcs)
  
  ###################################
  # 7. Saving final data for this replication
  ###################################
  
  # creating an object with all item parameters
  item_pool_m$item <- sprintf("m%03d", item_pool_m$item)
  item_pool_r$item <- sprintf("r%03d", item_pool_r$item)
  item_pool_s$item <- sprintf("s%03d", item_pool_s$item)
  
  genparms <- rbind(item_pool_m, item_pool_r, item_pool_s)
  
  # Store data in lists for all replications
  if (rep == 1) {
    # Initialize lists on first replication
    all_simdat <- list()
    all_bqdum <- list()
    all_pcbq <- list()
    all_genparms <- list()
    all_summaries <- list()
  }I
  
  # Store current replication data
  all_simdat[[rep]] <- simdat
  all_bqdum[[rep]] <- bqdum
  all_pcbq[[rep]] <- pcbq
  all_genparms[[rep]] <- genparms
  
  # Save each replication separately (individual files)
  save(list = c("simdat", "bqdum", "pcbq", "genparms"), 
       file = paste0(path, "data/replications/alldata_rep", sprintf("%03d", rep), ".RData"))
  
  # Optional: Save a summary of key statistics for each replication
  rep_summary <- list(
    replication = rep,
    n_subjects = n,
    theta_stats = th_stats,
    theta_correlations = overall_correlation,
    item_counts = c(math = Im, reading = Ir, science = Is)
  )
  
  all_summaries[[rep]] <- rep_summary
  
  cat(sprintf("Completed replication %d of %d\n", rep, n_reps))
}

# Save summary of all replications AND all data in lists
save(all_summaries, file = paste0(path, "data/replications/replication_summaries.RData"))
save(list = c("all_simdat", "all_bqdum", "all_pcbq", "all_genparms", "all_summaries"), 
     file = paste0(path, "data/all_replications.RData"))

cat(sprintf("All %d replications completed successfully!\n", n_reps))
cat("Individual files saved in:", paste0(path, "data/replications/"), "\n")
cat("Individual files follow naming pattern: alldata_rep001.RData, alldata_rep002.RData, etc.\n")
cat("All replications combined saved as:", paste0(path, "data/all_replications.RData"), "\n")
cat("To access replication 5 simdat: load the file and use all_simdat[[5]]\n")