fit_wemix_by_group_pv <- function(data, 
                                  grouping_var, 
                                  model_formula, 
                                  weight_vars, 
                                  pv_vars = NULL,
                                  pv_root_name = NULL,
                                  pv_value_name = NULL,
                                  benchmark_vars = NULL,
                                  benchmark_root_name = NULL,
                                  benchmark_value_name = NULL,
                                  max_iter = 10, 
                                  nQuad = 13L, 
                                  cores = parallel::detectCores() - 10,
                                  pool_results = TRUE,
                                  debug = FALSE) {
  
  # IMMEDIATE DEBUG OUTPUT
  cat("=== FUNCTION START ===\n")
  cat("Debug mode:", debug, "\n")
  cat("Data dimensions:", nrow(data), "x", ncol(data), "\n")
  cat("Grouping variable:", grouping_var, "\n")
  cat("Grouping variable exists:", grouping_var %in% names(data), "\n")
  if (grouping_var %in% names(data)) {
    cat("Unique groups:", length(unique(data[[grouping_var]])), "\n")
  }
  
  library(WeMix)
  library(pbapply)
  library(parallel)
  library(tidyr)
  
  cat("Libraries loaded successfully\n")
  
  # Input validation
  cat("Starting input validation...\n")
  
  if (!is.null(pv_vars) && !is.null(pv_root_name)) {
    stop("Specify either pv_vars OR pv_root_name, not both")
  }
  
  if (!is.null(benchmark_vars) && !is.null(benchmark_root_name)) {
    stop("Specify either benchmark_vars OR benchmark_root_name, not both")
  }
  
  cat("Input validation passed\n")
  
  # Helper function to identify PV columns
  identify_pv_columns <- function(data, pv_vars = NULL, pv_root_name = NULL) {
    if (!is.null(pv_vars)) {
      return(pv_vars)
    } else if (!is.null(pv_root_name)) {
      # Find columns that start with the root name
      pv_cols <- grep(paste0("^", pv_root_name), names(data), value = TRUE)
      if (length(pv_cols) == 0) {
        stop(paste("No columns found starting with", pv_root_name))
      }
      return(pv_cols)
    }
    return(NULL)
  }
  
  # Identify PV and benchmark columns
  cat("Identifying PV and benchmark columns...\n")
  cat("pv_vars:", paste(pv_vars, collapse = ", "), "\n")
  cat("pv_root_name:", pv_root_name, "\n")
  cat("benchmark_vars:", paste(benchmark_vars, collapse = ", "), "\n")
  cat("benchmark_root_name:", benchmark_root_name, "\n")
  
  pv_columns <- identify_pv_columns(data, pv_vars, pv_root_name)
  benchmark_columns <- identify_pv_columns(data, benchmark_vars, benchmark_root_name)
  
  cat("PV columns found:", paste(pv_columns, collapse = ", "), "\n")
  cat("Benchmark columns found:", paste(benchmark_columns, collapse = ", "), "\n")
  
  # Debug: Print initial country counts
  cat("=== INITIAL DATA CHECK ===\n")
  cat("Initial countries in data:", length(unique(data[[grouping_var]])), "\n")
  if (debug) {
    cat("Countries:", sort(unique(data[[grouping_var]])), "\n")
  }
  
  # Prepare data based on PV and benchmark specifications
  cat("=== DATA PREPARATION ===\n")
  cat("Has PV columns:", !is.null(pv_columns), "\n")
  cat("Has benchmark columns:", !is.null(benchmark_columns), "\n")
  if (!is.null(pv_columns) || !is.null(benchmark_columns)) {
    
    # Create tall format for PVs
    if (!is.null(pv_columns)) {
      # Determine the root name for the new column
      if (!is.null(pv_value_name)) {
        final_pv_name <- pv_value_name
      } else if (!is.null(pv_root_name)) {
        final_pv_name <- pv_root_name
      } else {
        # Extract common root from pv_vars
        final_pv_name <- gsub("\\d+.*$", "", pv_columns[1])
      }
      
      data_tall <- pivot_longer(data, 
                                cols = all_of(pv_columns),
                                names_to = 'pv_name',
                                values_to = final_pv_name)
    } else {
      data_tall <- data
      data_tall$pv_name <- 'single_pv'
    }
    
    # Handle benchmarks similarly
    if (!is.null(benchmark_columns)) {
      # Determine the root name for the new column
      if (!is.null(benchmark_value_name)) {
        final_benchmark_name <- benchmark_value_name
      } else if (!is.null(benchmark_root_name)) {
        final_benchmark_name <- benchmark_root_name
      } else {
        # Extract common root from benchmark_vars
        final_benchmark_name <- gsub("\\d+.*$", "", benchmark_columns[1])
      }
      
      data_tall <- pivot_longer(data_tall, 
                                cols = all_of(benchmark_columns),
                                names_to = 'benchmark_name',
                                values_to = final_benchmark_name)
    } else {
      data_tall$benchmark_name <- 'single_benchmark'
    }
    
    # Create combinations for analysis
    if (!is.null(pv_columns) && !is.null(benchmark_columns)) {
      # Both PVs and benchmarks - create paired combinations (pv1 with benchmark1, etc.)
      # Check that we have equal numbers
      if (length(pv_columns) != length(benchmark_columns)) {
        stop("Number of PV columns must equal number of benchmark columns for paired analysis")
      }
      
      # Create paired combinations by extracting the number/index from column names
      pv_numbers <- gsub(".*?(\\d+).*", "\\1", pv_columns)
      benchmark_numbers <- gsub(".*?(\\d+).*", "\\1", benchmark_columns)
      
      # Create pairing identifier
      data_tall$pv_number <- gsub(".*?(\\d+).*", "\\1", data_tall$pv_name)
      data_tall$benchmark_number <- gsub(".*?(\\d+).*", "\\1", data_tall$benchmark_name)
      
      # Keep only matching pairs
      data_tall <- data_tall[data_tall$pv_number == data_tall$benchmark_number, ]
      
      # Debug: Check data after pairing
      if (debug) {
        cat("Countries after pairing:", length(unique(data_tall[[grouping_var]])), "\n")
        missing_countries <- setdiff(unique(data[[grouping_var]]), unique(data_tall[[grouping_var]]))
        if (length(missing_countries) > 0) {
          cat("Missing countries after pairing:", missing_countries, "\n")
        }
      }
      
      # Create combination identifier
      data_tall$pv_benchmark_combo <- paste("pv", data_tall$pv_number, sep = "")
      dat_list <- split(data_tall, data_tall$pv_benchmark_combo)
    } else if (!is.null(pv_columns)) {
      # Only PVs
      dat_list <- split(data_tall, data_tall$pv_name)
    } else {
      # Only benchmarks
      dat_list <- split(data_tall, data_tall$benchmark_name)
    }
    
  } else {
    # No PVs or benchmarks - use original data
    dat_list <- list(original = data)
  }
  
  # Fixed version of the run_wemix_analysis function
  run_wemix_analysis <- function(dataset) {
    # Identify unique groups
    groups <- unique(dataset[[grouping_var]])
    
    # Debug: Print group information
    if (debug) {
      cat("Processing", length(groups), "groups:", groups, "\n")
    }
    
    # Function to run single group analysis
    analyze_single_group <- function(g) {
      subset_data <- subset(dataset, dataset[[grouping_var]] == g)
      out <- list(group = g, fixed = NULL, varcomp = NULL, fit = NULL, convergence = FALSE, error = NA, n_obs = nrow(subset_data))
      
      # Check if subset has sufficient data
      if (nrow(subset_data) == 0) {
        out$error <- "No data for this group"
        out$convergence <- FALSE
        return(out)
      }
      
      # Check for minimum observations
      if (nrow(subset_data) < 10) {
        out$error <- paste("Insufficient observations:", nrow(subset_data))
        out$convergence <- FALSE
        return(out)
      }
      
      tryCatch({
        model <- WeMix::mix(
          formula = model_formula,
          data = subset_data,
          weights = weight_vars,
          max_iteration = max_iter,
          nQuad = nQuad,
          run = TRUE,
          verbose = FALSE
        )
        
        # Check if model actually converged
        if (is.null(model) || is.null(model$coef)) {
          out$error <- "Model returned NULL or no coefficients"
          out$convergence <- FALSE
          return(out)
        }
        
        # Extract outputs
        out$fixed <- data.frame(term = names(model$coef), estimate = model$coef, SE = model$SE)
        out$varcomp <- model$varDF[, c("grp", "vcov")]
        out$fit <- data.frame(logLik = model$lnl, ICC = model$ICC, sigma = model$sigma)
        out$convergence <- TRUE
        out$error <- NA
        
      }, error = function(e) {
        # Capture the full error message
        out$error <- paste("Error:", conditionMessage(e))
        out$convergence <- FALSE
      })
      
      return(out)
    }
    
    # Try parallel processing with better error handling
    tryCatch({
      # Create cluster with reduced cores and better setup
      cl <- makeCluster(min(cores, length(groups), parallel::detectCores() - 1))
      
      # Export only necessary variables and ensure clean environment
      clusterExport(cl, varlist = c("model_formula", "weight_vars", "max_iter", "nQuad", "grouping_var"), envir = environment())
      clusterEvalQ(cl, {
        library(WeMix)
      })
      
      # Run analysis in parallel
      results <- pblapply(groups, cl = cl, function(g) {
        # Re-subset data within worker to avoid serialization issues
        subset_data <- subset(dataset, dataset[[grouping_var]] == g)
        out <- list(group = g, fixed = NULL, varcomp = NULL, fit = NULL, convergence = FALSE, error = NA, n_obs = nrow(subset_data))
        
        if (nrow(subset_data) == 0) {
          out$error <- "No data for this group"
          return(out)
        }
        
        if (nrow(subset_data) < 10) {
          out$error <- paste("Insufficient observations:", nrow(subset_data))
          return(out)
        }
        
        tryCatch({
          model <- WeMix::mix(
            formula = model_formula,
            data = subset_data,
            weights = weight_vars,
            max_iteration = max_iter,
            nQuad = nQuad,
            run = TRUE,
            verbose = FALSE
          )
          
          if (is.null(model) || is.null(model$coef)) {
            out$error <- "Model returned NULL or no coefficients"
            return(out)
          }
          
          out$fixed <- data.frame(term = names(model$coef), estimate = model$coef, SE = model$SE)
          out$varcomp <- model$varDF[, c("grp", "vcov")]
          out$fit <- data.frame(logLik = model$lnl, ICC = model$ICC, sigma = model$sigma)
          out$convergence <- TRUE
          out$error <- NA
          
        }, error = function(e) {
          out$error <- paste("WeMix Error:", conditionMessage(e))
          out$convergence <- FALSE
        })
        
        return(out)
      })
      
      stopCluster(cl)
      
    }, error = function(e) {
      # If parallel processing fails, fall back to sequential processing
      if (exists("cl")) {
        tryCatch(stopCluster(cl), error = function(e) {})
      }
      
      cat("Parallel processing failed, switching to sequential processing...\n")
      cat("Error was:", conditionMessage(e), "\n")
      
      # Sequential processing as fallback
      results <- lapply(groups, analyze_single_group)
    })
    
    # Compile output lists - preserve all error messages
    convergence_df <- do.call(rbind, lapply(results, function(x) {
      data.frame(
        group = x$group,
        convergence = x$convergence,
        n_obs = x$n_obs,
        error = if (!is.null(x$error) && !is.na(x$error)) x$error else NA,
        stringsAsFactors = FALSE
      )
    }))
    
    list(
      fixed_effects = do.call(rbind, lapply(results, function(x) if (!is.null(x$fixed)) cbind(group = x$group, x$fixed))),
      variance_components = do.call(rbind, lapply(results, function(x) if (!is.null(x$varcomp)) cbind(group = x$group, x$varcomp))),
      fit_stats = do.call(rbind, lapply(results, function(x) if (!is.null(x$fit)) cbind(group = x$group, x$fit))),
      convergence_log = convergence_df
    )
  }
  
  # Run analysis on each dataset
  cat("=== RUNNING ANALYSIS ===\n")
  cat("Number of datasets to analyze:", length(dat_list), "\n")
  cat("Dataset names:", names(dat_list), "\n")
  
  all_results <- lapply(dat_list, run_wemix_analysis)
  
  cat("Analysis completed\n")
  cat("Results structure:", str(all_results), "\n")
  
  # FIXED: Better pooling logic
  if (pool_results && length(all_results) > 1) {
    
    # Create a simple pooled summary instead of relying on WeMix's summary method
    pooled_results <- list()
    
    # Pool fixed effects across all PV/benchmark combinations
    all_fixed <- do.call(rbind, lapply(names(all_results), function(combo_name) {
      if (!is.null(all_results[[combo_name]]$fixed_effects)) {
        cbind(pv_combo = combo_name, all_results[[combo_name]]$fixed_effects)
      }
    }))
    
    # Pool variance components
    all_varcomp <- do.call(rbind, lapply(names(all_results), function(combo_name) {
      if (!is.null(all_results[[combo_name]]$variance_components)) {
        cbind(pv_combo = combo_name, all_results[[combo_name]]$variance_components)
      }
    }))
    
    # Pool fit statistics
    all_fit <- do.call(rbind, lapply(names(all_results), function(combo_name) {
      if (!is.null(all_results[[combo_name]]$fit_stats)) {
        cbind(pv_combo = combo_name, all_results[[combo_name]]$fit_stats)
      }
    }))
    
    # Combine convergence logs
    all_convergence <- do.call(rbind, lapply(names(all_results), function(combo_name) {
      if (!is.null(all_results[[combo_name]]$convergence_log)) {
        cbind(pv_combo = combo_name, all_results[[combo_name]]$convergence_log)
      }
    }))
    
    # Create pooled summary statistics using Rubin's Rules - BASE R APPROACH
    if (!is.null(all_fixed) && nrow(all_fixed) > 0) {
      tryCatch({
        cat("=== FIXED EFFECTS POOLING DEBUG ===\n")
        cat("Total rows in all_fixed:", nrow(all_fixed), "\n")
        
        # Ensure numeric types
        all_fixed$estimate <- as.numeric(as.character(all_fixed$estimate))
        all_fixed$SE <- as.numeric(as.character(all_fixed$SE))
        
        # Get unique combinations
        unique_combos <- unique(all_fixed[, c("group", "term")])
        cat("Number of unique group-term combinations:", nrow(unique_combos), "\n")
        
        # Initialize results list
        pooled_list <- list()
        
        # Process each group-term combination separately using base R
        for (i in 1:nrow(unique_combos)) {
          grp <- unique_combos$group[i]
          trm <- unique_combos$term[i]
          
          # Subset data for this combination
          subset_data <- all_fixed[all_fixed$group == grp & all_fixed$term == trm, ]
          
          # Remove any NA values
          subset_data <- subset_data[!is.na(subset_data$estimate) & !is.na(subset_data$SE), ]
          
          if (nrow(subset_data) > 0) {
            # Calculate components
            m <- nrow(subset_data)
            estimates <- subset_data$estimate
            ses <- subset_data$SE
            
            # Point estimate
            pooled_estimate <- mean(estimates)
            
            # Within-imputation variance
            within_var <- mean(ses^2)
            
            # Between-imputation variance - manual calculation
            if (m > 1) {
              between_var <- var(estimates)  # This should work in base R
            } else {
              between_var <- 0
            }
            
            # Debug output for problematic case
            if (grp == "7842" && trm == "scale(mathach)") {
              cat("=== DEBUG FOR 7842 scale(mathach) ===\n")
              cat("Number of estimates:", m, "\n")
              cat("Estimates:", paste(estimates, collapse = ", "), "\n")
              cat("Mean estimate:", pooled_estimate, "\n")
              cat("Variance of estimates:", between_var, "\n")
              cat("SEs:", paste(ses, collapse = ", "), "\n")
              cat("Within variance:", within_var, "\n")
            }
            
            # Total variance using Rubin's rules
            total_var <- within_var + between_var + (between_var / m)
            
            # Pooled SE
            pooled_SE <- sqrt(total_var)
            
            # Degrees of freedom
            if (m <= 1 || between_var == 0) {
              df <- Inf
            } else {
              df <- (m - 1) * (1 + within_var / ((1 + 1/m) * between_var))^2
            }
            
            # Store results
            pooled_list[[i]] <- data.frame(
              group = grp,
              term = trm,
              estimate = pooled_estimate,
              SE = sqrt(within_var),
              within_var = within_var,
              between_var = between_var,
              total_var = total_var,
              pooled_SE = pooled_SE,
              df = df,
              m = m,
              stringsAsFactors = FALSE
            )
          }
        }
        
        # Combine results
        pooled_fixed <- do.call(rbind, pooled_list)
        
        cat("Sample of pooled results (first 5 rows):\n")
        print(head(pooled_fixed[, c("group", "term", "estimate", "SE", "within_var", "between_var", "total_var", "pooled_SE", "m")], 5))
        
        # Check for problematic cases
        na_pooled <- sum(is.na(pooled_fixed$pooled_SE))
        zero_between <- sum(pooled_fixed$between_var == 0, na.rm = TRUE)
        nonzero_between <- sum(pooled_fixed$between_var > 0, na.rm = TRUE)
        
        cat("NA pooled SEs:", na_pooled, "\n")
        cat("Zero between variance:", zero_between, "\n") 
        cat("Non-zero between variance:", nonzero_between, "\n")
        
        # Show distribution of between variance
        cat("Summary of between variance:\n")
        print(summary(pooled_fixed$between_var))
        
        pooled_results$pooled_fixed_effects <- pooled_fixed
        
        cat("Fixed effects pooled using Rubin's rules (base R)\n")
        cat("Number of parameter estimates:", nrow(pooled_fixed), "\n")
        
      }, error = function(e) {
        cat("ERROR in pooling fixed effects:", conditionMessage(e), "\n")
        traceback()
        pooled_results$pooled_fixed_effects <- NULL
      })
    }
    
    if (!is.null(all_varcomp) && nrow(all_varcomp) > 0) {
      tryCatch({
        # Apply Rubin's Rules to variance components
        library(dplyr)
        pooled_varcomp <- all_varcomp %>%
          group_by(group, grp) %>%
          summarise(
            # Point estimate: average across imputations
            vcov = mean(vcov, na.rm = TRUE),
            
            # Between-imputation variance
            between_var = var(vcov, na.rm = TRUE),
            
            # Number of imputations  
            m = n(),
            
            .groups = 'drop'
          ) %>%
          mutate(
            # Handle case where between_var is NA
            between_var = ifelse(is.na(between_var) | m == 1, 0, between_var),
            
            # For variance components, we don't have within-imputation SEs from WeMix
            # So we use a simplified approach: SE = sqrt(between_var + between_var/m)
            pooled_SE = sqrt(between_var + (between_var / m)),
            
            # Set to NA if no between-imputation variance
            pooled_SE = ifelse(between_var == 0, NA, pooled_SE)
          ) %>%
          select(group, grp, vcov, pooled_SE, between_var, m) %>%
          as.data.frame()
        
        pooled_results$pooled_variance_components <- pooled_varcomp
        
        cat("Variance components pooled using modified Rubin's rules\n")
        cat("Number of groups with variance components:", nrow(pooled_varcomp), "\n")
        
      }, error = function(e) {
        cat("Warning: Could not pool variance components:", conditionMessage(e), "\n")
        pooled_results$pooled_variance_components <- NULL
      })
    }
    
    if (!is.null(all_fit) && nrow(all_fit) > 0) {
      tryCatch({
        # Use a safer approach for pooling fit statistics
        library(dplyr)
        pooled_fit <- all_fit %>%
          group_by(group) %>%
          summarise(
            logLik = mean(logLik, na.rm = TRUE),
            ICC = mean(ICC, na.rm = TRUE),
            sigma = mean(sigma, na.rm = TRUE),
            n_imputations = n(),
            .groups = 'drop'
          ) %>%
          as.data.frame()
        
        pooled_results$pooled_fit_stats <- pooled_fit
      }, error = function(e) {
        cat("Warning: Could not pool fit statistics:", conditionMessage(e), "\n")
        pooled_results$pooled_fit_stats <- NULL
      })
    }
    
    # Summary of convergence across all combinations
    if (!is.null(all_convergence) && nrow(all_convergence) > 0) {
      tryCatch({
        library(dplyr)
        convergence_summary <- all_convergence %>%
          group_by(group) %>%
          summarise(
            n_converged = sum(convergence, na.rm = TRUE),
            n_total = n(),
            prop_converged = mean(convergence, na.rm = TRUE),
            .groups = 'drop'
          ) %>%
          as.data.frame()
        
        pooled_results$convergence_summary <- convergence_summary
      }, error = function(e) {
        cat("Warning: Could not create convergence summary:", conditionMessage(e), "\n")
        pooled_results$convergence_summary <- NULL
      })
    }
    
    # Return both individual and pooled results
    return(list(
      individual_results = all_results,
      pooled_results = pooled_results,
      raw_combined = list(
        all_fixed = all_fixed,
        all_varcomp = all_varcomp,
        all_fit = all_fit,
        all_convergence = all_convergence
      )
    ))
    
  } else {
    # Return individual results without pooling
    if (length(all_results) == 1) {
      return(all_results[[1]])
    } else {
      return(all_results)
    }
  }
}

# Example usage function
example_usage <- function() {
  # Example 1: Using PV vars with custom column names
  # results <- fit_wemix_by_group_pv(
  #   data = pisa_data,
  #   grouping_var = "country",
  #   model_formula = mathinterest ~ scale(achievement) + benchmark_bin + (1|schoolid),
  #   weight_vars = c('w_fstuwt', 'w_fschwt'),
  #   pv_vars = c("pvm1", "pvm2", "pvm3", "pvm4", "pvm5"),
  #   pv_value_name = "achievement",
  #   benchmark_vars = c("bmm1_bin", "bmm2_bin", "bmm3_bin", "bmm4_bin", "bmm5_bin"),
  #   benchmark_value_name = "benchmark_bin",
  #   pool_results = TRUE,
  #   debug = TRUE
  # )
  
  # Example 2: Using root names (will auto-detect column names)
  # results <- fit_wemix_by_group_pv(
  #   data = pisa_data,
  #   grouping_var = "country",
  #   model_formula = mathinterest ~ scale(pv) + bmm + (1|schoolid),
  #   weight_vars = c('w_fstuwt', 'w_fschwt'),
  #   pv_root_name = "pv",
  #   benchmark_root_name = "bmm",
  #   pool_results = TRUE,
  #   debug = TRUE
  # )
}