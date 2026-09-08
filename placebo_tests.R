# Funktioner til at udføre og visualisere placebo-tests for RDD-analyse

# 1. Placebo-tærskel test - tester effekter ved falske tærskelværdier
run_placebo_threshold_test <- function(data, 
                                       outcome_var = "miljø_afhængig",
                                       running_var = "centered_lagged_pervote_samlet", 
                                       true_cutpoint = 0,
                                       placebo_range = c(-5, 5),
                                       placebo_step = 0.5,
                                       polynomial = 1,
                                       bandwidth = NULL,
                                       covariates = NULL) {
  
  # Generér placebo-tærskelværdier, undtagen den sande tærskel
  placebo_thresholds <- seq(placebo_range[1], placebo_range[2], by = placebo_step)
  placebo_thresholds <- placebo_thresholds[abs(placebo_thresholds - true_cutpoint) > 0.01]
  
  # Håndtering af kovariater
  covs_matrix <- NULL
  if (!is.null(covariates)) {
    # Forsøg at bygge kovariat-matrix
    covs_list <- list()
    for (i in 1:length(covariates)) {
      col_name <- covariates[i]
      if (col_name %in% names(data)) {
        covs_list[[i]] <- data[[col_name]]
      } else {
        warning(paste("Kolonne", col_name, "findes ikke i datasættet."))
        covs_list[[i]] <- NULL
      }
    }
    # Fjern NULL værdier
    covs_list <- covs_list[!sapply(covs_list, is.null)]
    
    # Hvis der er nogen kovariater tilbage, lav en matrix
    if (length(covs_list) > 0) {
      covs_matrix <- do.call(cbind, covs_list)
    }
  }
  
  # Tom dataframe til resultater
  placebo_results <- data.frame()
  
  # Kør analysen for hver placebo-tærskel
  for (c in placebo_thresholds) {
    # Prøv at køre rdrobust med fejlhåndtering
    tryCatch({
      # Kør rdrobust 
      if (is.null(covs_matrix)) {
        rd <- rdrobust(
          y = data[[outcome_var]],
          x = data[[running_var]],
          c = c,  # Placebo tærskel
          p = polynomial,
          h = bandwidth
        )
      } else {
        rd <- rdrobust(
          y = data[[outcome_var]],
          x = data[[running_var]],
          c = c,  # Placebo tærskel
          p = polynomial,
          h = bandwidth,
          covs = covs_matrix
        )
      }
      
      # Gem estimate og standardfejl
      temp_results <- data.frame(
        threshold = c,
        effect = rd$coef[1],
        se = rd$se[1],
        t_stat = rd$coef[1] / rd$se[1],
        p_value = rd$pv[1],
        ci_lower = rd$coef[1] - 1.96 * rd$se[1],
        ci_upper = rd$coef[1] + 1.96 * rd$se[1]
      )
      
      placebo_results <- rbind(placebo_results, temp_results)
      
    }, error = function(e) {
      warning(paste("Fejl ved kørsel af rdrobust for tærskel =", c, ":", e$message))
    })
  }
  
  # Kør analysen med den sande tærskel
  tryCatch({
    if (is.null(covs_matrix)) {
      true_rd <- rdrobust(
        y = data[[outcome_var]],
        x = data[[running_var]],
        c = true_cutpoint,
        p = polynomial,
        h = bandwidth
      )
    } else {
      true_rd <- rdrobust(
        y = data[[outcome_var]],
        x = data[[running_var]],
        c = true_cutpoint,
        p = polynomial,
        h = bandwidth,
        covs = covs_matrix
      )
    }
    
    true_effect <- data.frame(
      threshold = true_cutpoint,
      effect = true_rd$coef[1],
      se = true_rd$se[1],
      t_stat = true_rd$coef[1] / true_rd$se[1],
      p_value = true_rd$pv[1],
      ci_lower = true_rd$coef[1] - 1.96 * true_rd$se[1],
      ci_upper = true_rd$coef[1] + 1.96 * true_rd$se[1]
    )
    
    # Tilføj den sande effekt til resultaterne
    placebo_results <- rbind(placebo_results, true_effect)
    
  }, error = function(e) {
    warning(paste("Fejl ved kørsel af rdrobust for sand tærskel:", e$message))
  })
  
  # Check om der er resultater
  if (nrow(placebo_results) == 0) {
    stop("Ingen resultater blev produceret. Tjek data og parametre.")
  }
  
  # Beregn empirisk p-værdi (andel af placebo-effekter med større absolut værdi end sand effekt)
  placebo_only <- placebo_results[placebo_results$threshold != true_cutpoint, ]
  true_only <- placebo_results[placebo_results$threshold == true_cutpoint, ]
  
  if (nrow(true_only) > 0) {
    emp_p_value <- mean(abs(placebo_only$effect) >= abs(true_only$effect))
  } else {
    emp_p_value <- NA
  }
  
  # Tilføj markering for den sande tærskel
  placebo_results$is_true <- placebo_results$threshold == true_cutpoint
  
  return(list(
    results = placebo_results,
    empirical_p_value = emp_p_value,
    true_effect = true_only$effect
  ))
}

# 2. Visualisering af placebo-tærskel resultater
plot_placebo_threshold_test <- function(placebo_results, 
                                        title = "Placebo-tærskel test") {
  
  # Hent resultater og empirical p-værdi
  results <- placebo_results$results
  emp_p_value <- placebo_results$empirical_p_value
  
  # Find den sande effekt
  true_effect <- results[results$is_true, ]
  
  # Lav plot af effektstørrelser
  p1 <- ggplot(results, aes(x = threshold, y = effect)) +
    geom_point(aes(color = is_true, size = is_true)) +
    geom_errorbar(aes(ymin = effect - 1.96*se, ymax = effect + 1.96*se, 
                      color = is_true, size = is_true), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
    scale_color_manual(values = c("gray50", "red")) +
    scale_size_manual(values = c(1, 2)) +
    labs(title = title,
         subtitle = paste("Empirisk p-værdi:", round(emp_p_value, 3)),
         x = "Placebo-tærskel",
         y = "Estimeret effekt") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 12),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "lightgray")
    )
  
  # Tilføj annotation for sand effekt
  if (nrow(true_effect) > 0) {
    true_effect_val <- true_effect$effect
    p1 <- p1 + 
      geom_text(
        data = data.frame(
          x = true_effect$threshold * 1.1,
          y = 0.05,
          label = paste("Sand effekt:", round(true_effect_val, 3))
        ),
        aes(x = x, y = y, label = label),
        color = "red",
        hjust = 0,
        inherit.aes = FALSE
      )
  }
  
  # Lav plot af t-statistikker
  p2 <- ggplot(results, aes(x = threshold, y = t_stat)) +
    geom_point(aes(color = is_true, size = is_true)) +
    geom_hline(yintercept = c(-1.96, 0, 1.96), 
               linetype = c("dotted", "dashed", "dotted"), 
               color = "darkgray") +
    scale_color_manual(values = c("gray50", "red")) +
    scale_size_manual(values = c(1, 2)) +
    labs(title = "Placebo-tærskel t-statistikker",
         x = "Placebo-tærskel",
         y = "t-statistik") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "lightgray")
    )
  
  # Lav density plot af placebo-effekter
  p3 <- ggplot(results, aes(x = effect)) +
    geom_density(fill = "gray80", alpha = 0.7) +
    geom_vline(data = true_effect, aes(xintercept = effect), 
               color = "red", size = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
    labs(title = "Fordeling af placebo-effekter",
         x = "Estimeret effekt",
         y = "Tæthed") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "lightgray")
    )
  
  # Tilføj annotation for sand effekt
  if (nrow(true_effect) > 0) {
    true_effect_val <- true_effect$effect
    p3 <- p3 + 
      geom_text(
        data = data.frame(
          x = true_effect_val * 1.1,
          y = max(density(results$effect)$y) * 0.5,
          label = paste("Sand effekt:", round(true_effect_val, 3))
        ),
        aes(x = x, y = y, label = label),
        color = "red",
        hjust = 0,
        inherit.aes = FALSE
      )
  }
  
  # Rankplot (absolut værdi)
  results_ranked <- results %>%
    mutate(abs_effect = abs(effect)) %>%
    arrange(abs_effect) %>%
    mutate(rank = row_number())
  
  true_rank <- results_ranked %>%
    filter(is_true) %>%
    pull(rank)
  
  p4 <- ggplot(results_ranked, aes(x = rank, y = abs_effect)) +
    geom_point(aes(color = is_true, size = is_true)) +
    geom_hline(data = results_ranked %>% filter(is_true), 
               aes(yintercept = abs_effect), 
               color = "red", linetype = "dashed") +
    labs(title = "Rangering af effektstørrelser (absolut værdi)",
         x = "Rang",
         y = "Absolut effektstørrelse") +
    scale_color_manual(values = c("gray50", "red")) +
    scale_size_manual(values = c(1, 2)) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "lightgray")
    )
  
  # Tilføj annotation for empirisk p-værdi
  if (!is.na(emp_p_value) && nrow(true_effect) > 0) {
    true_effect_abs <- abs(true_effect$effect)
    p4 <- p4 + 
      geom_text(
        data = data.frame(
          x = nrow(results_ranked) * 0.8,
          y = true_effect_abs * 1.1,
          label = paste("Empirisk p-værdi:", round(emp_p_value, 3))
        ),
        aes(x = x, y = y, label = label),
        color = "red",
        hjust = 0,
        inherit.aes = FALSE
      )
  }
  
  # Kombiner de fire plots, hvis gridExtra er tilgængelig
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    combined_plot <- gridExtra::grid.arrange(
      p1, p2, p3, p4, 
      ncol = 2, 
      top = grid::textGrob(title, gp = grid::gpar(fontsize = 16, fontface = "bold"))
    )
    return(combined_plot)
  } else {
    # Hvis gridExtra ikke er tilgængelig, returner kun det første plot
    return(p1)
  }
}

# 3. Outcome placebo test - tester effekter for variable der ikke burde påvirkes
run_outcome_placebo_test <- function(data, 
                                     real_outcome_var = "miljø_afhængig",
                                     placebo_outcome_vars,
                                     running_var = "centered_lagged_pervote_samlet", 
                                     cutpoint = 0,
                                     polynomial = 1,
                                     bandwidth = NULL,
                                     covariates = NULL) {
  
  # Kombinér den rigtige outcome variabel med placebo variables
  all_outcomes <- c(real_outcome_var, placebo_outcome_vars)
  
  # Håndtering af kovariater
  covs_matrix <- NULL
  if (!is.null(covariates)) {
    # Forsøg at bygge kovariat-matrix
    covs_list <- list()
    for (i in 1:length(covariates)) {
      col_name <- covariates[i]
      if (col_name %in% names(data)) {
        covs_list[[i]] <- data[[col_name]]
      } else {
        warning(paste("Kolonne", col_name, "findes ikke i datasættet."))
        covs_list[[i]] <- NULL
      }
    }
    # Fjern NULL værdier
    covs_list <- covs_list[!sapply(covs_list, is.null)]
    
    # Hvis der er nogen kovariater tilbage, lav en matrix
    if (length(covs_list) > 0) {
      covs_matrix <- do.call(cbind, covs_list)
    }
  }
  
  # Tom dataframe til resultater
  placebo_results <- data.frame()
  
  # Kør analysen for hver outcome variabel
  for (y_var in all_outcomes) {
    # Prøv at køre rdrobust med fejlhåndtering
    tryCatch({
      # Kør rdrobust 
      if (is.null(covs_matrix)) {
        rd <- rdrobust(
          y = data[[y_var]],
          x = data[[running_var]],
          c = cutpoint,
          p = polynomial,
          h = bandwidth
        )
      } else {
        rd <- rdrobust(
          y = data[[y_var]],
          x = data[[running_var]],
          c = cutpoint,
          p = polynomial,
          h = bandwidth,
          covs = covs_matrix
        )
      }
      
      # Gem estimate og standardfejl
      temp_results <- data.frame(
        outcome_var = y_var,
        effect = rd$coef[1],
        se = rd$se[1],
        t_stat = rd$coef[1] / rd$se[1],
        p_value = rd$pv[1],
        ci_lower = rd$coef[1] - 1.96 * rd$se[1],
        ci_upper = rd$coef[1] + 1.96 * rd$se[1]
      )
      
      placebo_results <- rbind(placebo_results, temp_results)
      
    }, error = function(e) {
      warning(paste("Fejl ved kørsel af rdrobust for variabel =", y_var, ":", e$message))
    })
  }
  
  # Check om der er resultater
  if (nrow(placebo_results) == 0) {
    stop("Ingen resultater blev produceret. Tjek data og parametre.")
  }
  
  # Markér hvilken variabel der er den rigtige outcome
  placebo_results$is_real <- placebo_results$outcome_var == real_outcome_var
  
  return(placebo_results)
}

# 4. Visualisering af outcome placebo resultater
plot_outcome_placebo_test <- function(placebo_results, 
                                      title = "Outcome Placebo Test") {
  
  # Sørg for at outcome_var er en faktor med passende rækkefølge
  placebo_results$outcome_var <- factor(placebo_results$outcome_var, 
                                        levels = unique(placebo_results$outcome_var))
  
  # Lav coefficient plot
  p <- ggplot(placebo_results, aes(x = outcome_var, y = effect, 
                                   color = is_real, fill = is_real)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                  width = 0.2, position = position_dodge(0.9)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
    scale_color_manual(values = c("gray50", "red"), 
                       labels = c("Placebo", "Rigtig")) +
    scale_fill_manual(values = c("gray80", "red"), 
                      labels = c("Placebo", "Rigtig")) +
    labs(title = title,
         x = "Outcome variabel",
         y = "Estimeret effekt") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "lightgray"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}

# 5. Permutations test - tilfældig reassignering af treatment-status
run_permutation_test <- function(data, 
                                 outcome_var = "miljø_afhængig",
                                 running_var = "centered_lagged_pervote_samlet", 
                                 treatment_var = "lagged_i_parlament",
                                 cutpoint = 0,
                                 polynomial = 1,
                                 bandwidth = NULL,
                                 n_permutations = 1000,
                                 covariates = NULL) {
  
  # Håndtering af kovariater
  covs_matrix <- NULL
  if (!is.null(covariates)) {
    # Forsøg at bygge kovariat-matrix
    covs_list <- list()
    for (i in 1:length(covariates)) {
      col_name <- covariates[i]
      if (col_name %in% names(data)) {
        covs_list[[i]] <- data[[col_name]]
      } else {
        warning(paste("Kolonne", col_name, "findes ikke i datasættet."))
        covs_list[[i]] <- NULL
      }
    }
    # Fjern NULL værdier
    covs_list <- covs_list[!sapply(covs_list, is.null)]
    
    # Hvis der er nogen kovariater tilbage, lav en matrix
    if (length(covs_list) > 0) {
      covs_matrix <- do.call(cbind, covs_list)
    }
  }
  
  # Estimér den faktiske effekt
  tryCatch({
    # Kør rdrobust med de faktiske data
    if (is.null(covs_matrix)) {
      actual_rd <- rdrobust(
        y = data[[outcome_var]],
        x = data[[running_var]],
        c = cutpoint,
        p = polynomial,
        h = bandwidth
      )
    } else {
      actual_rd <- rdrobust(
        y = data[[outcome_var]],
        x = data[[running_var]],
        c = cutpoint,
        p = polynomial,
        h = bandwidth,
        covs = covs_matrix
      )
    }
    
    actual_effect <- actual_rd$coef[1]
    actual_t_stat <- actual_rd$coef[1] / actual_rd$se[1]
    
  }, error = function(e) {
    stop(paste("Fejl ved estimation af faktisk effekt:", e$message))
  })
  
  # Kør permutationer
  permutation_effects <- numeric(n_permutations)
  permutation_t_stats <- numeric(n_permutations)
  
  for (i in 1:n_permutations) {
    # Lav en kopi af data med permuteret treatment
    perm_data <- data
    
    # Permutér treatment-variablen
    perm_data[[treatment_var]] <- sample(perm_data[[treatment_var]])
    
    # Omdefiner running-variablen ift. den permuterede treatment
    # (Dette afhænger af din specifikke RDD-setup - justér efter behov)
    
    # Alternativt, for en sharp RDD tilgang, kan man permutere around the cutpoint
    perm_data$perm_above <- as.numeric(perm_data[[running_var]] >= cutpoint)
    
    # Estimér effekten med permuterede data
    tryCatch({
      perm_model <- ivreg(
        formula = as.formula(paste(
          outcome_var, "~", treatment_var, "|", running_var, ">= 0"
        )),
        data = perm_data
      )
      
      permutation_effects[i] <- coef(perm_model)[treatment_var]
      permutation_t_stats[i] <- coef(summary(perm_model))[treatment_var, "t value"]
      
    }, error = function(e) {
      warning(paste("Fejl i permutation", i, ":", e$message))
      # Sæt missing værdier hvis fejl opstår
      permutation_effects[i] <- NA
      permutation_t_stats[i] <- NA
    })
  }
  
  # Fjern missing værdier
  permutation_effects <- permutation_effects[!is.na(permutation_effects)]
  permutation_t_stats <- permutation_t_stats[!is.na(permutation_t_stats)]
  
  # Beregn empirisk p-værdi
  p_value_effect <- mean(abs(permutation_effects) >= abs(actual_effect))
  p_value_t_stat <- mean(abs(permutation_t_stats) >= abs(actual_t_stat))
  
  # Returner resultater
  return(list(
    actual_effect = actual_effect,
    actual_t_stat = actual_t_stat,
    permutation_effects = permutation_effects,
    permutation_t_stats = permutation_t_stats,
    p_value_effect = p_value_effect,
    p_value_t_stat = p_value_t_stat
  ))
}

# 6. Visualisering af permutationstest
plot_permutation_test <- function(permutation_results, 
                                  title = "Permutationstest") {
  
  # Lav dataframe af permuterede effekter
  perm_df <- data.frame(
    effect = permutation_results$permutation_effects
  )
  
  # Lav density plot
  p1 <- ggplot(perm_df, aes(x = effect)) +
    geom_density(fill = "lightblue", alpha = 0.7) +
    geom_vline(xintercept = permutation_results$actual_effect, 
               color = "red", size = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
    labs(title = title,
         subtitle = paste("Baseret på", length(permutation_results$permutation_effects), "permutationer"),
         x = "Estimeret effekt",
         y = "Tæthed") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 12),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "lightgray")
    )
  
  # Tilføj annotation for faktisk effekt og p-værdi
  p1 <- p1 + 
    geom_text(
      data = data.frame(
        x = permutation_results$actual_effect * 1.2,
        y = max(density(perm_df$effect)$y) * 0.5,
        label = paste("Faktisk effekt:", 
                      round(permutation_results$actual_effect, 3),
                      "\np-værdi:", 
                      round(permutation_results$p_value_effect, 3))
      ),
      aes(x = x, y = y, label = label),
      color = "red", 
      hjust = 0,
      inherit.aes = FALSE
    )
  
  # Lav t-stat density plot
  perm_t_df <- data.frame(
    t_stat = permutation_results$permutation_t_stats
  )
  
  p2 <- ggplot(perm_t_df, aes(x = t_stat)) +
    geom_density(fill = "lightgreen", alpha = 0.7) +
    geom_vline(xintercept = permutation_results$actual_t_stat, 
               color = "red", size = 1) +
    geom_vline(xintercept = c(-1.96, 0, 1.96), 
               linetype = c("dotted", "dashed", "dotted"), 
               color = "darkgray") +
    labs(title = "Permutationstest (t-statistikker)",
         x = "t-statistik",
         y = "Tæthed") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "lightgray")
    )
  
  # Tilføj annotation for faktisk t-statistik og p-værdi
  p2 <- p2 + 
    geom_text(
      data = data.frame(
        x = permutation_results$actual_t_stat * 1.2,
        y = max(density(perm_t_df$t_stat)$y) * 0.5,
        label = paste("Faktisk t-statistik:", 
                      round(permutation_results$actual_t_stat, 3),
                      "\np-værdi:", 
                      round(permutation_results$p_value_t_stat, 3))
      ),
      aes(x = x, y = y, label = label),
      color = "red", 
      hjust = 0,
      inherit.aes = FALSE
    )
  
  # Kombiner plots, hvis gridExtra er tilgængelig
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    combined_plot <- gridExtra::grid.arrange(
      p1, p2, 
      ncol = 1
    )
    return(combined_plot)
  } else {
    return(p1)
  }
}