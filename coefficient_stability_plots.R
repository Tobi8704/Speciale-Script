# Funktion til at generere coefficient stability plots
create_coefficient_stability_plot <- function(data, 
                                              outcome_var = "miljø_afhængig",
                                              running_var = "centered_lagged_pervote_samlet", 
                                              bw_list = seq(1, 3, by = 0.5),
                                              polynomials = c(1, 2, 3),
                                              cutpoint = 0,
                                              covariates = NULL,
                                              title = "Coefficient Stability Plot") {
  
  results <- data.frame()
  covs_matrix <- NULL
  
  if (!is.null(covariates)) {
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
    covs_list <- covs_list[!sapply(covs_list, is.null)]
    if (length(covs_list) > 0) {
      covs_matrix <- do.call(cbind, covs_list)
    }
  }
  
  for (p in polynomials) {
    for (bw in bw_list) {
      tryCatch({
        if (is.null(covs_matrix)) {
          rd <- rdrobust(
            y = data[[outcome_var]],
            x = data[[running_var]],
            c = cutpoint,
            p = p,
            h = bw
          )
        } else {
          rd <- rdrobust(
            y = data[[outcome_var]],
            x = data[[running_var]],
            c = cutpoint,
            p = p,
            h = bw,
            covs = covs_matrix
          )
        }
        
        results <- rbind(results, data.frame(
          bandwidth = bw,
          polynomial = p,
          coefficient = rd$coef[1],
          se = rd$se[1],
          ci_lower = rd$coef[1] - 1.96 * rd$se[1],
          ci_upper = rd$coef[1] + 1.96 * rd$se[1],
          p_value = rd$pv[1]
        ))
        
      }, error = function(e) {
        warning(paste("Fejl ved kørsel af rdrobust med p =", p, "og bw =", bw, ":", e$message))
      })
    }
  }
  
  if (nrow(results) == 0) {
    stop("Ingen resultater blev produceret. Tjek data og parametre.")
  }
  
  results$polynomial <- factor(results$polynomial, 
                               levels = polynomials,
                               labels = paste0("Polynomial ", polynomials))
  
  poly_colors <- c("blue", "red", "green", "purple", "orange")[1:length(polynomials)]
  
  p <- ggplot(results, aes(x = bandwidth, y = coefficient, 
                           color = polynomial, group = polynomial)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = polynomial), 
                alpha = 0.2, color = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
    labs(title = title,
         subtitle = paste("Estimeret effekt ved forskellige båndbredder og polynomiske grader"),
         x = "Båndbredde",
         y = "Estimeret behandlingseffekt",
         color = "Polynomisk grad",
         fill = "Polynomisk grad") +
    scale_color_manual(values = poly_colors) +
    scale_fill_manual(values = poly_colors) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 12),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "lightgray")
    )
  
  signifikante_resultater <- results %>%
    dplyr::filter(p_value < 0.05) %>%
    dplyr::group_by(polynomial) %>%
    dplyr::summarize(
      min_bw = min(bandwidth, na.rm = TRUE),
      max_bw = max(bandwidth, na.rm = TRUE),
  )

  
  if (nrow(signifikante_resultater) > 0) {
    for (i in 1:nrow(signifikante_resultater)) {
      poly_level <- signifikante_resultater$polynomial[i]
      min_bw <- signifikante_resultater$min_bw[i]
      max_bw <- signifikante_resultater$max_bw[i]
          
      }
    }
  
  
  return(list(plot = p, sig_intervals = signifikante_resultater))
}  


# Funktion til at sammenligne flere modeller i samme coefficient stability plot
compare_coefficient_stability <- function(data_list, 
                                          labels,
                                          outcome_var = "miljø_afhængig",
                                          running_var = "centered_lagged_pervote_samlet", 
                                          bw_list = seq(1, 3, by = 0.5),
                                          polynomial = 1,  # Fixed polynomial for comparison
                                          cutpoint = 0,
                                          covariates = NULL,
                                          title = "Comparison of Models") {
  
  # Tom dataframe til resultater
  all_results <- data.frame()
  
  # Loop over hver model/dataset
  for (i in 1:length(data_list)) {
    current_data <- data_list[[i]]
    current_label <- labels[i]
    
    # Håndtering af kovariater
    covs_matrix <- NULL
    if (!is.null(covariates)) {
      # Forsøg at bygge kovariat-matrix
      covs_list <- list()
      for (j in 1:length(covariates)) {
        col_name <- covariates[j]
        if (col_name %in% names(current_data)) {
          covs_list[[j]] <- current_data[[col_name]]
        } else {
          warning(paste("Kolonne", col_name, "findes ikke i datasættet", i))
          covs_list[[j]] <- NULL
        }
      }
      # Fjern NULL værdier
      covs_list <- covs_list[!sapply(covs_list, is.null)]
      
      # Hvis der er nogen kovariater tilbage, lav en matrix
      if (length(covs_list) > 0) {
        covs_matrix <- do.call(cbind, covs_list)
      }
    }
    
    # Loop over båndbredder
    for (bw in bw_list) {
      # Prøv at køre rdrobust med fejlhåndtering
      tryCatch({
        # Kør rdrobust 
        if (is.null(covs_matrix)) {
          rd <- rdrobust(
            y = current_data[[outcome_var]],
            x = current_data[[running_var]],
            c = cutpoint,
            p = polynomial,
            h = bw
          )
        } else {
          rd <- rdrobust(
            y = current_data[[outcome_var]],
            x = current_data[[running_var]],
            c = cutpoint,
            p = polynomial,
            h = bw,
            covs = covs_matrix
          )
        }
        
        # Gem estimate og standardfejl
        temp_results <- data.frame(
          model = current_label,
          bandwidth = bw,
          coefficient = rd$coef[1],  # Conventional estimate
          se = rd$se[1],
          ci_lower = rd$coef[1] - 1.96 * rd$se[1],
          ci_upper = rd$coef[1] + 1.96 * rd$se[1],
          p_value = rd$pv[1]
        )
        
        all_results <- rbind(all_results, temp_results)
        
      }, error = function(e) {
        warning(paste("Fejl ved kørsel af rdrobust for model", i, "med bw =", bw, ":", e$message))
      })
    }
  }
  
  # Check om man fik resultater
  if (nrow(all_results) == 0) {
    stop("Ingen resultater blev produceret. Tjek data og parametre.")
  }
  
  # Konverter model til faktor
  all_results$model <- factor(all_results$model, levels = labels)
  
  # Opret farveskala
  model_colors <- c("blue", "red", "green", "purple", "orange", 
                    "cyan", "magenta", "brown", "pink", "gray")[1:length(labels)]
  
  # Lav plot
  p <- ggplot(all_results, aes(x = bandwidth, y = coefficient, 
                               color = model, group = model)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = model), 
                alpha = 0.15, color = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
    labs(title = title,
         subtitle = paste0("Sammenligning med polynomisk grad ", polynomial),
         x = "Båndbredde",
         y = "Estimeret behandlingseffekt",
         color = "Model",
         fill = "Model") +
    scale_color_manual(values = model_colors) +
    scale_fill_manual(values = model_colors) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 12),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "lightgray")
    )
  
  # Find signifikante båndbredder for hver model
  signifikante_resultater <- all_results %>%
    dplyr::filter(p_value < 0.05) %>%
    dplyr::group_by(model) %>%
    dplyr::summarize(
      min_bw = min(bandwidth, na.rm = TRUE),
      max_bw = max(bandwidth, na.rm = TRUE)
    )
  
  # Hvis der er signifikante resultater, tilføj semitransparente rektangler
  if (nrow(signifikante_resultater) > 0) {
    for (i in 1:nrow(signifikante_resultater)) {
      model_level <- signifikante_resultater$model[i]
      min_bw <- signifikante_resultater$min_bw[i]
      max_bw <- signifikante_resultater$max_bw[i]
      

      
    }
  }
  
  return(p)
}

# Dot-whisker plot for compact representation of multiple models
create_dotwhisker_plot <- function(data_list, 
                                   labels,
                                   outcome_var = "miljø_afhængig",
                                   running_var = "centered_lagged_pervote_samlet", 
                                   bw_values = c(1, 2, 3),  # Selected bandwidths
                                   polynomial = 1,
                                   cutpoint = 0,
                                   covariates = NULL,
                                   title = "Dot-whisker Plot of Model Estimates") {
  
  # Tom dataframe til resultater
  all_results <- data.frame()
  
  # Loop over hver model/dataset
  for (i in 1:length(data_list)) {
    current_data <- data_list[[i]]
    current_label <- labels[i]
    
    # Håndtering af kovariater
    covs_matrix <- NULL
    if (!is.null(covariates)) {
      # Forsøg at bygge kovariat-matrix
      covs_list <- list()
      for (j in 1:length(covariates)) {
        col_name <- covariates[j]
        if (col_name %in% names(current_data)) {
          covs_list[[j]] <- current_data[[col_name]]
        } else {
          warning(paste("Kolonne", col_name, "findes ikke i datasættet", i))
          covs_list[[j]] <- NULL
        }
      }
      # Fjern NULL værdier
      covs_list <- covs_list[!sapply(covs_list, is.null)]
      
      # Hvis der er nogen kovariater tilbage, lav en matrix
      if (length(covs_list) > 0) {
        covs_matrix <- do.call(cbind, covs_list)
      }
    }
    
    # Loop over udvalgte båndbredder
    for (bw in bw_values) {
      # Prøv at køre rdrobust med fejlhåndtering
      tryCatch({
        # Kør rdrobust 
        if (is.null(covs_matrix)) {
          rd <- rdrobust(
            y = current_data[[outcome_var]],
            x = current_data[[running_var]],
            c = cutpoint,
            p = polynomial,
            h = bw
          )
        } else {
          rd <- rdrobust(
            y = current_data[[outcome_var]],
            x = current_data[[running_var]],
            c = cutpoint,
            p = polynomial,
            h = bw,
            covs = covs_matrix
          )
        }
        
        # Gem estimate og standardfejl
        temp_results <- data.frame(
          model = current_label,
          bandwidth = paste0("h = ", bw),
          coefficient = rd$coef[1],
          se = rd$se[1],
          ci_lower = rd$coef[1] - 1.96 * rd$se[1],
          ci_upper = rd$coef[1] + 1.96 * rd$se[1],
          p_value = rd$pv[1]
        )
        
        all_results <- rbind(all_results, temp_results)
        
      }, error = function(e) {
        warning(paste("Fejl ved kørsel af rdrobust for model", i, "med bw =", bw, ":", e$message))
      })
    }
  }
  
  # Check om man fik resultater
  if (nrow(all_results) == 0) {
    stop("Ingen resultater blev produceret. Tjek data og parametre.")
  }
  
  # Konverter til faktorer for bedre visualisering
  all_results$model <- factor(all_results$model, levels = labels)
  all_results$bandwidth <- factor(all_results$bandwidth, 
                                  levels = paste0("h = ", bw_values))
  
  # Create a grouping variable for faceting
  all_results$group <- paste(all_results$model, all_results$bandwidth, sep = " / ")
  all_results$group <- factor(all_results$group, 
                              levels = unique(all_results$group)[order(all_results$model, all_results$bandwidth)])
  
  # Markér signifikante estimater
  all_results$significant <- all_results$p_value < 0.05
  
  # Lav plot
  p <- ggplot(all_results, aes(x = model, y = coefficient)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper, color = bandwidth, 
                        shape = significant),
                    position = position_dodge(width = 0.5), size = 0.7) +
    labs(title = title,
         subtitle = paste0("Polynomisk grad = ", polynomial),
         x = NULL,
         y = "Estimeret behandlingseffekt",
         color = "Båndbredde",
         shape = "p < 0.05") +
    coord_flip() +
    scale_shape_manual(values = c(1, 16)) +  # Open circle for non-significant, filled for significant
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 12),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "lightgray")
    )
  
  return(p)
}