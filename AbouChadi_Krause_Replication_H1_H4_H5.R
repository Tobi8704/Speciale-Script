# =====================================================================
# H1, H4, H5 -- Direct replication of Abou-Chadi & Krause (2020, BJPS)
# =====================================================================
# Mirrors their published method on our own data as closely as possible,
# for all three hypotheses. Outcome = miljø_afhængig (mainstream party's
# environmental position, a MARPOR log-ratio built the same way as their
# radical-right contagion variable), running variable =
# centered_lagged_pervote_samlet (Green party's lagged vote share minus the
# country's electoral threshold), treatment = lagged_i_parlament (Green
# party in parliament; sharp cutoff at 0).
#
# Nothing added beyond their specification: IK bandwidth, triangular kernel
# weighting, local linear (p = 1) regression with country and election-date
# fixed effects, two-way clustered SEs by party and election date (their
# actual clust1/clust2, confirmed from their published source code
# rrp_rdd.R / rrp_rdd_functions.R). No robustness battery, no bandwidth
# sensitivity check, no bias-corrected CI -- this script reports exactly
# their Table 2/3/4-style estimates and nothing more, run end to end.

pacman::p_load(dplyr, rdd, fixest)

final_dataset <- readRDS("final_dataset_publisering.rds")

# ---- Shared estimator: exactly Abou-Chadi & Krause's specification ----
# IK bandwidth (rdd::IKbandwidth), triangular kernel weights (rdd::kernelwts),
# local linear regression with country + election-date FE, two-way clustered
# SEs by party and edate (their clust1/clust2). Returns one row: their
# Table-style estimate for whichever subgroup is passed in.
ack_replicate <- function(data, label) {
  bw <- IKbandwidth(
    X = data$centered_lagged_pervote_samlet,
    Y = data$miljø_afhængig,
    cutpoint = 0
  )

  data$kernel_weight <- kernelwts(
    X = data$centered_lagged_pervote_samlet,
    center = 0,
    bw = bw,
    kernel = "triangular"
  )
  data <- data[data$kernel_weight > 0, ]

  # Sharp design: lagged_i_parlament is a deterministic function of the
  # running variable's sign, so their ivreg() (sharp treatment as its own
  # instrument) is numerically identical to weighted OLS -- feols()
  # reproduces the same point estimate and SE given the same weights and
  # clustering.
  model <- feols(
    miljø_afhængig ~ lagged_i_parlament + centered_lagged_pervote_samlet +
      lagged_i_parlament:centered_lagged_pervote_samlet | country + edate,
    data = data,
    weights = ~kernel_weight,
    cluster = ~party + edate
  )

  cat("\n====", label, "====\n")
  print(summary(model))

  data.frame(
    Group     = label,
    Bandwidth = round(bw, 3),
    LATE      = round(coef(model)["lagged_i_parlament"], 3),
    StdError  = round(fixest::se(model)["lagged_i_parlament"], 3),
    pvalue    = round(fixest::pvalue(model)["lagged_i_parlament"], 4),
    N_below_c = sum(data$centered_lagged_pervote_samlet < 0),
    N_above_c = sum(data$centered_lagged_pervote_samlet >= 0)
  )
}

# =====================================================================
# H1 -- Pooled effect across all mainstream parties
# =====================================================================
h1_data <- final_dataset[complete.cases(final_dataset[, c(
  "miljø_afhængig", "centered_lagged_pervote_samlet", "lagged_i_parlament",
  "country", "edate", "party"
)]), ]

h1_result <- ack_replicate(h1_data, "H1 (Pooled)")

# =====================================================================
# H4 -- Moderation by mainstream party ideology (left vs. right)
# =====================================================================
h4_data <- final_dataset[complete.cases(final_dataset[, c(
  "miljø_afhængig", "centered_lagged_pervote_samlet", "lagged_i_parlament",
  "left_party_lag1", "country", "edate", "party"
)]), ]

h4_right <- subset(h4_data, left_party_lag1 == 0)
h4_left  <- subset(h4_data, left_party_lag1 == 1)

h4_right_result <- ack_replicate(h4_right, "H4 (Right-Wing)")
h4_left_result  <- ack_replicate(h4_left,  "H4 (Left-Wing)")

# =====================================================================
# H5 -- Moderation by mainstream party size (large vs. small)
# =====================================================================
h5_data <- final_dataset[complete.cases(final_dataset[, c(
  "miljø_afhængig", "centered_lagged_pervote_samlet", "lagged_i_parlament",
  "large_party_lag1", "country", "edate", "party"
)]), ]

h5_large <- subset(h5_data, large_party_lag1 == 1)
h5_small <- subset(h5_data, large_party_lag1 == 0)

h5_large_result <- ack_replicate(h5_large, "H5 (Large Parties)")
h5_small_result <- ack_replicate(h5_small, "H5 (Small Parties)")

# =====================================================================
# Combined summary table
# =====================================================================
ack_results <- rbind(
  h1_result,
  h4_right_result, h4_left_result,
  h5_large_result, h5_small_result
)

cat("\n\n==== Abou-Chadi & Krause-style replication: all results ====\n")
print(ack_results)

dir.create("Publisering Analyse/AC&K Replication", showWarnings = FALSE)
write.csv(ack_results, "Publisering Analyse/AC&K Replication/ack_replication_results.csv", row.names = FALSE)
