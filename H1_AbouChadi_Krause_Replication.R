# =====================================================================
# H1 -- Direct replication of Abou-Chadi & Krause (2020, BJPS) specification
# =====================================================================
# Mirrors their published method on our own data as closely as possible:
# outcome = miljø_afhængig (mainstream party's environmental position,
# a MARPOR log-ratio built the same way as their radical-right contagion
# variable), running variable = centered_lagged_pervote_samlet (Green
# party's lagged vote share minus the country's electoral threshold),
# treatment = lagged_i_parlament (Green party in parliament; sharp cutoff
# at 0).
#
# Nothing added beyond their specification: IK bandwidth, triangular kernel
# weighting, local linear (p = 1) regression with country and election-date
# fixed effects, two-way clustered SEs by party and election date (their
# actual clust1/clust2, confirmed from their published source code
# rrp_rdd.R / rrp_rdd_functions.R). No robustness battery, no bandwidth
# sensitivity check, no bias-corrected CI -- this script reports exactly
# their Table 2-style estimate and nothing more.

pacman::p_load(dplyr, rdd, fixest)

final_dataset <- readRDS("final_dataset_publisering.rds")

h1_data <- final_dataset[complete.cases(final_dataset[, c(
  "miljø_afhængig", "centered_lagged_pervote_samlet", "lagged_i_parlament",
  "country", "edate", "party"
)]), ]

# ---- 1. Bandwidth: Imbens-Kalyanaraman (2012), as AC&K use ----
# rdd::IKbandwidth() is what Abou-Chadi & Krause's own code calls. The 'rdd'
# package was archived by CRAN in 2025; if install.packages("rdd") fails,
# install the last CRAN release directly instead:
#   remotes::install_version("rdd", version = "0.57", repos = "http://cran.r-project.org")
bw <- IKbandwidth(
  X = h1_data$centered_lagged_pervote_samlet,
  Y = h1_data$miljø_afhængig,
  cutpoint = 0
)
cat("IK bandwidth:", round(bw, 3), "\n")

# ---- 2. Triangular kernel weights (their rdd::kernelwts(), kernel = "triangular") ----
h1_data$kernel_weight <- kernelwts(
  X = h1_data$centered_lagged_pervote_samlet,
  center = 0,
  bw = bw,
  kernel = "triangular"
)
h1_data <- h1_data[h1_data$kernel_weight > 0, ]

# ---- 3. Local linear regression (p = 1), country + election-date FE ----
# Sharp design: lagged_i_parlament is a deterministic function of the running
# variable's sign, so their ivreg() (with the sharp treatment as its own
# instrument) is numerically identical to weighted OLS -- feols() reproduces
# the same point estimate and SE given the same weights and clustering.
#
# Clustering deliberately kept as their clust1 = 'party', clust2 = 'edate'
# (not the country_election_id used elsewhere in our own pipeline to avoid
# cross-country date collisions) -- this script mirrors their code as-is,
# nothing improved on top of it.
h1_model <- feols(
  miljø_afhængig ~ lagged_i_parlament + centered_lagged_pervote_samlet +
    lagged_i_parlament:centered_lagged_pervote_samlet | country + edate,
  data = h1_data,
  weights = ~kernel_weight,
  cluster = ~party + edate
)

summary(h1_model)

h1_result <- data.frame(
  Bandwidth = round(bw, 3),
  LATE      = round(coef(h1_model)["lagged_i_parlament"], 3),
  StdError  = round(fixest::se(h1_model)["lagged_i_parlament"], 3),
  pvalue    = round(fixest::pvalue(h1_model)["lagged_i_parlament"], 4),
  N_below_c = sum(h1_data$centered_lagged_pervote_samlet < 0),
  N_above_c = sum(h1_data$centered_lagged_pervote_samlet >= 0)
)
print(h1_result)
