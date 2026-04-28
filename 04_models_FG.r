############################################################
# FINE-GRAY MODELS (COMPETING RISKS FRAMEWORK) FOR TMA DATA 
############################################################
# ----- PACKAGES -----
library(here)       # relative file importing
library(tidyverse)  # syntax simplicity 
library(survival)   # survival analysis 
library(cmprsk)     # competing risks survival analysis
library(gt)         # output table building

# ----- DATA IMPORT -----
tma_clean <- readRDS(here("data", "processed.rds"))



# ----- DESIGN MATRIX -----
X <- model.matrix(
  ~ sex + atgcampathgp + condgp + racegp +
    diseasegp + kidfunc + gvhd_gp + priautogp,
  data = tma_clean
)[, -1]  # remove intercept



# ----- MULTIVARIABLE FINE–GRAY MODEL -----
fg_model <- cmprsk::crr(
  ftime   = tma_clean$compete_time,
  fstatus = tma_clean$compete_event,
  cov1    = X
)
summary(fg_model)



# Build multivariable fine and gray model output table
fg_results <- tibble(
  variable = rownames(summary(fg_model)$coef),
  sHR      = summary(fg_model)$coef[, "exp(coef)"],
  CI_lower = summary(fg_model)$conf.int[, "2.5%"],
  CI_upper = summary(fg_model)$conf.int[, "97.5%"],
  p_value  = summary(fg_model)$coef[, "p-value"]
)

fg_output <- fg_results %>%
  mutate(
    p_value = case_when(
      p_value < 0.001 ~ "<0.001",
      TRUE ~ format(round(p_value, 3), nsmall = 3)
    )) %>%
  gt() %>%
  tab_header(
    title = "Multivariable Fine-Gray Competing Risks Regression Results"
    ) %>%
  fmt_number(columns = c(sHR, CI_lower, CI_upper), decimals = 2) %>% # round
  # Fix names 
  cols_label(
    variable = "Covariate Level",
    sHR      = "Subdistribution HR",
    CI_lower = "CI Lower",
    CI_upper = "CI Upper",
    p_value  = "p-value"
  ) %>%
  # Bold significant p-values 
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = p_value,
      rows    = fg_results$p_value < 0.05
    )
  ) %>%
  tab_source_note(
    source_note = "Reference categories: sex = Female; racegp = White; atgcampathgp = No ATG or CAMPATH; condgp = Myeloablative non-TBI; diseasegp = ALL; kidfunc = Normal kidney function; priautogp = No prior autologous HCT; gvhd_gp = Post-cy +/- others."
  )
fg_output

gtsave(fg_output, here("outputs/tables", "fg_output.html"))



# ----- MODEL DIAGNOSTICS -----

# Event distribution
table(tma_clean$compete_event)
table(tma_clean$sex,       tma_clean$compete_event)
table(tma_clean$racegp,    tma_clean$compete_event)
table(tma_clean$diseasegp, tma_clean$compete_event)
table(tma_clean$kidfunc,   tma_clean$compete_event)
table(tma_clean$gvhd_gp,   tma_clean$compete_event)

# Events per parameter
n_events <- sum(tma_clean$compete_event == 1)
n_params <- ncol(X)
n_events / n_params

# Subdistribution PH assumption check (manual) 
X_time <- X * log(tma_clean$compete_time)

fg_time_check <- cmprsk::crr(
  ftime   = tma_clean$compete_time,
  fstatus = tma_clean$compete_event,
  cov1    = cbind(X, X_time)
)

summary(fg_time_check)

# Linearity assumption: Not formally assessed (all covariates are categorical)
