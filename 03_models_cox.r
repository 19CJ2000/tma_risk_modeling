############################################################
# COX MODELS (REPLICATION + EXTENSION OF ORIGINAL STUDY)
   # Cause-specific Cox model for TA-TMA
   # Competing risk of death handled separately via Fine & Gray
############################################################
# ----- PACKAGES -----
library(here)      # relative file importing
library(tidyverse) # syntax simplicity 
library(survival)  # survival analysis 
library(gt)        # output table building


# ----- DATA IMPORT -----
tma_clean <- readRDS(here("data", "processed.rds"))


# ----- UNIVARIATE COX MODELS -----

# Select covariates
cox_vars <- c(
  "sex",
  "racegp",
  "atgcampathgp",
  "condgp",
  "diseasegp",
  "kidfunc",
  "priautogp",
  "gvhd_gp"
)


# Build univariate model function
fit_cox <- function(var, df) {
  
  form <- as.formula(paste0("Surv(compete_time, tma_event) ~ ", var))
  
  model <- coxph(form, data = df, ties = "efron")
  
  list(
    model = model,
    ph_test = cox.zph(model)
  )
}


# Fit all univariate models
cox_models <- setNames(
  map(cox_vars, ~ fit_cox(.x, tma_clean)),
  cox_vars
)


# PH assumption checks
ph_results <- map(cox_models, function(x) {
  ph <- cox.zph(x$model)
  as.data.frame(ph$table)
})
ph_results


# Build extraction function for univariate model results 
extract_cox <- function(model_obj, var_name) {
  
  m <- model_obj$model
  s <- summary(m)
  
  hr <- exp(coef(m))
  ci <- exp(confint(m))
  
  # Clean level names properly
  clean_levels <- names(hr) %>%
    stringr::str_remove(paste0("^", var_name)) %>%  # remove predictor name
    stringr::str_remove("^=") %>%                   # safety cleanup
    stringr::str_trim()
  
  tibble(
    model = var_name,
    level = clean_levels,
    HR = as.numeric(hr),
    CI_lower = ci[, 1],
    CI_upper = ci[, 2],
    p_value = s$coefficients[, "Pr(>|z|)"]
  )
}


# Build univariate cox model output table
univ_results <- imap_dfr(cox_models, extract_cox)    # Store univariate results in tibble
univ_results

uni_cox_output <- univ_results %>%
  mutate(
    p_value = case_when(   # ensure no p-vals "0.0000"
      p_value < 0.001 ~ "<0.001",
      TRUE ~ format(round(p_value, 3), nsmall = 3)
      )
    ) %>%
  gt() %>%
  # Header / Title
  tab_header(
    title = "Univariate Cox Proportional Hazards Model Results"
    ) %>%
  fmt_number(columns = c(HR, CI_lower, CI_upper), decimals = 2) %>% # rounding 
  fmt_number(columns = p_value, decimals = 4) %>%  # rounding 
  # Column Renaming 
  cols_label(
    model = "Model",
    level = "Level (Ref excluded)",
    HR = "Hazard Ratio",
    CI_lower = "CI Lower",
    CI_upper = "CI Upper",
    p_value = "p-value"
    ) %>%
  # Bold p-value Formating 
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = p_value,
      rows = p_value < 0.05
      )
    ) %>% 
  # Footer 
  tab_source_note(
    source_note = "Reference categories: sex = Female; racegp = White; atgcampathgp = No ATG or CAMPATH; condgp = Myeloablative TBI; diseasegp = AML; kidfunc = Normal kidney function; priautogp = No prior autologous HCT; gvhd_gp = Post-cy +/- others."
  )
uni_cox_output

# Save final uni_cox_output table
gtsave(uni_cox_output, here("outputs/tables", "univar_cox_output.html"))



# ----- MULTIVARIABLE COX MODEL -----
# Fit model
full_model <- coxph(
  Surv(compete_time, tma_event) ~ 
    sex +
    atgcampathgp +
    condgp +
    racegp +
    diseasegp +
    kidfunc +
    priautogp +
    strata(gvhd_gp),  # GVHD prophylaxis violates proportional hazards assumption, thus included as strata
  data = tma_clean,
  ties = "efron"
)


# PH assumption check
cox.zph(full_model)  # not significant overall 


# Build multivariable cox model output table
multi_summary <- summary(full_model) # Store multivariable results in tibble

multi_results <- tibble(
  variable = rownames(multi_summary$coefficients),
  HR = exp(multi_summary$coefficients[, "coef"]),
  CI_lower = exp(confint(full_model)[, 1]),
  CI_upper = exp(confint(full_model)[, 2]),
  p_value = multi_summary$coefficients[, "Pr(>|z|)"]
)
multi_results

multi_cox_output <- multi_results %>%
  mutate(
    p_value = case_when(
      p_value < 0.001 ~ "<0.001",
      TRUE ~ format(round(p_value, 3), nsmall = 3)
    )) %>% 
  gt() %>%
  tab_header(
    title = "Multivariable Cox Proportional Hazards Model Results"
  ) %>%
  fmt_number(
    columns = c(HR, CI_lower, CI_upper),
    decimals = 2
  ) %>%
  cols_label(
    variable = "CovariateLevel",
    HR = "Hazard Ratio",
    CI_lower = "CI Lower",
    CI_upper = "CI Upper",
    p_value = "p-value"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = p_value,
      rows = multi_results$p_value < 0.05
    )
  ) %>%
  tab_source_note(
    source_note = "GVHD prophylaxis (gvhd_gp) was included as a stratification variable due to violation of proportional hazards assumption; therefore no hazard ratio is estimated for this covariate. Reference categories: sex = Female; racegp = White; atgcampathgp = No ATG or CAMPATH; condgp = Myeloablative TBI; diseasegp = AML; kidfunc = Normal kidney function; priautogp = No prior autologous HCT; gvhd_gp = Post-cy +/- others."
  )
multi_cox_output

# Save final multi_cox_output table
gtsave(multi_cox_output, here("outputs/tables", "multi_cox_output.html"))

