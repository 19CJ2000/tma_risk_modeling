############################################################
# SURVIVAL RANDOM FOREST (SRF) MODEL FOR TMA DATA 
############################################################
# ----- PACKAGES -----
library(here)              # relative file importing
library(tidyverse)         # syntax simplicity 
library(survival)          # survival analysis 
library(randomForestSRC)   # survival random forests
library(gt)                # output table building


# ----- DATA IMPORT -----
tma_clean <- readRDS(here("data", "processed.rds"))


# ----- SRF MODEL FIT -----

set.seed(416)  # reproducibility

srf_model <- rfsrc(
  formula = Surv(compete_time, compete_event) ~ sex + atgcampathgp + condgp + racegp +
    diseasegp + kidfunc + gvhd_gp + priautogp,
  data = tma_clean,
  ntree = 150,             # number of trees 
  importance = TRUE,       # for VIMP 
  mtry = 3,         
  nodesize = 15,           # minimum number of patients per leaf node 
  splitrule = "logrankCR", # competing risk (death)
  cause = 1                # 1 = TMA = event of interest 
)

print(srf_model)

saveRDS(srf_model, here("outputs/models", "srf_model.rds"))  # import for CIF plot in 06_CIF_plots.r 


# ----- OUT-OF-BAG (OOB) PREDICTION ERROR -----
tail(srf_model$err.rate, 1)   # tail = final OOBs 
#    TMA: ~ ~0.38 # moderate
#    Death: ~0.44 # weak


# ----- VIMP -----

# Extract variable importance matrix
vimp <- srf_model$importance

# Assign variable names as rownames
rownames(vimp) <- colnames(srf_model$xvar)

# Convert to data frame
vimp_df <- as.data.frame(vimp) %>%
  tibble::rownames_to_column("Variable") %>%
  tidyr::pivot_longer(
    cols = everything()[-1],       # all columns except Variable
    names_to = "Event",    # event.1 = TMA, event.2 = competing event 
    values_to = "VIMP"
  ) %>%
  dplyr::arrange(desc(VIMP))

vimp_df
# For Event 1 (TMA) 
   # GVHD: Most important predictor of TMA. Changes in GVHD treatment/group affect TMA risk.
   # Prior autologous transplant strongly predicts TMA

# For Event 2 (Competing) 
   # Disease type is most important for the competing event (not TMA).

# Some variables (like racegp) have low or even slightly negative VIMP, meaning they contribute little.



# -----TUNE HYPERPARAMTERES -----  
tune <- tune.rfsrc(
  Surv(compete_time, compete_event) ~ sex + atgcampathgp + condgp + racegp + diseasegp + kidfunc + gvhd_gp + priautogp,
  data = tma_clean,
  ntreeTry = 75,   
  stepFactor = 1.5,
  improve = 0.01,
  nodesizeTry = c(15, 20),
  mtryStart = 3,
  nsplitTry = 5
)
print(tune)


# Apply tuned choices 
srf_model_tuned <- rfsrc(
  Surv(compete_time, compete_event) ~ sex + atgcampathgp + condgp + racegp +
    diseasegp + kidfunc + gvhd_gp + priautogp,
  data = tma_clean,
  ntree = 150,
  importance = TRUE,
  mtry = tune$optimal["mtry"],
  nodesize = tune$optimal["nodesize"],
  splitrule = "logrankCR",
  cause = 1
)

# Compare original vs tuned model
srf_model$err.rate
srf_model_tuned$err.rate


