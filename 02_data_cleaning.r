############################################################
# IMPORT AND CLEANING OF RAW TMA DATA 
############################################################
# ----- PACKAGES -----
library(here)      # relative file importing
library(readxl)    # excel file importing
library(tidyverse) # syntax simplicity 
library(naniar)    # missing data exploration 


# ----- DATA IMPORT -----
raw <- read_excel(here("data", "raw.xlsx"),
                  skip = 1,  # remove label row from excel file
                  na = ""    # blank cells are taken as NA
                 )


# ----- DATA GLANCE -----
glimpse(raw)
str(raw)           


# ----- DATA CLEANING -----

# Select variables of interest
processed <- raw %>% 
  select(tma_event, intxtma, intxsurv, dead, sex, atgcampathgp, condgp, 
         racegp, diseasegp, kidfunc, gvhd_gp, priautogp)

# Transform them
processed <- processed %>% 
  mutate(
    
    tma_event = case_when(  
      tma_event == "1" ~ 1,
      tma_event == "0" ~ 0,
      tma_event == "99" ~ NA, # convert 99's to NA 
      TRUE ~ NA),
    
    tma_time = case_when(
      tma_event == 1 ~ intxtma,
      tma_event == 0 ~ Inf,   #  TMA event does not occur within the observable time window.
      TRUE ~ NA_real_
    ),
    
    surv_time = intxsurv,  # time of death or last follow up
    
    compete_time = pmin(tma_time, surv_time, na.rm = TRUE),
    
    compete_event = case_when(                      # Will be used for F&G models (assumes non-informative censoring except for death as competing event)
      tma_event == 1 & tma_time <= surv_time ~ 1,   # TMA
      dead == 1 & surv_time < tma_time       ~ 2,   # Death before TMA
      TRUE                                   ~ 0    # Censored
    ),
      
    sex = as.factor(case_when(
      sex == "1" ~ "Male",
      sex == "2" ~ "Female",
      sex == "8" ~ NA, 
      sex == "9" ~ NA,
      is.na(sex) ~ NA,
      TRUE ~ NA)),
    
    atgcampathgp = factor(
      case_when(
        atgcampathgp %in% c("1", "2") ~ "ATG +/- CAMPATH",
        atgcampathgp == "3" ~ "CAMPATH alone",
        atgcampathgp == "4" ~ "No ATG or CAMPATH",
        TRUE ~ NA
      ),
      levels = c("No ATG or CAMPATH", "ATG +/- CAMPATH", "CAMPATH alone")), 
      
    condgp = as.factor(case_when(
      condgp == "1" ~ "Myeloablative TBI",
      condgp == "2" ~ "Myeloablative Bu based",
      condgp == "3" ~ "Other myeloablative",
      condgp == "4" ~ "RIC/NMA",
      condgp == "99" ~ NA,
      is.na(condgp) ~ NA,
      TRUE ~ NA)),
    
    racegp = factor(case_when(
      racegp == "1" ~ "White",
      racegp == "2" ~ "Black",
      racegp == "3" ~ "Other",
      TRUE ~ NA), 
      levels = c("White", "Black", "Other")),
    
    diseasegp = factor(case_when(
      diseasegp == "1" ~ "AML",
      diseasegp == "2" ~ "ALL",
      diseasegp == "3" ~ "CLL",
      diseasegp == "4" ~ "CML",
      diseasegp == "5" ~ "MDS",
      diseasegp == "6" ~ "MPS",
      diseasegp == "7" ~ "Other acute leukemia",
      diseasegp == "8" ~ "Lymphoma",
      diseasegp == "9" ~ "Multiple myeloma/PCD",
      diseasegp == "10" ~ "Other malignant diseases",
      diseasegp == "11" ~ "Immune disorders",
      diseasegp == "12" ~ "Inborn errors of metabolism",
      diseasegp == "13" ~ "Aplastic Anemia",
      diseasegp == "14" ~ "Hemoglobinopathy",
      diseasegp == "15" ~ "Inherited bone marrow failure",
      diseasegp == "16" ~ "Other nonmalignant disorders",
      diseasegp == "17" ~ "Other diseases",
      TRUE ~ NA)),
      
    kidfunc = factor(case_when(
      kidfunc == "0" ~ "Decreased kidney function",
      kidfunc == "1" ~ "Normal kidney function",
      TRUE ~ NA), 
      levels = c("Normal kidney function", "Decreased kidney function")),
    
    gvhd_gp = factor(case_when(
      gvhd_gp == "1" ~ "Post-cy +- others",
      gvhd_gp == "2" ~ "CNI + siro (no ptcy)",
      gvhd_gp == "3" ~ "Siro (without CNI or ptcy)",
      gvhd_gp == "4" ~ "Tac (no siro or ptcy)",
      gvhd_gp == "5" ~ "Csa (no tac, siro or ptcy)",
      gvhd_gp == "6" ~ "Other",
      TRUE ~ NA), levels = c("Post-cy +- others",
                             "Tac (no siro or ptcy)",
                             "Siro (without CNI or ptcy)",
                             "Csa (no tac, siro or ptcy)",
                             "CNI + siro (no ptcy)",
                             "Other")),
  
    priautogp = as.factor(case_when(
      priautogp == "0" ~ "No", 
      priautogp == "1" ~ "Yes", 
      priautogp == "99" ~ NA,
      TRUE ~ NA)),
    )

# Drop NAs (Complete Case Analysis) and duplicate/renamed variables
processed <- processed %>%  
  drop_na() %>% 
  select(-c(intxtma, intxsurv))    
  

# Quality Check
glimpse(processed)
str(processed)   

# Save cleaned data
saveRDS(processed, here("data", "processed.rds"))





