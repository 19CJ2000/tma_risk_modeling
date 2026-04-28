############################################################
# CUMULATIVE INCIDENCE FUNCTION (CIF) PLOTS OF TMA DATA 
#    CIF vs KM: KM curves assume independent censoring / dont account for competing risk of death
#    F&G vs SRF: 
#       Fine-Gray assumes a parametric form (subdistribution hazard)
#       RSF learns patterns from the data without assuming a hazard model
############################################################
# ----- PACKAGES -----
library(here)       # relative file importing
library(tidyverse)  # syntax simplicity 
library(cmprsk)     # competing risks survival analysis
library(gt)         # output table building

# ----- DATA IMPORT -----
tma_clean <- readRDS(here("data", "processed.rds"))


# ----- RAW TA-TMA INCIDENCE -----
tma_inc <- tma_clean %>%
  summarise(
    n_total = n(),
    n_tma = sum(compete_event == 1, na.rm = TRUE),
    n_no_tma = sum(compete_event %in% c(0,2), na.rm = TRUE),
    incidence = n_tma / n_total * 100
  )
tma_inc  # rare ~ 3% 



# ----- CIF PLOTS FOR F&G MODEL -----

# CIF PLOT 1
cif_tma <- cmprsk::cuminc(           # data prep 
  ftime = tma_clean$compete_time,   
  fstatus = tma_clean$compete_event, # 0=censor, 1=TMA, 2=death
  cencode = 0  
)


# CIF Plot 1.1: base plot 
png("outputs/figures/cif_plot1_full.png", width = 2000, height = 1600, res = 300) # open graphics
plot(cif_tma, 
     main = "Cumulative Incidence of TA-TMA \n(Accounting for Competing Risks)",
     xlab = "Time from transplant (months)",
     ylab = "Cumulative Incidence",
     curvlab = c("CIF of TMA for those at risk", "CIF of Death for those at risk"),
     col = c("red", "blue"),
     lwd = 2)
dev.off()  # close graphics 


# CIF Plot 1.2: zoomed
png("outputs/figures/cif_plot1_zoomed.png", width = 2000, height = 1600, res = 300)
plot(cif_tma, 
     main = "Cumulative Incidence of TA-TMA \n(Accounting for Competing Risks)",
     xlab = "Time from transplant (months)",
     ylab = "Cumulative Incidence",
     curvlab = c("CIF of TMA for those at risk", "CIF of Death for those at risk"),
     col = c("red", "blue"),
     lwd = 2,
     xlim = c(0, 24),
     ylim = c(0, 0.05),
     yaxt = "n")
axis(2, at = seq(0, 0.1, by = 0.01), las = 1) # Y-axis ticks
grid(nx = 25, ny = 30, col = "gray92", lty = "dotted") # Minor grid lines
grid(nx = 6, ny = 6, col = "gray80", lty = "solid") # Major grid lines
dev.off()



# CIF PLOT 2 (KIDNEY FUNCTION STRATIFIED)
cif_by_sex <- cmprsk::cuminc(         # data prep 
  ftime   = tma_clean$compete_time,   
  fstatus = tma_clean$compete_event,  
  group   = tma_clean$kidfunc,
  cencode = 0
)

# CIF Plot 2.1: base plot 
png("outputs/figures/cif_plot2_full.png", width = 2000, height = 1600, res = 300)
plot(cif_by_sex,   
     main = "Cumulative Incidence of TA-TMA \n(Accounting for Competing Risks) By Kidney Function",
     xlab = "Time from transplant (months)",
     ylab = "Cumulative Incidence",
     curvlab = c("CIF of TMA (Normal Kidney Function)", 
                 "CIF of TMA (Decreased Kidney Function)",
                 "CIF of Death (Normal Kidney Function)",
                 "CIF of Death (Decreased Kidney Function)"
     ),
     col = c("red", "blue", "green", "orange"),# Green/Orange = Death, Red/Blue = TMA
     lwd = 2)
dev.off()

# CIF Plot 2.2: zoomed 
png("outputs/figures/cif_plot2_zoomed.png", width = 2000, height = 1600, res = 300)
plot(cif_by_sex,  
     main = "Cumulative Incidence of TA-TMA \n(Accounting for Competing Risks) By Kidney Function",
     xlab = "Time from transplant (months)",
     ylab = "Cumulative Incidence",
     curvlab = c("CIF of TMA (Normal Kidney Function)", 
                 "CIF of TMA (Decreased Kidney Function)",
                 "CIF of Death (Normal Kidney Function)",
                 "CIF of Death (Decreased Kidney Function)"
     ),
     col = c("red", "blue", "green", "orange"), 
     lwd = 2,
     xlim = c(0, 24),
     ylim = c(0, 0.05),
     yaxt = "n")
axis(2, at = seq(0, 0.1, by = 0.01), las = 1) # Y-axis ticks
grid(nx = 25, ny = 30, col = "gray92", lty = "dotted") # Minor grid lines
grid(nx = 6, ny = 6, col = "gray80", lty = "solid") # Major grid lines
dev.off()



# ----- CIF PLOTS FROM SRF MODEL -----

srf_model <- readRDS(here("outputs/models", "srf_model.rds")) # import RSF model from 05_model_SRF.r 

cif <- predict(srf_model, newdata = tma_clean, cause = 1) # Examine CIF 
cif_tma <- cif$cif[, , 1]  # select 3rd dim = 1, because 1 = Cause = TMA
dim(cif_tma)               # should be 21825 x 150 (reduced from 3 to 2 dimensions)


# CIF Plot 3: first 10 patients using SRF data 
png("outputs/figures/cif_plot3_srf.png", width = 2000, height = 1600, res = 300)
matplot(
  x = cif$time.interest,   # length 150
  y = t(cif_tma[1:10, ]),  # now 150 x 10
  type = "l",
  lty = 1:10,
  col = 1:10,
  main = "Cumulative Incidence of TA-TMA \n(Accounting for Competing Risks), Using SRF Model",
  xlab = "Time",
  ylab = "Cumulative incidence of TMA"
)
legend("topright", legend = paste("Patient", 1:10), col = 1:10, lty = 1:10)
dev.off()

# EXPLANATION: CIF gives the probability of a specific event happening over time, accounting for competing risks.
#              In classical survival analysis, we often estimate one CIF for the whole population.
#              In RSF, each patient has their own predicted CIF, because RSF produces individualized predictions based on the patient’s covariates.




