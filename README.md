# Competing Risks Analysis of Transplant-Associated Thrombotic Microangiopathy (TA-TMA)

This project re-analyses a large allogeneic hematopoietic cell transplant (allo-HCT) cohort (2008–2016) to evaluate risk factors for transplant-associated thrombotic microangiopathy (TA-TMA), with a focus on appropriate handling of competing risks.

The analysis extends a published study that primarily used Cox proportional hazards models by incorporating competing risks methodology and machine learning approaches.

---

## Objective

To investigate the incidence and risk factors of TA-TMA following allo-HCT, and to compare traditional survival analysis with competing risks and machine learning approaches.

Specifically, this project:

- Reproduces Cox proportional hazards models from the original study
- Implements competing risks analysis using Fine & Gray models
- Accounts for death without TA-TMA as a competing event
- Applies Survival Random Forests for non-parametric prediction
- Compares model outputs and variable importance across approaches

---

## Data

- Source: Published allo-HCT cohort (2008–2016)
- Sample size: ~23,000+ patients
- Population: Patients receiving allogeneic hematopoietic cell transplant

### Outcomes
- Primary outcome: Time to TA-TMA
- Competing event: Death without TA-TMA
- Censoring: last follow-up

---

## Key Variables

Clinical covariates included:

- Sex
- Race
- Disease type (e.g., AML, ALL, aplastic anemia)
- Conditioning regimen
- GVHD prophylaxis strategy
- Kidney function
- Prior autologous transplant
- Use of ATG/Campath

---

## Methods

### 1. Data Processing

- Recoding of clinical variables into analysis-ready factors
- Construction of:
  - Event indicator for TA-TMA
  - Competing risk indicator (death before TA-TMA)
- Missing data explored using graphical and tabular summaries

---

### 2. Classical Survival Analysis

- Cox proportional hazards models used to estimate risk factors for TA-TMA

---

### 3. Competing Risks Analysis

- Fine & Gray subdistribution hazard models used to account for competing risk of death

---

### 4. Machine Learning Approach

- Survival Random Forests used for non-parametric time-to-event modeling

- Model:
  - logrankCR split rule (accounts for competing risk)
  - cause-specific event modeling
- Hyperparameter Tuning
- Evaluation: 
  - prediction error (OOB)
  - variable importance (VIMP) 
  
---

## Key Findings (from modeling results)

- TA-TMA incidence rare (~3% over 3 years), consistent with prior literature
- Risk is strongly associated with:
  - GVHD-related factors
  - Kidney dysfunction
  - Prior autologous transplant
  - Conditioning intensity and immunosuppressive strategy

- Competing risks models show that ignoring death leads to biased incidence estimates

- Survival Random Forests identified similar key predictors to regression-based models, supporting robustness of findings

---

## Reproducibility

To reproduce this analysis:

```r
renv::restore()
```

---

## Software

- R
- survival
- cmprsk
- randomForestSRC
- tidyverse
- ggplot2
- here
- DT

---

## Notes

- All preprocessing steps are fully reproducible
- Competing risks framework is used throughout to avoid biased incidence estimation
- Machine learning model used for complementary predictive insight, not causal inference