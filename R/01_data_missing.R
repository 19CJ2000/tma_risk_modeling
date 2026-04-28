############################################################
# MISSING DATA EXPLORATION 
############################################################
# ----- PACKAGES -----
library(here)      # relative file importing
library(tidyverse) # syntax simplicity 
library(naniar)    # missing data exploration 
library(DT)        # output table creation 


# ----- DATA IMPORT -----
raw <- read_excel(here("data", "raw.xlsx"),
                  skip = 1,  # remove label row from excel file
                  na = ""    # blank cells are taken as NA
                  )


#  ----- PERCENT MISSINGNESS -----
pct_miss_case(raw)    # % of rows with at least one NA value (only 5.56%, CCA justified)
pct_miss(raw)         # % of dataset that is NA 


# ----- TABULAR OVERVIEW -----
miss_var_summary(raw) %>%  arrange(desc(pct_miss)) %>%
  mutate(
    pct_miss = round(pct_miss, 2),
    variable = recode(variable,
                      racegp = "racegp (Patient Race)",
                      priautogp = "priautogp (Prior autologous transplant)",
                      kidfunc = "kidfunc (Baseline Kidney Fucntion)",
                      condgp = "condgp (Conditioning Regimen)"
    )) %>%
  rename(
    "Variable (description)" = variable,
    "Missing (n)" = n_miss,
    "Missing (%)" = pct_miss
  ) %>%
  datatable(             # from DT package (interactable tables)
    options = list(
      pageLength = 10,    # make table only 5 rows long 
      scrollY = TRUE,    # enable horizontal scrolling 
      paging = TRUE,     # shows page count
      info = TRUE        # showing x to y of z entries
      
    )) 


#  ----- VISUAL OVERIVIEWS -----
# Plot 1: Bar Chart 
gg_miss_var(raw) +
  labs(title = "Missingness Across the Combined Dataset")

# Plot 2: Heat Map 
vis_miss(raw) +
  labs(title = "Heat Map of Missingness Across the Combined Dataset") +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8, angle = 90, hjust = 0.1)
  )

# Plot 3: Upset Plot
gg_miss_upset(raw)



