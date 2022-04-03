library(tidyverse)
library(magrittr)

# Data import ----

data1 = read_csv("data/ADNIMERGE_May15.2014.csv")
data2 = read_csv("data/UPENN_CSF Biomarkers_baseline_May15.2014.csv")

data1 %<>% 
  filter(VISCODE == "bl")

## Combining both datasets and selecting relevant variables
dataset %<>% 
  inner_join(data1, data2, by = "RID") %>% 
  select(RID, FDG, ABETA, PTAU, APOE4, PTGENDER, AGE, PTEDUCAT, DX) %>%
  na.omit()

# Data description ----