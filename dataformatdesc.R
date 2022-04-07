library(tidyverse)
library(magrittr)
library(cluster)

# Data import ----

data1 = read_csv("data/ADNIMERGE_May15.2014.csv")
data2 = read_csv("data/UPENN_CSF Biomarkers_baseline_May15.2014.csv")

data1 %<>% 
  filter(VISCODE == "bl")

## Combining both datasets and selecting relevant variables
dataset <- inner_join(data1, data2, by = "RID") %>% 
  select(RID, FDG, ABETA, PTAU, APOE4, PTGENDER, AGE, PTEDUCAT, DX.bl) %>%
  na.omit() %>% 
  mutate(APOE4 = as.factor(APOE4),
         PTGENDER = as.factor(PTGENDER),
         DX = as.factor(DX.bl)) %>% 
  mutate(DX = fct_collapse(DX, MCI = c("EMCI", "LMCI"))) %>% 
  select(-DX.bl) 

# Saving dataset for further use
write_csv(dataset, "data/datasetADNI.csv")

# Data description ----

str(dataset)
summary(dataset)

dataset %>% 
  select_if(is.numeric) %>% 
  map(sd)

dataset %>% 
  select_if(is.factor) %>% 
  map(table) %>% map(prop.table)

# Cluster analysis ----

?agnes
datasetn <- dataset %>% 
  select_if(is.numeric)
agnclus <- agnes(datasetn, metric = "euclidean", 
                 stand = FALSE, method="complete")
# stand = FALSE porque no queremos estandarizar
summary(agnclus)
plot(agnclus, main=paste("Agnes:",agnclus$method,sep=""))
agnclus$ac  
