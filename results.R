#-------------------------------------------------------------------------------
# analysis for chapter 1: temporal population genetics


# libraries --------------------------------------------------------------------
library(tidyverse)
library(dartR)
library(emR)

# load data --------------------------------------------------------------------
peaks <- emR::em.data.peaks
fst <- emR::em.data.ph_fst_grids[[4]] # min sample size 4
fis <- emR::em.data.ph_fis_grids[[4]] 
glimpse(peaks)
head(peaks)
