
library(coda)
library(readxl)
library(readr)

ParPMMH_data1 <- read_csv("TDIST_HMC.csv", col_names = FALSE)

A=effectiveSize(ParPMMH_data1)
IACT_A=length(ParPMMH_data1$X1)/A
write.csv(IACT_A,"IACT_TDIST_HMC.csv")
write.csv(A,"ESS_TDIST_HMC.csv")






