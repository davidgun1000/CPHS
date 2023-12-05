
library(coda)
library(readxl)
library(readr)

STANDARDGARCH_N1000 <- read_csv("STANDARDGARCH1000.csv", col_names = FALSE)

A=effectiveSize(STANDARDGARCH_N1000)
IACT_A=10000/A
write.csv(IACT_A,"IACT_STANDARDGARCH_N1000.csv")
