

library(metaDigitise)
options(scipen = 100)# setting the options to not show scientific notation

dir <- paste0(getwd(),"/paper_digitization") 


data <- metaDigitise(dir = dir)


View(data)



# calculating some effects sizes from Chi-square data
library(esc)

esc_chisq(22.67, totaln = 600, es.type = "d")
esc_f()
