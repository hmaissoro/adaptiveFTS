library(data.table)
library(matrixStats)

all_func <- list.files("./R/")
lapply(all_func, function(func){
if (func != "zzz.R") source(file = file.path("./R", func), echo = FALSE)
})

