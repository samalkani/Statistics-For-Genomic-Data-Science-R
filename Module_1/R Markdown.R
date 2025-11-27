# R Markdown

## Figures
x = rnorm(100)
plot(x, col=3, pch=19)


## Figures no code
x = rnorm(100)
plot(x, col=3, pch=19)

##  Using cache=FALSE the first time it will run slow (delay for 10 seconds)
Sys.sleep(10)

### Using cache=TRUE the second time it will run fast (recompile)
Sys.sleep(10)

# Session Info & Date
sessionInfo()
BiocManager::install("devtools")
library(devtools)
devtools::session_info()
Sys.Date()

# Time and date: `r Sys.Date()`
