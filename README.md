# Outlier detection methods for SIR models (**outlierSIR**)

**outlierSIR** is a sampling method allowing to detect outlier or extreme individuals.

## Installation

There is currently only one way to install **outlierSIR**

  * From the under development repository from GitHub thanks to `devtools`
  
  ```r
  # install.packages("outlierSIR")
  devtools::install_github("hlorenzo/outlierSIR")
  ```

Once that package is installed, you can test it on the ozone dataset

  ```r
  data(ozone)
  # MONO method
  output <- outlierSIR(ozone[,-c(1,12,13)],ozone[,1],method="MONO")
  # TTR method
  # output <- outlierSIR(ozone[,-c(1,12,13)],ozone[,1],method="TTR",
  # plotBandwidthEst = TRUE,kernel="normal",nb.replications = 100)
  # BOOT method
  # output <- outlierSIR(ozone[,-c(1,12,13)],ozone[,1],method="BOOT",
  # plotBandwidthEst = TRUE,kernel="normal",nb.replications = 100)
  ```
