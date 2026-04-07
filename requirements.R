# MetaDrugSig — R package dependencies
# Run this script before launching the portal

install.packages(c(
  "shiny",
  "data.table",
  "dplyr",
  "ggplot2",
  "plotly",
  "metafor",
  "patchwork",
  "DT"
))

# msigdbr from CRAN
install.packages("msigdbr")

# Verify installations
required <- c("shiny","data.table","dplyr","ggplot2",
               "plotly","metafor","patchwork","msigdbr","DT")
missing  <- required[!required %in% installed.packages()[,"Package"]]

if (length(missing) == 0) {
  cat("All packages installed successfully.\n")
} else {
  cat("Missing packages:", paste(missing, collapse=", "), "\n")
}
