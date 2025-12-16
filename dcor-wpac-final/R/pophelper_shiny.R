install.packages(c("ggplot2","gridExtra","label.switching","tidyr","remotes",
                   "colourpicker","DT","highcharter","htmlwidgets","magrittr",
                   "markdown","RColorBrewer","shiny","shinyAce","shinyBS",
                   "shinythemes","shinyWidgets","viridisLite","writexl"),
                 repos = "http://cran.us.r-project.org")

# install pophelper package from GitHub
remotes::install_github('royfrancis/pophelper', force=TRUE, dependencies = TRUE)
library(pophelper)

# install the package from GitHub
remotes::install_github('royfrancis/pophelperShiny', force=TRUE)
library(pophelperShiny)
library(shiny)
# launch app
runPophelper()
