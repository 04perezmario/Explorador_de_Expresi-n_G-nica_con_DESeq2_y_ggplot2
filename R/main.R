# main.R

# Establecer el directorio de trabajo (si es necesario)
setwd("/home/labs/lslab/studentlslab/Desktop/ShinyAPPMario")

# Cargar los scripts necesarios
source("R/load_data.R")
source("R/deseq_analysis.R")
source("R/visualization.R")
source("R/shiny_app.R")

# Cargar el paquete Shiny
library(shiny)

# Ejecutar la aplicación Shiny
shinyApp(ui = ui, server = server)
