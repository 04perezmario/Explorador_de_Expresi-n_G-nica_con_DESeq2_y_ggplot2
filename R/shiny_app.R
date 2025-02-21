# shiny_app.R

library(shiny)
library(DESeq2)
library(ggplot2)

ui <- fluidPage(
  titlePanel("Análisis de Expresión de Genes"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("readcounts_file", "Sube el archivo de Read Counts:", accept = c(".csv", ".txt")),
      fileInput("samplesheet_file", "Sube el archivo de Sample Sheet:", accept = c(".csv", ".txt")),
      fileInput("deseqres_file", "Sube el archivo de DESeq Results:", accept = c(".csv", ".txt")),
      
      textInput("gene", "Introduce el nombre del gen:", value = "MYC"),
      actionButton("update", "Actualizar Gráfico")
    ),
    
    mainPanel(
      plotOutput("genePlot"),
      verbatimTextOutput("debugOutput")  # Para mostrar la depuración
    )
  )
)

server <- function(input, output, session) {
  
  # Función para detectar la columna que contiene los símbolos de genes
  detectar_columna_genes <- function(df) {
    for (col in colnames(df)) {
      valores <- as.character(df[[col]])
      
      # Verifica si al menos un valor tiene letras y números (ejemplo: MYC1, ASIC1)
      if (any(grepl("[A-Za-z]", valores) & grepl("[0-9]", valores), na.rm = TRUE)) {
        return(col)  # Retorna el nombre de la columna si cumple con la condición
      }
    }
    return(NULL)  # Si no encuentra, devuelve NULL
  }
  
  # Cargar archivos una vez el usuario los suba
  readcounts <- reactive({
    req(input$readcounts_file)
    
    # Leer el archivo TXT asegurando que los separadores sean correctos
    data <- read.delim(input$readcounts_file$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Mostrar los primeros 5 datos cargados para depuración
    output$debugOutput <- renderPrint({
      # Mostrar los primeros registros de la columna de genes y sus nombres
      gene_column <- detectar_columna_genes(data)
      if (!is.null(gene_column)) {
        head(data[[gene_column]])
      } else {
        "No se encontró una columna de genes"
      }
    })
    
    # Detectar la columna que contiene los nombres de genes
    gene_column <- detectar_columna_genes(data)
    
    if (is.null(gene_column)) {
      stop("No se encontró una columna con identificadores de genes (nombres con letras y números).")
    }
    
    # Establecer la columna de genes como nombres de fila
    rownames(data) <- trimws(data[[gene_column]])  # Trim para evitar problemas de espacios
    data <- data[, !colnames(data) %in% gene_column]  # Eliminar la columna de genes del dataframe
    
    return(data)
  })
  
  samplesheet <- reactive({
    req(input$samplesheet_file)
    read.csv(input$samplesheet_file$datapath, header = TRUE, row.names = 1)
  })
  
  deseqres <- reactive({
    req(input$deseqres_file)
    read.csv(input$deseqres_file$datapath, header = TRUE, row.names = 1)
  })
  
  # Renderizar el gráfico del gen seleccionado
  output$genePlot <- renderPlot({
    req(input$gene)
    data <- readcounts()
    
    # Convertir el gen ingresado a mayúsculas para compararlo de manera insensible al caso
    gene_input <- toupper(input$gene)
    
    # Depuración: Verificar si el gen existe en los datos
    if (gene_input %in% toupper(rownames(data))) {
      # Detectar las clases dinámicamente basadas en los nombres de las columnas
      class_labels <- colnames(data)
      
      # Organizar los datos para hacer un boxplot por clase
      gene_expression <- as.numeric(data[gene_input,])
      
      # Crear un dataframe para usar en ggplot
      gene_df <- data.frame(
        Expression = gene_expression,
        Class = factor(class_labels)
      )
      
      # Usar ggplot para crear el gráfico de caja por clase
      ggplot(gene_df, aes(x = Class, y = Expression, fill = Class)) +
        geom_boxplot() +
        labs(title = paste0(input$gene, " - Expresión por Clase"),
             y = "Expresión del Gen",
             x = "Clase") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
    } else {
      plot.new()
      text(0.5, 0.5, paste("El gen", input$gene, "no se encuentra en los datos."), cex = 1.5)
    }
  })
}

shinyApp(ui = ui, server = server)

