library(shiny)       # Carga la librería Shiny, que se usa para crear aplicaciones web interactivas en R.
library(ggplot2)     # Carga la librería ggplot2, que permite crear gráficos avanzados en R.
library(devtools)    # Carga la librería devtools, usada para instalar paquetes desde GitHub.

# Intentamos instalar EnhancedVolcano desde GitHub si no está disponible
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {  # Verifica si el paquete EnhancedVolcano está instalado.
  install_github("kevinblighe/EnhancedVolcano")              # Si no está instalado, lo descarga e instala desde GitHub.
}
library(EnhancedVolcano)  # Carga la librería EnhancedVolcano, usada para generar Volcano Plots.

# Definición de la interfaz de usuario (UI) con fluidPage
ui <- fluidPage(
  titlePanel("Análisis de Expresión de Genes"),  # Define el título principal de la aplicación.
  
  sidebarLayout(  # Crea un diseño con un panel lateral y un panel principal.
    sidebarPanel(  # Panel lateral con controles de entrada.
      fileInput("readcounts_file", "Sube el archivo de Read Counts:", accept = c(".csv", ".txt")),  
      # Permite subir un archivo con los datos de expresión génica.
      
      fileInput("samplesheet_file", "Sube el archivo de Sample Sheet:", accept = c(".csv", ".txt")),  
      # Permite subir un archivo con información de las muestras.
      
      fileInput("deseqres_file", "Sube el archivo de DESeq Results (opcional):", accept = c(".csv", ".txt")),  
      # Permite subir un archivo con los resultados de análisis diferencial de expresión.
      
      textInput("gene", "Introduce el nombre del gen:", value = "MYC"),  
      # Campo de texto donde el usuario introduce un gen de interés (por defecto "MYC").
      
      selectInput("plot_type", "Selecciona el tipo de gráfico:", 
                  choices = c("Expresión por Grupo" = "group", 
                              "Expresión por Clase" = "class")),  
      # Menú desplegable para seleccionar el tipo de gráfico a visualizar.
      
      #actionButton("update", "Actualizar Gráfico"),  
      # Botón que permite actualizar el gráfico según la selección del usuario.
      
      hr(),  # Línea horizontal para separar secciones.
      
      checkboxInput("show_volcano", "Mostrar Volcano Plot", FALSE)  
      # Casilla de verificación para activar o desactivar la visualización del Volcano Plot.
    ),
    
    #
    mainPanel(  # Panel principal donde se muestran los gráficos.
      plotOutput("genePlot"),  # Área donde se mostrará el gráfico de expresión del gen seleccionado.
      
      plotOutput("volcanoPlot", height = "500px"),  
      # Área donde se mostrará el Volcano Plot (con una altura de 500 píxeles).
      
      verbatimTextOutput("debugOutput")  
      # Área para mostrar mensajes de depuración y ayudar en la solución de errores.
    )
  )
)

# Definición del servidor
server <- function(input, output, session) {
  
  #=============================#
  # Cargar archivo de Read Counts
  #=============================#
  readcounts <- reactive({
    
    # Asegura que el archivo se haya subido antes de ejecutar el código.
    req(input$readcounts_file)  
    
    # Lee el archivo de Read Counts en formato tabular (TSV o CSV).
    data <- read.delim(input$readcounts_file$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)  
    
    #=============================#
    # Procesamiento de la tabla
    #=============================#
    
    # Verifica si la columna "Gene_Symbol" está presente en los datos.
    if ("Gene_Symbol" %in% colnames(data)) {  
      rownames(data) <- data$Gene_Symbol  # Usa "Gene_Symbol" como nombres de fila.
    } else {
      rownames(data) <- data[, 1]  # Si no hay "Gene_Symbol", usa la primera columna como nombres de fila.
    }
    
    # Se eliminan las dos primeras columnas:
    # - "Gene_ID" porque es solo un identificador sin valores de expresión.
    # - "Gene_Symbol" porque ya se ha usado como nombres de fila.
    data <- data[, -c(1, 2)]  
    
    # Renombra las columnas eliminando el prefijo "N_" si estuviera presente.
    colnames(data) <- gsub("^N_", "", colnames(data))  
    
    return(data)  # Devuelve la tabla procesada.
  })
  
  #=============================#
  # Cargar archivo Sample Sheet
  #=============================#
  samplesheet <- reactive({
    # Asegura que el archivo se haya subido antes de ejecutar el código.
    req(input$samplesheet_file)  
    
    # Lee el archivo Sample Sheet en formato tabular.
    sheet_data <- read.delim(input$samplesheet_file$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)  
    
    # Renombra la columna "Sample.Group" a "Group" para mayor claridad.
    colnames(sheet_data)[colnames(sheet_data) == "Sample.Group"] <- "Group"  
    
    return(sheet_data)  # Devuelve la tabla procesada.
  })
  
  #=============================#
  # Cargar resultados de DESeq2
  #=============================#
  deseqres <- reactive({
    # Asegura que el archivo se haya subido antes de ejecutar el código.
    req(input$deseqres_file)  
    
    # Lee el archivo con los resultados de DESeq2.
    data <- read.delim(input$deseqres_file$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)  
    
    # Elimina filas con valores NA en cualquiera de sus columnas.
    data <- na.omit(data)  
    
    #=============================#
    # Clasificación de genes significativos
    #=============================#
    
    # Crea una nueva columna "Significativo":
    # - "Sí" si el gen tiene un p-adj < 0.05 y un cambio en expresión mayor a 1 en valor absoluto.
    # - "No" en caso contrario.
    data$Significativo <- ifelse(data$padj < 0.05 & abs(data$log2FoldChange) > 1, "Sí", "No")  
    
    return(data)  # Devuelve la tabla procesada.
  })
  




  # ================================================================================================================= 
  
                        # Generar gráfico de expresión del gen por grupo
  
  # ================================================================================================================= 
  output$genePlot <- renderPlot({ # Inicia la renderización del gráfico de expresión en la interfaz de usuario.
   
     req(input$gene, readcounts(), samplesheet())  # Asegura que se proporcionen los datos necesarios: el gen seleccionado, los datos de lectura y la hoja de muestras.
    
    data <- readcounts()  # Asigna los datos de lectura a la variable 'data'.
    samplesheet_data <- samplesheet()  # Asigna los datos de la hoja de muestras a la variable 'samplesheet_data'.
    gene_input <- trimws(toupper(input$gene))  # Toma el gen ingresado, lo convierte a mayúsculas y elimina espacios en blanco.
    
    available_genes <- toupper(rownames(data))  # Convierte los nombres de los genes en los datos de lectura a mayúsculas para compararlos.
    if (!(gene_input %in% available_genes)) {  # Verifica si el gen ingresado está en los datos disponibles.
      plot.new()  # Crea una nueva trama vacía si el gen no está en los datos.
      text(0.5, 0.5, paste("El gen", input$gene, "no se encuentra en los datos."), cex = 1.5)  # Muestra un mensaje indicando que el gen no se encuentra en los datos.
      return()  # Sale de la función si el gen no se encuentra.
    }
    
    
    gene_expression <- as.numeric(data[which(available_genes == gene_input),])  # Extrae la expresión del gen seleccionado de los datos de lectura.
    
    gene_df <- data.frame(  # Crea un dataframe con la expresión del gen y los nombres de las muestras.
      Sample = colnames(data),
      Expression = gene_expression
    )
    # Gráfico de expresión por grupo
    if (input$plot_type == "group") {  # Verifica si se seleccionó el tipo de gráfico "group" para visualizar la expresión por grupo.
      gene_df <- merge(gene_df, samplesheet_data, by.x = "Sample", by.y = "Sample.ID", all.x = TRUE)  # Combina los datos de expresión con la hoja de muestras por la columna 'Sample'.
      
      if (any(is.na(gene_df$Group))) {  # Verifica si alguna muestra no tiene un grupo asignado.
        stop("Algunas muestras no tienen un grupo asignado. Revisa los nombres en Sample Sheet y Read Counts.")  # Detiene la ejecución si hay muestras sin grupo asignado.
      }
      ggplot(gene_df, aes(x = Group, y = Expression, fill = Group)) +  # Crea un gráfico de caja utilizando ggplot, con 'Group' en el eje x y la expresión del gen en el eje y.
        geom_boxplot() +  # Agrega los elementos del gráfico de caja.
        labs(title = paste(input$gene, " - Expresión por Grupo (CPM)"),  # Define el título del gráfico, que incluye el nombre del gen.
             y = "Expresión del Gen (CPM)",  # Etiqueta del eje y.
             x = "Grupo") +  # Etiqueta del eje x.
        theme_minimal() +  # Aplica un tema minimalista al gráfico.
        theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Ajusta la rotación de las etiquetas del eje x para mayor legibilidad.
    } 
    
    
    # ================================================================================================================= 
    
                             # Gráfico de expresión por clase
    
    # ================================================================================================================= 
    else {  # Si el tipo de gráfico no es "group", se ejecuta este bloque (probablemente cuando el tipo de gráfico es "class").
      class_labels <- colnames(data)  # Asigna los nombres de las columnas de los datos a la variable 'class_labels', que representan las clases de las muestras.
      gene_df <- data.frame(  # Crea un dataframe con la expresión del gen y las clases asociadas.
        Expression = gene_expression,  # Columna de expresión del gen.
        Class = factor(class_labels)  # Columna de clase, convertida en factor (categoría) para su uso en el gráfico.
      )
      
      ggplot(gene_df, aes(x = Class, y = Expression, fill = Class)) +  # Crea un gráfico de caja con las clases en el eje x y la expresión en el eje y.
        geom_boxplot() +  # Agrega el gráfico de caja al gráfico.
        labs(title = paste(input$gene, " - Expresión por Clase (CPM)"),  # Define el título del gráfico, que incluye el nombre del gen.
             y = "Expresión del Gen (CPM)",  # Etiqueta del eje y.
             x = "Clase") +  # Etiqueta del eje x.
        theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Ajusta la rotación de las etiquetas del eje x para mayor legibilidad.
    }
    
  })
  
  # ================================================================================================================= 
  
            # Generar Volcano Plot con EnhancedVolcano
  
  # ================================================================================================================= 

  output$volcanoPlot <- renderPlot({
    # Aseguramos que se haya activado el botón 'show_volcano' y que los resultados de deseqres() estén disponibles
    req(input$show_volcano, deseqres())
    
    # Asignamos los resultados de deseqres() a la variable 'volcano_data'
    volcano_data <- deseqres()
    
    # Filtramos los valores NA para asegurarnos de no tener datos incompletos (log2FoldChange, pvalue, padj)
    volcano_data <- na.omit(volcano_data[, c("symbol", "log2FoldChange", "pvalue", "padj")])
    
    # Generamos el Volcano Plot usando la librería EnhancedVolcano
    EnhancedVolcano(volcano_data,
                    lab = volcano_data$symbol,  # Etiquetas de los genes usando la columna 'symbol'
                    x = 'log2FoldChange',        # Eje X basado en 'log2FoldChange'
                    y = 'pvalue',                # Eje Y basado en 'pvalue'
                    pCutoff = 0.05,              # Umbral de significancia para el p-value (0.05)
                    FCcutoff = 1,                # Umbral para el log2FoldChange (cambio de expresión mayor que 1 o menor que -1)
                    pointSize = 3.0,             # Tamaño de los puntos en el gráfico
                    colAlpha = 0.5,              # Transparencia de los puntos (0.5)
                    legendPosition = 'top',      # Posición de la leyenda en la parte superior
                    legendLabSize = 12,          # Tamaño de la fuente de la leyenda
                    caption = 'Volcano Plot con EnhancedVolcano',  # Título del gráfico
                    drawConnectors = TRUE,       # Activar líneas de conexión entre los puntos significativos
                    widthConnectors = 0.5        # Grosor de las líneas de conexión
    ) 
  })
  
  
}
shinyApp(ui = ui, server = server)


