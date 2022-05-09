
library(shiny)
library(bslib)
library(ggplot2)
library(colourpicker)


deseq_choices <-
  c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
# Define UI for application that draws a histogram
ui <- fluidPage(
  theme = bs_theme(bootswatch = "minty"),
  # Application title

  markdown(paste0("To use this application, download the CSV `deseq_res.csv`",
                  " from the data directory of this app's repository.")),
  
  # Sidebar with a slider input for p-value magnitude
  sidebarLayout(
    sidebarPanel(
      fileInput(
        "file",
        label = "Load differential expression results",
        accept = c("text/csv",
                   "text/comma-separated-values,text/plain",
                   ".csv"),
        placeholder = "deseq_res.csv"
      ),
      HTML(
        paste(
          "A volcano plot can be generated with",
          "<b>\"log<sub>2</sub> fold-change\"</b> on the x-axis and",
          "<b>\"p-adjusted\"</b> on the y-axis.<br>"
        )
      ),
      br(),
      radioButtons(
        inputId = "x_axis",
        label = "Choose the column for the x-axis",
        choices = deseq_choices,
        selected = "log2FoldChange"
      ),
      radioButtons(
        inputId = "y_axis",
        label = "Choose the column for the y-axis",
        choices = deseq_choices,
        selected = "padj"
      ),
      colourInput(
        inputId = "base",
        label = "Base point color",
        value = "#22577A",
        closeOnClick = T
      ),
      colourInput(
        inputId = "highlight",
        label = "Highlight point color",
        value = "#FFCF56",
        closeOnClick = T
      ),
      sliderInput(
        "slider",
        "Select the magnitude of the p adjusted coloring:",
        min = -300,
        max = 0,
        value = -150,
        step = 1
      ),
      submitButton(
        text = "Plot",
        icon = icon("car-crash"),
        width = "100%"
      )
    ),
    # Show the volcano plot
    mainPanel(tabsetPanel(
      tabPanel("Plot", {
        plotOutput("volcano")
      }),
      tabPanel("Table",
               tableOutput("table"))
    ))
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  

  load_data <- reactive({
    df <- read.csv(input$file$datapath)
    colnames(df)[1] <- "gene"
    return(df)
  })
  

  volcano_plot <-
    function(dataf, x_name, y_name, slider, color1, color2) {
      p <- ggplot(dataf, aes(x = !!sym(x_name),
                             y = -log10(!!sym(y_name)))) +
        geom_point(aes(color = !!sym(y_name) < 1 * 10 ^ (as.numeric(slider)))) +
        theme_bw() +
        scale_color_manual(values = c(color1, color2)) +
        theme(legend.position = "bottom") +
        labs(color = paste0(y_name, " < 1 Ã 10^", slider))
      return(p)
    }
  
  #' Draw and filter table

  draw_table <- function(dataf, slider) {
    df_out <- dataf[which(dataf$padj < 1 * 10 ^ (as.numeric(slider))),]
    df_out$pvalue <- formatC(df_out$pvalue, digits = -2)
    df_out$padj <- formatC(df_out$padj, digits = -2)
    return(df_out)
  }
  
  #' These outputs aren't really functions, so they don't get a full skeleton, 
  #' but use the renderPlot() and renderTabel() functions to return() a plot 
  #' or table object, and those will be displayed in your application.
  output$volcano <- renderPlot({
    req(input$file)
    df <- load_data()
    p <-volcano_plot(df,
                     input$x_axis,
                     input$y_axis,
                     input$slider,
                     input$base,
                     input$highlight)
    return(p)
  }, height = 700)
  
  # Same here, just return the table as you want to see it in the web page
  output$table <- renderTable({
    req(input$file)
    table <- load_data()
    colnames(table)[1] <- "gene"
    return(draw_table(dataf = table, slider = input$slider))
  }, striped = T)
}

# Run the application
shinyApp(ui = ui, server = server)
