#Tab2
library(tidyverse)
library(shiny)
library(ggplot2)
library(DT)

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      fileInput(inputId = 'normcount', label = 'Load normalized counts  CSV'),
      sliderInput(inputId = 'percvar', label = 'Select minimum percentile variance of genes', min = 0, max = 100, value = 80),
      sliderInput(inputId = 'nonzero', label = 'Select  minimum number of non-zero samples', min = 0, max = 100, value = 10),
      submitButton(text='Submit')
    ),
    mainPanel(
      tabsetPanel(
        tabPanel('Summary', tableOutput(outputId = 'normalised')),
        tabPanel('Diagnostic Plots', plotOutput(outputId = 'medvariance'), plotOutput(outputId = 'medzero')),
        tabPanel('Heatmap', plotOutput(outputId = "heat_map")),
        tabPanel('PCA', selectInput(inputId = "comp1", label="Select X-axis", choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")),
                 selectInput(inputId = "comp2", label="Select Y-axis", choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"), selected = "PC2"),
                 plotOutput(outputId = "PCAplot"))
      ))
  ))
server <- function(input, output, session){
  options(shiny.maxRequestSize=30*1024^2)

  load_counts <- reactive({
    if (!is.null(input$normcount)){
      meta <- read_csv(input$normcount$datapath)
      return(meta)}
    else{
      return(NULL)
    }
  })

  summary_count <- function(count, pv, ngenes){
    if (!is.null(input$normcount)){
      tot_samp <- ncol(count)-1  
      tot_genes <- nrow(count)
      count <- count %>% mutate(variance = apply(count[-1], MARGIN = 1, FUN = var))  
      perc_val <- quantile(count$variance, probs = pv/100)   
      count <- filter(count, variance >= perc_val)  
      count <- na_if(count, 0)    
      count$non_zero <- tot_samp-rowSums(is.na(count))  
      count <- filter(count, non_zero >= ngenes)  
      filt_genes <- nrow(count)  
      perc_pass_genes <- filt_genes/tot_genes*100
      fail_genes <- tot_genes-filt_genes
      perc_fail_genes <- fail_genes/tot_genes*100
      
      summ_tib <- tibble('Measure' = c('Number of Samples', 'Number of Genes', 'Number of genes passing filters', " (%) passing filter", 'No. of genes not passing filters', ' (%) not passing filter'),
                         'Value' = c(tot_samp, tot_genes, filt_genes, perc_pass_genes, fail_genes, perc_fail_genes))
      return(summ_tib)}
    else{return(NULL)}
  }
  
  med_vs_var <- function(count, pv){
    if (!is.null(input$normcount)){
    
      plot_tib <- count%>%
        mutate(Median = apply(count[-1], MARGIN = 1, FUN = median), 
               Variance = apply(count[-1], MARGIN = 1, FUN = var))
      perc_val <- quantile(plot_tib$Variance, probs = pv/100)   
      plot_tib <- plot_tib %>% mutate(thresh = case_when(Variance >= perc_val ~ "TRUE", TRUE ~ "FALSE"))
      
      cols <- c("FALSE" = "red", "TRUE" = "black")
      scatter <- ggplot(plot_tib, aes(Median, Variance))+
        geom_point(aes(color=thresh), alpha=0.75)+
        scale_color_manual(values = cols)+
        labs(title = 'Plot of Median vs Variance.', subtitle = "Genes filtered out are in red.")+
        scale_y_log10()+
        scale_x_log10()+
        theme_bw()
      
      return(scatter)}
    else{return(NULL)}
  }
  
  
  # median vs non-zero samples plot
  med_vs_nz <- function(count, ngenes){
    if (!is.null(input$normcount)){
      tot_samp <- ncol(count)-1  
      
      plot_tib <- count %>%   
        mutate(Median = apply(count[-1], MARGIN = 1, FUN = median)) %>% na_if(0)  
      plot_tib$no_zeros <- rowSums(is.na(plot_tib))  
      plot_tib <- plot_tib %>% mutate(thresh = case_when(no_zeros <= ngenes ~ "TRUE", TRUE ~ "FALSE"))
      
      
      #plot scatter plot
      cols <- c("FALSE" = "red", "TRUE" = "black")
      scatter <- ggplot(plot_tib, aes(Median, no_zeros))+
        geom_point(aes(color=thresh), alpha=0.75)+
        scale_color_manual(values = cols)+
        scale_x_log10()+
        labs(title = 'Plot of Median vs Number of Non-Zero genes', subtitle = "Genes filtered out are in red.")+
        ylab('Number of samples with zero count')
      return(scatter)}
    else{return(NULL)}
  }
  
  # PCA plot
  plot_pca <- function(count, pv, comp1, comp2){
    if (!is.null(input$normcount)){
    
      filt_tib <- count %>% 
        mutate(variance = apply(count[-1], MARGIN = 1, FUN = var), .after = gene)
      perc_val <- quantile(filt_tib$variance, probs = pv/100, na.rm = TRUE)   
      filt_tib <- filter(filt_tib, variance >= perc_val) 
      pca_res <- prcomp(t(filt_tib[-c(1,2)]), scale = FALSE) 
    
      variance <- summary(pca_res)$importance[2,]
      x <- round(variance[comp1]*100, 2)
      y <- round(variance[comp2]*100, 2)
     
      
       # PCA plot
      plot_tib <- tibble(PC1 = pca_res$x[,comp1], PC2=pca_res$x[,comp2])
      pca <- ggplot(plot_tib, aes(PC1, PC2))+
        geom_point()+
        labs(title="Princple Component Analysis Plot")+
        xlab(str_c(comp1, x, "% variance", sep=" "))+
        ylab(str_c(comp2, y, "% variance", sep=" "))+
        theme_bw()
      return(pca)}
    else{return(NULL)}
  }
  
  #HEATMAPS
  plot_heatmap <- function(count, pv, ngenes){
    if (!is.null(input$normcount)){
      count <- na_if(count, 0)
      count$no_zeros <- rowSums(is.na(count))  
      count <- filter(count, no_zeros <= ngenes)
      count <- log10(count[,!colnames(count) %in% c("gene", "no_zeros")])  
  
      plot_tib <- count %>% 
        mutate(variance = apply(count, MARGIN = 1, FUN = var)) 
      perc_val <- quantile(plot_tib$variance, probs = pv/100, na.rm = TRUE)   
      plot_tib <- filter(plot_tib, variance >= perc_val) 
      heat_map <- heatmap(as.matrix(plot_tib[-ncol(plot_tib)]), scale = "row")
      return(heat_map)}
    else{return(NULL)}
  }
  
  output$normalised <- renderTable({
    summary_count(load_counts(), input$percvar, input$nonzero)
  })
  output$medvariance <- renderPlot({
    med_vs_var(load_counts(), input$percvar)
  })
  output$medzero <- renderPlot({
    med_vs_nz(load_counts(), input$nonzero)
  })
  output$heat_map <- renderPlot({
    plot_heatmap(load_counts(), input$percvar, input$nonzero)
  })
  output$PCAplot <- renderPlot({
    plot_pca(load_counts(), input$percvar, input$comp1, input$comp2)
  })
}

shinyApp(ui=ui, server = server)