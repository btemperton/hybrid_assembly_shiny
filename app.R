library(shiny)
library(tidyverse)
library(cowplot)
library(DT)
library(stringr)
library(ggthemes)
library(ggsci)
library(scales)

df = read.csv('short.hybrid.cluster.data.csv', header=TRUE,
              colClasses=c('character', 'numeric', 'factor', 'numeric', 'character', 'factor',
                           'numeric', 'character', 'character','numeric','numeric','numeric','numeric','numeric'))
df$virsorter_category <- as.factor(df$virsorter_category)

cluster_df <- read.csv('clusters.csv', header=TRUE,
                       colClasses=c('character', 'numeric', 'numeric', 'character', 'character'))

source('functionality.R')
source('contigsPanel.R')
source('comparisonPanel.R')

ui <- navbarPage("Comparison of Hybrid Assemblies",
                 load.contigs.panel(),
                 load.comparison.panel(cluster_df)
  
)
server <- function(input, output) {
  

  output$contigs <- renderDataTable(prepare.contig.dataframe(input, df))
  
  output$clusters <- DT::renderDataTable(prepare.cluster.dataframe(input, cluster_df),
                                         selection='single')
  
  output$cluster_summary <- renderPlot(plot.cluster.summary(input, df, cluster_df))
  
  ranges <- reactiveValues(x = NULL, y = NULL)
  output$contig_plot <- renderPlot(plot.GC.vs.coverage(input, df, ranges))
  
  output$contig_info <- renderTable({
    
    nearPoints(apply.filters(input, df),
               input$contig_plot_click, threshold = 10, maxpoints = 15,
               addDist = TRUE)
  })
  
  output$summary_contig_plot <- renderPlot(plot.contig.summary(input, df))
  
  # Allows for zooming on the plot
  observeEvent(input$contig_plot_dblclick, {
    brush <- input$contig_plot_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })

}
shinyApp(ui = ui, server = server)