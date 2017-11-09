library(shiny)
library(tidyverse)
library(cowplot)
library(ggsci)
library(DT)
library(stringr)
library(ggthemes)
library(scales)

df = read.csv('short.hybrid.cluster.data.csv', header=TRUE,
              colClasses=c('character', 'numeric', 'factor', 'numeric', 'character', 'factor',
                           'numeric', 'character', 'character','numeric','numeric','numeric','numeric','numeric'))
df$virsorter_category <- as.factor(df$virsorter_category)

cluster_df <- read.csv('clusters.csv', header=TRUE,
                       colClasses=c('character', 'numeric', 'numeric', 'character', 'character'))

source('functionality.R')

ui <- navbarPage("Comparison of Hybrid Assemblies",
                 tabPanel("Contigs",
                          fluidRow(
                            
                            column(2,
                   sliderInput(inputId = "min_length",
                     label="Minimum contig length:",
                      min = 2500,
                      max = 500000,
                      step=1000,
                      value = 10000),
                  checkboxInput(inputId = 'VirFinder',
                    label='VirFinder viral'),
  
                  checkboxInput(inputId = 'VirSorter',
                    label='VirSorter viral'),
  
                  checkboxInput(inputId = 'VirSorter_circular',
                    label='VirSorter circular'),
  
                  selectInput(inputId = 'read_type',
                    label = 'Contig Type',
                    choices=c('Short-only' = 'short', 'Hybrid' = 'hybrid', 'Both' = 'both'),
                    selected = 'both')),
                  
                          column(10,
  
                  plotOutput(outputId = 'contig_plot',
                    dblclick = "contig_plot_dblclick",
                    click = "contig_plot_click",
                    brush = brushOpts(id = "contig_plot_brush", resetOnNew = TRUE))),
                  fluidRow(
                    column(12,
                  tableOutput(outputId = "contig_info"))))),
  
              tabPanel("Clusters",
  
                selectInput(inputId = 'hybrid_rep',
                  label='Cluster Type:',
                  choices=unique(cluster_df$type)),
  
                DT::dataTableOutput(outputId = 'clusters'),
  
                plotOutput(outputId = 'cluster_summary'))
  
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