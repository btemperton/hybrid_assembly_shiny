library(shiny)
library(tidyverse)
library(cowplot)
library(DT)
library(stringr)
library(ggthemes)
library(scales)

df = read.csv('short.hybrid.cluster.data.csv', header=TRUE,
              colClasses=c('character', 'numeric', 'factor', 'numeric', 'character', 'factor',
                           'numeric', 'character', 'character','numeric','numeric','numeric','numeric','numeric'))
df$virsorter_category <- as.factor(df$virsorter_category)

cluster_df <- read.csv('clusters.csv', header=TRUE,
                       colClasses=c('character', 'numeric', 'numeric', 'character', 'character', 'numeric', 'numeric', 'numeric'))

gene_map <- read.csv('gene.map.representatives.csv', header=TRUE,
                     colClasses=c('character', 'character', 'numeric', 'numeric', 'character', 'character'))

vcontact_map <-read.csv('vcontact.membership.csv', header=TRUE, 
                        colClasses = c('character', 'character', 'character', 'numeric'))


source('functionality.R')
source('contigsPanel.R')
source('comparisonPanel.R')

ui <- navbarPage("Comparison of Hybrid Assemblies",
                 load.contigs.panel(),
                 load.comparison.panel(cluster_df)
  
)
server <- function(input, output) {
  

  output$contigs <- DT::renderDataTable(prepare.contig.dataframe(input, df))
  
  output$clusters <- DT::renderDataTable(
    datatable(prepare.cluster.dataframe(input, cluster_df) %>% select(-pos_cluster), 
              rownames=FALSE,
              class='compact',
              colnames=c('Cluster Name', 'Length', '# Members', 'Name of Representatitve', 'Type', 'Number of viral members', 'Viral ratio'),
              selection='single',
              options = list(pageLength = 10,
                             initComplete = JS(
                               "function(settings, json) {",
                               "$(this.api().table().header()).css({'background-color': '#268bd2', 'color': '#fff'});",
                               "}"),
                             columnDefs = list(list(className = 'dt-right', targets=c(3,4)),
                                               list(width = '60px', targets=c(0,1,2,5))))) %>% DT::formatRound(c(7)))
  
  output$cluster_counts <- renderPlot(plot.cluster.counts(input, cluster_df))
  
  output$cluster_summary <- renderPlot(plot.cluster.summary(input, df, cluster_df))
  

  output$cluster_rep_genes <-DT::renderDataTable(
    datatable(prepare.protein.dataframe(input, cluster_df, gene_map),
              rownames=FALSE,
              class='compact',
              colnames=c('Contig', 'start', 'end', 'strand', 'Best NR Hit', 'VOG'),
              selection = 'single',
              options = list(pageLength = 10,
                             columnDefs = list(list(className = 'dt-right', targets=c(1,2,3)),
                                               list(width = '100px', targets=c(0)),
                                               list(width = '60px', targets=c(1,2,3))),
                             initComplete = JS(
                               "function(settings, json) {",
                               "$(this.api().table().header()).css({'background-color': '#268bd2', 'color': '#fff'});",
                               "}"))))
  
  
  output$cluster_members <- DT::renderDataTable(
    datatable(prepare.cluster.members.dataframe(input, df, cluster_df),
              rownames = FALSE,
              class='compact',
              colnames=c('Contig Name', 'Length (bp)'),
              selection = 'none',
              options = list(pageLength = 10,
                             initComplete = JS(
                               "function(settings, json) {",
                               "$(this.api().table().header()).css({'background-color': '#756bb1', 'color': '#fff'});",
                               "}"))))
    
  output$vcontact_members <- DT::renderDataTable(
    datatable(prepare.vcontact.dataframe(input, cluster_df, vcontact_map),
                rownames=FALSE,
                class='compact',
                colnames=c('id', 'Name', 'Taxonomy'),
      selection = 'single',
      options = list(pageLength = 10,
                     initComplete = JS(
                       "function(settings, json) {",
                       "$(this.api().table().header()).css({'background-color': '#268bd2', 'color': '#fff'});",
                       "}"))))
    
  
  
  
  ranges <- reactiveValues(x = NULL, y = NULL)
  output$contig_plot <- renderPlot(plot.GC.vs.coverage(input, df, ranges))
  
  output$contig_info <- DT::renderDataTable({
    datatable(  
      nearPoints(apply.filters(input, df),
               input$contig_plot_click, threshold = 10, maxpoints = 15,
               addDist = FALSE) %>% select(-ref_start, -ref_finish, -percent_covered),
      options = list(pageLength = 10,
                     initComplete = JS(
                       "function(settings, json) {",
                        "$(this.api().table().header()).css({'background-color': '#268bd2', 'color': '#fff'});",
                       "}"),
                     columnDefs = list(list(className = 'dt-right', targets=c(2,5,9)))), colnames=c('Contig id', 
                                         'Contig Length', 
                                         'Type', 
                                         'VirFinder q-value', 
                                         'Cluster Name', 
                                         'Cluster Representative?', 
                                         '% identity to Representative',
                                         '% GC',
                                         'Mean contig coverage',
                                         'VirSorter Category',
                                         'VirSorter Circular'), rownames=FALSE,class='compact', selection='single') %>%
      DT::formatRound(c(4,7,8,9))
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