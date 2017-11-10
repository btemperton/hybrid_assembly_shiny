load.comparison.panel<-function(cluster_df){
  return(tabPanel("Clusters",
         
         selectInput(inputId = 'hybrid_rep',
                     label='Cluster Type:',
                     choices=unique(cluster_df$type)),
         
         DT::dataTableOutput(outputId = 'clusters'),
         
         plotOutput(outputId = 'cluster_summary')))
}