load.comparison.panel<-function(cluster_df){
  return(tabPanel("Clusters",
         
         selectInput(inputId = 'hybrid_rep',
                     label='Cluster Type:',
                     choices=c('Choose a type ...' = 0, unique(cluster_df$type))),
         
         DT::dataTableOutput(outputId = 'clusters'),
         
         plotOutput(outputId = 'cluster_summary')))
}

prepare.cluster.dataframe<-function(input, cluster_df){
  
  if(input$hybrid_rep >0){
    cluster_df<-cluster_df %>% filter(type==input$hybrid_rep)
  }  
  
  return(cluster_df)
}