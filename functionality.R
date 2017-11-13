apply.filters<-function(input, df){
  
  filtered <- df %>% filter(contig_len >= input$min_length)
  
  if(input$VirFinder){
    filtered <- filtered %>% filter(virfinder_qvalue <=0.05)
  }
  
  if(input$VirSorter){
    viral_values = c(1,2,4,5)
    filtered <- filtered %>% filter(virsorter_category %in% viral_values)
    
  }
  
  if (input$VirSorter_circular){
    filtered <- filtered %>% filter(is_circular ==1)
  }
  
  if (input$read_type=='short'){
    filtered <-filtered %>% filter(contig_type=='S')
  } else if (input$read_type=='hybrid'){
    filtered <-filtered %>% filter(contig_type=='H')
  }
  

  return(filtered)
}




