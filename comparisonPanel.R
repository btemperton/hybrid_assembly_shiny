load.comparison.panel<-function(cluster_df){
  return(tabPanel("Clusters",
            fluidRow(
              column(2,
                sliderInput(inputId = "cluster_min_length",
                                 label="Minimum cluster length:",
                                 min = 2500,
                                 max = 500000,
                                 step=500,
                                 value = 10000, 500000),
                
               selectInput(inputId = 'hybrid_rep',
                           label='Cluster Type:',
                           choices=c('Choose a type ...' = 0, unique(cluster_df$type))),
               radioButtons(inputId = "contig_color_scheme", "Color scheme:",
                            c("Per Contig" = "pc",
                              "Per Type" = "pt"))
               
              ),
            column(7,   
               DT::dataTableOutput(outputId = 'clusters')
            ),
            column(3,
                   plotOutput(outputId = 'cluster_counts'))),
            fluidRow(
               plotOutput(outputId = 'cluster_summary')))
            )
         
}


filter.clusters<-function(input, cluster_df){
  filtered <- cluster_df %>% filter(representative_length >= input$cluster_min_length)
  return(filtered)
}


prepare.cluster.dataframe<-function(input, cluster_df){
  
  filtered <- filter.clusters(input, cluster_df)
  
  if(input$hybrid_rep >0){
    filtered<-filtered %>% filter(type==input$hybrid_rep)
  }
  
  return(filtered)
}


plot.cluster.summary<-function(input, df, cluster_df){
  if(!is.null(input$clusters_rows_selected)){
    cluster_df<-prepare.cluster.dataframe(input, cluster_df)
    selected <- cluster_df[input$clusters_rows_selected, ]$cluster_name
    filtered <- df %>% filter(cluster_name == selected)
    p <- plot.cluster.alignment(filtered, input)
    return(p)
  }
}

plot.cluster.counts<-function(input, cluster_df){
  
  filtered <- filter.clusters(input, cluster_df)
  p<-ggplot(filtered, aes(type)) +
    geom_bar(aes(fill=type), color='black') +
    xlab('Cluster Type') + 
    ylab('Count') +
    scale_fill_solarized() +
    theme(legend.position='none',
          axis.text.x = element_text(angle = 45, hjust = 1, size=12))
  
  return(p)
}


plot.cluster.alignment<-function(cluster_df, input){
  cluster_name <-unique(cluster_df$cluster_name)
  cluster<- cluster_df %>% 
    select(contig_id, contig_len, ref_start,ref_finish, is_cluster_representative, contig_type) %>%
    arrange(desc(contig_len))
  
  
  representative <- cluster %>% filter(is_cluster_representative=='Y')
  representative$ref_start<-1
  representative$ref_finish<-representative$contig_len
  
  representative$ref_start<-as.numeric(representative$ref_start)
  representative$ref_finish<-as.numeric(representative$ref_finish)
  
  members<-cluster %>% filter(is_cluster_representative=='N')
  
  members_start<- members %>% select(-ref_finish) %>%
    transform(ref_start = strsplit(as.character(ref_start), split='|', fixed=TRUE)) %>%
    unnest()
  members_finish <- members %>% select(-ref_start) %>% 
    transform(ref_finish = strsplit(as.character(ref_finish), split='|', fixed=TRUE)) %>%
    unnest()
  
  members_start$ref_finish <- members_finish$ref_finish
  members_start$ref_start<-as.numeric(members_start$ref_start)
  members_start$ref_finish<-as.numeric(members_finish$ref_finish)
  
  
  #members_final <- members_start %>% group_by(contig_id) %>%
  #  filter(!(ref_start > min(ref_start) & ref_finish < max(ref_finish)))
  
  members_final<-members_start %>% bind_rows(representative) %>%
    mutate(ranking = dense_rank(contig_len), 
           full_length=max(contig_len))
  
  
  colourCount <- length(unique(members_final$contig_id))
  
  p<- ggplot(members_final) 
  if(input$contig_color_scheme=='pc'){
      p <- p + 
        geom_rect(aes(xmin=ref_start, 
                      xmax=ref_finish,
                      ymin=ranking-0.5, 
                      ymax=ranking+0.5,
                      fill=contig_id), size=1, color='black') + 
        scale_fill_manual(values = colorRampPalette(solarized_pal()(8))(colourCount)) 
  } else {
    p <- p + 
      geom_rect(aes(xmin=ref_start, 
                    xmax=ref_finish,
                    ymin=ranking-0.5, 
                    ymax=ranking+0.5,
                    fill=contig_type), size=1, color='black') + xlab('locus of representative (bp)') +
      scale_fill_solarized(accent='green') + 
      theme(legend.position='bottom',
            plot.title = element_text(face="bold", size=16, hjust=0),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
  }
  
  p <- p + xlab('locus of representative (bp)') +
        ggtitle(paste(cluster_name, 'alignments')) + 
    theme(legend.position='none',
          plot.title = element_text(face="bold", size=16, hjust=0),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  return(p)
}
