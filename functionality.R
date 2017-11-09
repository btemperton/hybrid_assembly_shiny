apply.filters<-function(input, df){
  
  filtered <- df %>% filter(contig_len >= input$min_length)
  
  if(input$VirFinder){
    filtered <- filtered %>% filter(virfinder_qvalue <=0.05)
  }
  
  if(input$VirSorter){
    viral_values = c(1,2,4,5)
    filtered <- filtered %>% filter(virsorter_category %in% viral_values)
    if (input$VirSorter_circular){
      filtered <- filtered %>% filter(is_circular ==1)
    }
  }
  
  if (input$read_type=='short'){
    filtered <-filtered %>% filter(contig_type=='S')
  } else if (input$read_type=='hybrid'){
    filtered <-filtered %>% filter(contig_type=='H')
  }
  

  return(filtered)
}

prepare.contig.dataframe<-function(input, df){
  filtered<-apply.filters(input, df)
  
  return(filtered)
}

prepare.cluster.dataframe<-function(input, cluster_df){
  
  if(!is.null(input$hybrid_rep)){
    cluster_df<-cluster_df %>% filter(type==input$hybrid_rep)
  }  
  
  return(cluster_df)
}


plot.GC.vs.coverage<-function(input, df, ranges){
  filtered<-apply.filters(input, df)
  p<-ggplot(filtered, aes(x=GC, y=mean_coverage, 
                          color=virsorter_category,
                          size=contig_len)) +
    geom_point() +
    xlab('GC %') + 
    ylab('Mean Coverage') +
    scale_y_log10() +
    scale_color_nejm() + 
    coord_cartesian(xlim = ranges$x, ylim = ranges$y)
  return(p)
}

plot.cluster.summary<-function(input, df, cluster_df){
  if(!is.null(input$clusters_rows_selected)){
    filtered<-apply.filters(input, df)
    cluster_df<-prepare.cluster.dataframe(input, cluster_df)
    selected <- cluster_df[input$clusters_rows_selected, ]$cluster_name
    filtered <- filtered %>% filter(cluster_name == selected)
    p <- plot.cluster.alignment(filtered)
    return(p)
  }
}

plot.cluster.alignment<-function(cluster_df){
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

  p<- ggplot(members_final) + 
    geom_segment(aes(x=ref_start, 
                     xend=ref_finish,
                     y=ranking, 
                     yend=ranking,
                     color=contig_id), size=8) + xlab('locus of representative (bp)') +
    ggtitle(paste(cluster_name, 'alignments')) + 
    scale_colour_manual(values = colorRampPalette(solarized_pal()(8))(colourCount)) + 
    theme(legend.position='none',
          plot.title = element_text(face="bold", size=16, hjust=0),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
   
  return(p)
}
