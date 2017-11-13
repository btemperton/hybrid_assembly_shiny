load.comparison.panel<-function(cluster_df){
  return(tabPanel("Clusters",
            fluidRow(column(9, h1("Cluster Analysis"), HTML("
<p>Here, you can examine the content of the clusters. Clusters were created by combining all short and hybrid scaffolds and then
clustering them at 95% identity and 80% length into 'viral OTUs'. 
You can click on a cluster to show how its members aligned against the reference sequence (the one at the top). The alignment plot is ordered
from top to bottom on scaffold length. The alignments can be colored by selecting an appropriate scheme on the left-hand side. You
can also filter the number of clusters by size and whether any members were identified as viral, either by VirSorter, or by VirFinder.
There are several types of cluster:</p>
<ul>
<li><b>Hybrid Extension</b> - The longest read in a cluster is from the hybrid assembly and typically denotes a cluster where short-read scaffolds 
have been increased in length (e.g. Cluster_1), or scaffolded by the long reads (e.g. Cluster_2 and Cluster_4)</li>
<li><b>Identical</b> - The cluster comprises of two members, each of identical length. In the alignment plot, these will appear as one block</li>
<li><b>Short longer than hybrid</b> - The longest member of the cluster is from the short-read assembly</li>
<li><b>Singleton S/H</b> - The cluster consists of only one member, either a short or hybrid scaffold</li>
<li><b>Unclassified</b> - A cluster that did not meet one of the above criteria</li>
</ul><hr/>
                                        "))),
            fluidRow(
              column(2,
                sliderInput(inputId = "cluster_min_length",
                                 label="Minimum cluster length:",
                                 min = 2500,
                                 max = 500000,
                                 step=500,
                                 value = 10000, 500000),
                
                checkboxInput(inputId = 'circ_rep',
                              label='Restrict to circular representatives'),
                
               selectInput(inputId = 'hybrid_rep',
                           label='Cluster Type:',
                           choices=c('Choose a type ...' = 0, unique(cluster_df$type))),
               selectInput(inputId = 'viral_ratio',
                           label='Viral Evidence:',
                           choices=c('Choose a type ...' = 'empty', 'Some viral members' = 'svm', 'All viral' = 'av', 'No viral' = 'nv')),
               radioButtons(inputId = "contig_color_scheme", "Color scheme:",
                            c("Per Contig" = "pc",
                              "Per Type" = "pt", 
                              "Per viralness" = "pv"))
               
              ),
            column(7,   
               DT::dataTableOutput(outputId = 'clusters')
            ),
            column(3,
                   plotOutput(outputId = 'cluster_counts'))),
            fluidRow(
               plotOutput(outputId = 'cluster_summary')),
            
            fluidRow(
              DT::dataTableOutput(outputId = 'cluster_rep_genes')
            )
            
        )
    )
  
     
}


filter.clusters<-function(input, cluster_df){
  filtered <- cluster_df %>% filter(representative_length >= input$cluster_min_length)
  if(input$viral_ratio=='svm'){
    filtered <- filtered %>% filter(viral_ratio >0)
  } else if (input$viral_ratio=='av'){
    filtered <- filtered %>% filter(viral_ratio == 1)
  } else if (input$viral_ratio=='nv'){
    filtered <- filtered %>% filter(viral_ratio == 0)
  }
  
  if(input$circ_rep){
    filtered <- filtered %>% filter(is_circular ==1)
  }
  
  return(filtered)
}


prepare.cluster.dataframe<-function(input, cluster_df){
  
  filtered <- filter.clusters(input, cluster_df)
  
  if(input$hybrid_rep >0){
    filtered<-filtered %>% filter(type==input$hybrid_rep)
  }
  
  return(filtered %>% select(-is_circular))
}

prepare.protein.dataframe<-function(input, protein_df){
  filtered <- protein_df
  if(!is.null(input$clusters_rows_selected)){
    cluster_df<-prepare.cluster.dataframe(input, cluster_df)
    selected <- cluster_df[input$clusters_rows_selected, ]$representative
    filtered <- filtered %>% filter(contig_id == selected)
  }
  return(filtered %>% select(-gene_id))
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
    select(contig_id, contig_len, ref_start,ref_finish, is_cluster_representative, contig_type, virfinder_qvalue, virsorter_category) %>%
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
           full_length=max(contig_len),
           is_viral=(virfinder_qvalue <= 0.05 | virsorter_category %in% c(1,2,4,5)))
    
  colourCount <- length(unique(members_final$contig_id))
  
  p<- ggplot(members_final) 
  if(input$contig_color_scheme=='pc'){
      p <- p + 
        geom_rect(aes(xmin=ref_start, 
                      xmax=ref_finish,
                      ymin=ranking-0.5, 
                      ymax=ranking+0.5,
                      fill=contig_id), size=1, color='black') + 
        scale_fill_manual(values = colorRampPalette(solarized_pal()(8))(colourCount))  + 
        theme(legend.position='none',
              plot.title = element_text(face="bold", size=16, hjust=0),
              axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())
  } else if (input$contig_color_scheme=='pt'){
    p <- p + 
      geom_rect(aes(xmin=ref_start, 
                    xmax=ref_finish,
                    ymin=ranking-0.5, 
                    ymax=ranking+0.5,
                    fill=contig_type), size=1, color='black') + xlab('locus of representative (bp)') +
      scale_fill_solarized(accent='green')  + 
      theme(legend.position='bottom',
            plot.title = element_text(face="bold", size=16, hjust=0),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
  } else {
    p <- p + 
      geom_rect(aes(xmin=ref_start, 
                    xmax=ref_finish,
                    ymin=ranking-0.5, 
                    ymax=ranking+0.5,
                    fill=is_viral), size=1, color='black') + xlab('locus of representative (bp)') +
      scale_fill_manual(values = c('skyblue4', 'skyblue1'))  + 
      theme(legend.position='bottom',
            plot.title = element_text(face="bold", size=16, hjust=0),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
  }
  
  p <- p + xlab('locus of representative (bp)') +
        ggtitle(paste(cluster_name, 'alignments'))
  return(p)
}
