load.contigs.panel<-function(){
  return(tabPanel("Contigs",
           fluidRow(
             
             column(2,
                    sliderInput(inputId = "min_length",
                                label="Minimum contig length:",
                                min = 2500,
                                max = 500000,
                                step=500,
                                value = 10000, 500000),
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
             
             column(7,
                    
                    plotOutput(outputId = 'contig_plot',
                               dblclick = "contig_plot_dblclick",
                               click = "contig_plot_click",
                               brush = brushOpts(id = "contig_plot_brush", resetOnNew = TRUE))),
             
             column(3,
                    plotOutput(outputId = 'summary_contig_plot')),
             
             fluidRow(
               column(12,
                      tableOutput(outputId = "contig_info"))))))
  
  
}

plot.contig.summary<-function(input, df){
  filtered<-apply.filters(input, df)
  p<-ggplot(filtered, aes(contig_type)) +
    geom_bar(aes(fill=contig_type), color='black') +
    xlab('Contig Type') + 
    ylab('Count') +
    scale_fill_solarized() +
    theme(legend.position='none')
    
  return(p)
}