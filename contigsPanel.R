load.contigs.panel<-function(){
  return(tabPanel("Contigs",
           fluidRow(column(9, h1("Contig Analysis"), p("Here, we can explore the short read only scaffolds (S) and hybrid scaffolds (H) from a marine virome
                                  from the Western English Channel. You can filter the data using the controls on the left hand side
                                  of the page. Clicking on a point on the graph will show the 15 nearest scaffolds to that point in the
                                  table below. You can also draw a box on the graph and double-click to zoom in.
                                 "))),
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
                      DT::dataTableOutput(outputId = "contig_info"))))))
  
  
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

plot.GC.vs.coverage<-function(input, df, ranges){
  filtered<-apply.filters(input, df)
  p<-ggplot(filtered, aes(x=GC, y=mean_coverage)) +
    geom_point(aes(fill=virsorter_category, size=contig_len), color='black', pch=21, stroke=1, alpha=0.5) +
    xlab('GC %') + 
    ylab('Mean Coverage') +
    scale_y_log10() +
    scale_colour_solarized() + 
    scale_size_continuous(range = c(1,12)) +
    coord_cartesian(xlim = ranges$x, ylim = ranges$y)
  return(p)
}

prepare.contig.dataframe<-function(input, df){
  filtered<-apply.filters(input, df)
  
  return(filtered)
}