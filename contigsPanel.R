load.contigs.panel<-function(){
  return(tabPanel("Contigs",
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
                      tableOutput(outputId = "contig_info"))))))
  
  
}