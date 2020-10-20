library(shinydashboard)
library(tidyverse)
library(plotly)
library(DT)

options(shiny.maxRequestSize=50*1024^2)

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "PARcos"),
  dashboardSidebar(
    br(),
    sidebarUserPanel(name = "PAR-CLIP data explorer", subtitle = HTML("&#169 Sarah Arcos, 2020")),
    br(),
    br(),
    fileInput('file1', 'Cluster CSV File',
              accept=c('text/csv', 
                       'text/comma-separated-values,text/plain',
                       '.csv')),
    fileInput('file2', 'SAM File',
              accept=c('text/tsv',
                       'text/tab-separated-values,text/plain',
                       '.sam'))
  ),
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
    ),
    fluidRow(
      box(title = "PAR-CLIP Cluster Table", 
          width = 12,
          status = "primary",
          solidHeader = TRUE,
          DT::dataTableOutput('table'))
    ),
    fluidRow(
      box(title = "Density Plot",
          width = 12, 
          status = "success",
          solidHeader = TRUE,
          plotlyOutput("test"))
    ),
    fluidRow(
      box(title = "Reads in Selected Cluster",
          width = 12,
          status = "success",
          solidHeader = TRUE,
          plotOutput("alignment"))
    )
  )
)

server <- function(input, output) {
  cluster <- reactive({
    if (is.null(input$file1))
      return(NULL)
    clusterFile <- input$file1
    read_csv(clusterFile$datapath, col_names = T)
  })
  
  output$table <- DT::renderDataTable({
    if (is.null(input$file1))
      return(NULL)
    f <- input$file1
    datatable(cluster(), selection='single', class = 'display compact',
              caption = paste("PAR-CLIP Clusters for", gsub(".clusters.csv", "", f$name, )),
              options = list(scrollX = TRUE))
  })
  
  output$test <- renderPlotly({
    if (is.null(input$file1))
      return(NULL)
    
    p <- ggplot(cluster(), aes(color = `Aligned to`, x = T2Cfraction)) +
      geom_density() +
      labs(title = "Distribution of T2C Fraction") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 14, hjust = 0.5),
            legend.position = "bottom")
    
    ggplotly(p)
  })
  
  sam_columns = c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "X12", "X13", "X14")
  
  sam <- reactive({
    if (is.null(input$file2))
      return(NULL)
    samFile <- input$file2
    read_tsv(samFile$datapath, comment = "@", col_names = sam_columns) %>%
      mutate(FLAG = case_when(FLAG == 16 ~ "-",
                              FLAG == 0 ~ "+"))
  })
  
  
  
  output$alignment <- renderPlot({
    req(input$table_rows_selected)
    
    if (is.null(input$file2))
      return(NULL)
    
    #Identify members of cluster
    s=input$table_rows_selected
    peak <- cluster()
    current_row <- peak[s, ]
    samdata <- sam()
    
    cluster_members <- samdata %>%
      filter(POS >= (current_row$Start - 40) &
               POS <= (current_row$End + 40) &
               RNAME==current_row$Chr &
               FLAG==current_row$Strand)
    
    #structure to hold cluster member info for alignment
    align_seq <- data_frame(
      sequence = cluster_members$SEQ,
      start = cluster_members$POS,
      strand = cluster_members$FLAG
    )
    
    #dataframe to hold parent cluster information
    cluster_sequence <- data_frame(
      letters = strsplit(current_row$ClusterSequence, "")[[1]],
      position = current_row$Start:(current_row$Start + length(letters) - 1),
      seq = current_row$ClusterSequence,
      conv = 'O'
    )
    
    #fix the positioning of the minus cluster sequences
    if(current_row$Strand =='-'){
      cluster_sequence$position <- rev(cluster_sequence$position)
    }
    
    #function to form aligned dataframe, called in for loop
    make_sequence_info <- function(sequence, start){
      data_frame(
        letters = strsplit(sequence, "")[[1]], #splits character string into an array of each of its letters.
        position = start:(start + length(letters) - 1), #places each individual letter along the genome.
        seq = sequence,
        conv = 'T'
      )
    }
    
    #set up plotting_data dataframe
    plotting_data <- data_frame(letters = character(), position = integer(), seq = character(), conv = character())
    
    #attach cluster_members to plotting_data, assign values
    for(i in 1:dim(align_seq)[1]){
      
      #complement the minus strand
      if(align_seq$strand[i] == '-'){
        align_seq$sequence[i] <- chartr("ATGC","TACG", align_seq$sequence[i])
      }
      plotting_data <- rbind(
        plotting_data,
        make_sequence_info(align_seq$sequence[i], align_seq$start[i])
      )
    }#end for loop
    
    #add cluster sequence to plot
    #merge instead of rbind so that positions can be compared easily
    plotting_data <- merge(plotting_data, cluster_sequence, by ="position", all=T)
    
    #Conversion assignment, using position comparisons
    for(i in 1:length(plotting_data$position)){
      if(is.na(plotting_data$letters.x[i]) | is.na(plotting_data$letters.y[i])){
        next
      }#this works
      if(plotting_data$letters.x[i] =="C" & plotting_data$letters.y[i] == "T"){
        plotting_data$conv.x[i] <- 'C'
      }
    }
    
    #remove unnecessary caolumns, rename so rbind works
    plotting_data$letters.y=NULL
    plotting_data$seq.y=NULL
    plotting_data$conv.y=NULL
    names(cluster_sequence) <- c("letters.x", "position", "seq.x", "conv.x")
    #browser()
    #plotting_data <- rbind(cluster_sequence, plotting_data)
    #Draw
    ggplot(plotting_data, aes(x = position, y = seq.x,
                              family = "Courier")) +
      coord_fixed() +
      geom_text(aes(label = letters.x, color = conv.x), size = 6) +
      # geom_blank(aes(x = position, y = "Z", label = NULL)) +
      geom_blank(aes(x = position, y = "Y", label = NULL)) +
      geom_blank(aes(x = position, y = "X", label = NULL)) +
      geom_blank(aes(x = position, y = "W", label = NULL)) +
      scale_x_continuous(limits = c(current_row$ModeLocation - 30, current_row$ModeLocation + 30),
                         breaks = c(current_row$ModeLocation - 30,
                                    current_row$ModeLocation - 15,
                                    current_row$ModeLocation,
                                    current_row$ModeLocation + 15,
                                    current_row$ModeLocation + 30)) +
      annotate("text", label = cluster_sequence$letters.x, 
               x = cluster_sequence$position, 
               y = "X", vjust = 0.5, color = "black",
               family = "Courier", size = 8) +
      scale_color_manual(values = c("#ef8a62", "#999999"), name = "T2C Conversion") +
      geom_vline(xintercept=current_row$ModeLocation,
                 color = "blue", size=6, alpha = 0.2) +
      theme(
        legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 16),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(fill = 'white')
      ) +
      labs(title = paste("Reads from cluster: ", current_row$ClusterID, sep = ""),
           subtitle = paste("Aligned to: ", current_row$GeneName, ", ", current_row$`Aligned to`, sep = ""),
           x = paste("Position on ", current_row$Chr, sep = ""))
  })
  
}

shinyApp(ui, server)