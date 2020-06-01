#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyFiles)
library(DT)

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    fluidRow(
        # Left menu shwing analysis and loaded files
        column(3, style = "background-color:lightblue;", 
               fluidRow("Analysis Options"),
               fluidRow(actionButton("newAnalysisactBut","Start New Analysis")),
               fluidRow(fileInput(label = "Load Analysis", inputId = "loadAnalysis")),
               fluidRow(shinyDirButton("dir", "Select Save Directory", "Upload"),
                        verbatimTextOutput("dir", placeholder = TRUE)  
               ),
               
               fluidRow(textInput(label = "Save Analysis Filename", inputId = "saveFn")),
               fluidRow(actionButton("saveAnalysis", "Save Analysis")),
               #fluidRow(tableOutput("cntTables"))
               fluidRow(fileInput(label = "Load Gene Count table", inputId = "loadCntTab")),

               # table showing count data
               fluidRow(h4("Loaded Gene Count tables")),
               #fluidRow(DTOutput("genCntTbl")),
               fluidRow(
                   h4('Uploaded Gene Count tables'),
                   DT::dataTableOutput('genCntTabTab'),
                   verbatimTextOutput('y12')
               ),
               fluidRow(actionButton("removeGnCnt", "Remove Gene Count Table"))

        )
        ,
        # Main menus showing analysis options
        column(9,
               tabsetPanel(
                   tabPanel("Quality Control", 
                            column(4),
                            # column(3, sliderInput("bins",
                            #                       "Number of bins:",
                            #                       min = 1,
                            #                       max = 50,
                            #                       value = 30))
                            ), # Check fastq qc
                   tabPanel("Normalization"),
                   tabPanel("Data Exploration"),
                   tabPanel("Differential Analysis"), # Check overlaps
                   tabPanel("Genelists/Sequences"),
                   tabPanel("Pathway Analysis"),
                   tabPanel("")
                   
               ),
               
        )
    )
)

# Server code ========================================================================================== ###
# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    
    
    # Reactive values ---------------------------------------------------------------------------------- # 
    reVals <- reactiveValues(geneSetDes = "", analysisOb = new("RNASeqAnalysis"), geneCntIn = NULL, selectedGnCnts = c())
    
    # Load/create analysis objects --------------------------------------------------------------------- # 
    # Create new analysis object to store data and results
    observeEvent(input$newAnalysisactBut, {
        
        analysisOb <- new("RNASeqAnalysis")
        # Need to clear all the variables
        
    })
    
    # Load a previous analysis
    observeEvent(input$loadAnalysis, {
        req(input$loadAnalysis$datapath)
        print(input$loadAnalysis$name)

        # Import the rdata object and assign the analysis object to analysisOb
        
        # Load fastq results
        
        # load merging results
        
    })
    
    # Save analysis
    shinyDirChoose(
        input,
        'dir',
        roots = c(home = "/data/RNASeqAnalysis"),
        #filetypes = c('', 'txt', 'bigWig', "tsv", "csv", "bw")
    )
    
    global <- reactiveValues(datapath = getwd())
    
    dir <- reactive(input$dir)
    
    output$dir <- renderText({
        "Select Save Path"
    })
    
    observeEvent(ignoreNULL = TRUE,
                 eventExpr = {
                     input$dir
                 },
                 handlerExpr = {
                     if (!"path" %in% names(dir())) return()
                     home <- normalizePath("~")
                     global$datapath <-
                         file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
                     #print(global$datapath)
                     output$dir <- renderText({
                         global$datapath
                     }) 
                 })

    
    observeEvent(input$saveAnalysis, {
        
        #saveRDS(analysisOb, )
        
    })
    
    
    # Modal dialog box for writing a description ------------------------------------------------- #
    dataModal <- function(failed = FALSE) {
        modalDialog(
            textInput("description", "Tag Data", 
                      placeholder = ""
            ),
            span('Write a brief description of data for your reference'),
            if (failed)
                div(tags$b("Invalid name of data object", style = "color: red;")),
             
            footer = tagList(
                modalButton("Cancel"),
                actionButton("ok", "OK")
            )
        )
    }
 
    # Load new datasets --------------------------------------------------------------------------- #
    # Load gene count tables
    observeEvent(input$loadCntTab, {
        #print(reVals$analysisOb)
        #print(input$loadCntTab$name)
        reVals$geneCntIn <- read.csv(input$loadCntTab$datapath, skip = 4, sep = "\t")
        colnames(reVals$geneCntIn) <- c("ENSEMBL_ID", "+", "-", "All")
        print(head(reVals$geneCntIn))
        #print(analysisOb)
        
        # show description dialog
        showModal(dataModal())

    })
    # close dialog box on clck ok
    observeEvent(input$ok, {
        #print(input$description)
        reVals$geneSetDes <- input$description
        
        # Write gene count table and meta data to analysisOb
        reVals$analysisOb <- newGeneCnts(isolate(reVals$analysisOb), reVals$geneCntIn, data.frame(Description = reVals$geneSetDes, FileName = input$loadCntTab$name, stringsAsFactors = F))
        reVals$geneSetDes <- "" # reset the description
        reVals$geneCntIn <- NULL
        print(reVals$analysisOb)
        removeModal()

        # Update the gene counts table
        output$genCntTabTab = DT::renderDataTable(reVals$analysisOb@GeneMeta, server = FALSE, options = list(dom = 't'))
        output$y12 = renderPrint(input$genCntTabTab_rows_selected)
        
    })
    
    observeEvent(input$removeGnCnt, {
      #  print(input$genCntTabTab_rows_selected)
        reVals$analysisOb <- rmGeneCnts(reVals$analysisOb, as.numeric(input$genCntTabTab_rows_selected)) # remove the counts
       # print(reVals$analysisOb)
        output$genCntTabTab = DT::renderDataTable(reVals$analysisOb@GeneMeta, server = FALSE, options = list(dom = 't')) # rerender table
    })

    
### ================================================================================================================ ###
### Output Objects ================================================================================================= ###
    
    
    
 
}

# Run the application 
shinyApp(ui = ui, server = server)
