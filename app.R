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
library(dplyr)
source("RNA_SeqSaveClass.R")

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
                   #verbatimTextOutput('y12')
               ),
               fluidRow(actionButton("removeGnCnt", "Remove Gene Count Table"))

        )
        ,
        # Main menus showing analysis options
        column(9,
               tabsetPanel(
                   tabPanel("Quality Control", 
                            column(4),
                            ), # Check fastq qc
                   tabPanel("Normalization",
                            column(1),
                            column(11,
                            fluidRow(

                                selectInput(label = "Choose Factors List", inputId = "selectFacTab", choices = NULL)
                                
                                ), 
                            fluidRow(
                                actionButton(label = "New factors table", inputId = "newFacTab", icon = icon("plus")),
                                actionButton(label = "Remove factors table", inputId = "rmFacTab", icon = icon("minus"))

                            ), 
                            fluidRow(
                                h4('Factors (double click on table cell to fill in)'),
                                DT::dataTableOutput('factors'),
                                ),
                            fluidRow(actionButton(label = "Update", inputId = "updateFactors")
                                     ),
                            
                            fluidRow(
                                actionButton(label = "Add Factor", inputId = "addFactor", icon = icon("plus")),
                                actionButton(label = "Add level", inputId = "addLevel", icon = icon("plus")),
                                actionButton(label = "Remove Factor", inputId = "rmFactor", icon = icon("minus")),
                                actionButton(label = "Remove level", inputId = "rmLevel", icon = icon("minus")),
                                     ),
                            fluidRow(
                                h4('Assign Factors'),
                                DT::dataTableOutput('assignfactors'),
                                
                            ),                            
                            ),
                   ),
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
 
    reVals <- reactiveValues(geneSetDes = "", analysisOb = new("RNASeqAnalysis"), geneCntIn = NULL, selectedGnCnts = c(), 
                             factorsTab = tibble(Factors = "", Level1 = "", Level2 = ""), curFacTabDx = 1
                             )
    # analysisOb is an s4 class with slots to store analysis
    # geneCntIn stores the temporary file location of the gene count table when it is loaded
    # selectGnCnts 
    # factorsTab is a tibble with the current factors
    # curFac stores the index of the currently in use factors table stored in the analysisOb
    
    
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
        #output$y12 = renderPrint(input$genCntTabTab_rows_selected)
        
    })
    
    observeEvent(input$removeGnCnt, {
      #  print(input$genCntTabTab_rows_selected)
        reVals$analysisOb <- rmGeneCnts(reVals$analysisOb, as.numeric(input$genCntTabTab_rows_selected)) # remove the counts
       # print(reVals$analysisOb)
        #output$genCntTabTab = DT::renderDataTable(reVals$analysisOb@GeneMeta, server = FALSE, options = list(dom = 't')) # rerender table
    })

### Normalization Tab ============================================================================================== ###
    
    ### Factor table ----------------------------------------------------------------------------------------------- ###
    
    output$factors = DT::renderDataTable(reVals$factorsTab, server = FALSE, options = list(dom = 't'), rownames = F, class = 'cell-border stripe', editable = "cell")
    
    # Add new factor
    observeEvent(input$addFactor, {
        
        reVals$factorsTab[nrow(reVals$factorsTab) + 1, ] <- as.list(rep("", length(reVals$factorsTab)))
        
    })
    
    # Remove selected factor
    observeEvent(input$rmFactor, {
        
        reVals$factorsTab <- reVals$factorsTab[!(1:nrow(reVals$factorsTab) %in% as.numeric(input$factors_rows_selected)), ]

    })
    
    # Add new level
    observeEvent(input$addLevel, {
        #browser()
        reVals$factorsTab[, paste("Level", as.character(ncol(reVals$factorsTab)))] <- rep("", nrow(reVals$factorsTab))
        
    })
    
    # Remove level
    observeEvent(input$rmLevel, {
        
        if (ncol(reVals$factorsTab) > 1){
        reVals$factorsTab <- reVals$factorsTab[, 1:ncol(reVals$factorsTab) -1]
        }
    })
    
    # Update the local factor tbl based on the user input to the table
    observeEvent(input$factors_cell_edit, {
        reVals$factorsTab[input$factors_cell_edit$row, input$factors_cell_edit$col + 1] <- input$factors_cell_edit$value
        #print(reVals$factorsTab)        
    })
    
    
    # Update factors table
    observeEvent(input$updateFactors, {
        #browser()
        reVals$analysisOb <- addFactorsTab(isolate(reVals$analysisOb), isolate(reVals$factorsTab), isolate(reVals$curFacTabDx))
        print(reVals$analysisOb)
        updateSelectInput(session, "selectFacTab", choices = as.character(1:length(reVals$analysisOb@factorsTab)), selected =
                              length(reVals$analysisOb@factorsTab) )
        

    })
    
    # Use drop down meunu to select the factors table
    observeEvent(input$selectFacTab, {
       # browser()
        if (!(input$selectFacTab == "")){
        reVals$curFacTabDx <- as.numeric(input$selectFacTab)
        reVals$factorsTab <- reVals$analysisOb@factorsTab[[reVals$curFacTabDx]]
        }
        print(reVals$factorsTab)



    })
    # Add a new factor table
    observeEvent(input$newFacTab, {
        #browser()
        reVals$curFacTabDx <- length(reVals$analysisOb@factorsTab) + 1
        reVals$factorsTab = tibble(Factors = "", Level1 = "", Level2 = "")

    })
    
    # Remove a factor table
    observeEvent(input$rmFacTab, {
        
        print(reVals$analysisOb)
        reVals$analysisOb <- rmFactorsTab(isolate(reVals$analysisOb), isolate(reVals$factorsTab), isolate(reVals$curFacTabDx))
        print(reVals$analysisOb)
        
        if (length(reVals$analysisOb@factorsTab) > 0){
        updateSelectInput(session, "selectFacTab", choices = as.character(1:length(reVals$analysisOb@factorsTab)), selected = 1)
            reVals$factorsTab <- reVals$analysisOb@factorsTab[[1]]
        }else
        {
            updateSelectInput(session, "selectFacTab", choices = NULL)
            reVals$factorsTab = tibble(Factors = "", Level1 = "", Level2 = "")
            
        }
            
        reVals$curFacTabDx <- 1
        
        
    })
    
    ### Exp Factor table ------------------------------------------------------------------------------------------- ###
    
### ================================================================================================================ ###
### Output Objects ================================================================================================= ###
    
    
    
 
}

# Run the application 
shinyApp(ui = ui, server = server)
