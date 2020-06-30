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
library(ggplot2)
library(plotly)
source("RNA_SeqSaveClass.R")
source("analysisFunctions.R")
source("nrmCntResults.R")
source("plottingFunctions.R")
options(shiny.maxRequestSize = 50*1024^2)
# Define UI for application that draws a histogram
ui <- fluidPage(
    
    fluidRow(
        # Left menu shwing analysis and loaded files
        column(3, style = "background-color:lightblue;", 
               fluidRow(h4("Analysis Options")),
               fluidRow(actionButton("newAnalysisactBut","Start New Analysis")),
               fluidRow(fileInput(label = "Load Analysis", inputId = "loadAnalysis")),
               fluidRow(h4('Save Analysis')),
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
                   
                   
                   DT::dataTableOutput('genCntTabTab'),
                   #verbatimTextOutput('y12')
               ),
               fluidRow(actionButton("removeGnCnt", "Remove Gene Count Table")),
               fluidRow(h4('Factors tables'), DT::dataTableOutput('ldedGnFacsTab')),
              

        )
        ,
        # Main menus showing analysis options
        column(9,
               tabsetPanel(
                   ## QC tab panel ---------------------------------------------------------------------------## 
                   tabPanel("Quality Control", 
                            column(4),
                            ), # Check fastq qc
                   
                   ## Factors tab panel ---------------------------------------------------------------------------## 
                   tabPanel("Factors",
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
                                h4('Factors (double click on table cell to fill in, put control factors in level 1)'),
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
                                h4('Experiment sample table'),
                                textInput("nameExpSamTab", "Experiment Sample table name", placeholder = "Experiment name"),
                                selectInput(inputId = "selectExpSmp", label = "Select Experiment sample table", choices = NULL),
                                DT::dataTableOutput('assignfactors'),
                                
                            ),
                            fluidRow(
                                
                                actionButton("newExpSmpTab", "New Sample Exp table"), # start a new assign factors table
                                actionButton("saveExpSmpTab", "Save Sample Exp table"),
                                actionButton("rmExpSmpTab", "Remove Sample Exp table")
                                
                            )
                            ),
                   ),
                   ## Normalization tab panel ---------------------------------------------------------------------------## 
                   tabPanel("Normalization", column(1),
                            column(11,
                            h3("Normalization using DESeq2 median of ratios (we will add other methods if needed)"),
                            fluidRow(textInput("nrmedCntsName", "Name normed counts"),
                                     h4("Select Design Factors"),
                                     DT::dataTableOutput("selectDesignFactors")),
                            fluidRow(
                            selectInput(inputId = "selectExpSmpNrm", label = "Select Experiment sample table for normalization", choices = NULL),
                            numericInput(inputId = "rmLowCnts", label = "Remove genes with total counts less than:", value = 10),
                            actionButton("newDESeqNrmCnts", "New DESeq2 Norm"), actionButton("rmNrmCnts", "Remove normalized counts"), 
                            
                            ),
                            #fluidRow(,
                            fluidRow(h3("Normed Count Matrices"),
                                     DT::dataTableOutput("nrmCntsDT")),
                            
                            )
                            
                            ),                        
                   
                   ## Data Exploration tab panel ---------------------------------------------------------------------------## 
                   tabPanel("Data Exploration"),
                   
                   ## Differential Analysis tab panel ---------------------------------------------------------------------------## 
                   tabPanel("Differential Analysis", # Check overlaps
                   column(1),
                   column(3,
                          fluidRow(actionButton(label = "Run DESeq2", inputId = "runDESeq2")),
                          fluidRow(selectInput(label = "Select Normalized Counts", inputId = "selectNrmCntsDiff", choices = NULL), 
                                   selectInput(label = "Select Contrast Factor", inputId = "selectNrmCntsCntrst", choices = NULL),
                                   selectInput(label = "Select constrast condition", inputId = "selectNrmCntsCnd", choices = NULL)),
                          fluidRow(actionButton(label = "Get DESeq2 results", inputId = "deseqResAct"), checkboxInput(label = "LFC", inputId = "LFCSelect")),
                          fluidRow(actionButton(label = "MA Plot", inputId = "plotMAAct"), actionButton(label = "LFC MA Plot", inputId = "lfcplotMAAct")), 
                          fluidRow(actionButton(label = "Volcano Plot", inputId = "plotVolcanoAct"), actionButton(label = "LFC Volcano Plot", inputId = "lfcplotVolcanoAct")), 
                          fluidRow(),
                          fluidRow(actionButton(label = "Save MA Plot", inputId = "saveMAPltAct"), actionButton(label = "Save LFC MA Plot", inputId = "saveLFCMAPltAct")), 
                          fluidRow(actionButton(label = "Save Volcano Plot", inputId = "saveVolcanoPltAct"), actionButton(label = "Save LFC Volcano Plot", inputId = "saveLFCVolcanoPltAct")), 
                   ),

                          column(7, 
                          fluidRow(mainPanel(plotlyOutput("MAPlot"))),
                          fluidRow(mainPanel(plotlyOutput("volPlot"))),
                          ),
                          
                          
                   
                   ),
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
   #output$MAPlot <- renderPlotly(plot_ly(data = iris, x = ~Sepal.Length, y = ~Petal.Length))
 
    # Reactive values ---------------------------------------------------------------------------------- # 
 
    reVals <- reactiveValues(geneSetDes = "", analysisOb = new("RNASeqAnalysis"), geneCntIn = NULL, selectedGnCnts = c(), 
                             factorsTab = tibble(Factors = "", Level1 = "", Level2 = ""), curFacTabDx = 1,
                             assignFactorsTab = tibble(`Gene Count Table` = "Nothing Selected"),
                             nrmedCntsTab = tibble(Name = "", `Exp Smp List` = "", `Norm Method` = "", `Design Factors` = "")
                             )
    
    #assignFactorsTab <- tibble(Gene_Count_tab = "Nothing Selected") # tibble ot store the assign factors table
    
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
        #browser()
        
        
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
                     home <- normalizePath("/data/RNASeqAnalysis")
                     global$datapath <-
                         file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
                     #print(global$datapath)
                     output$dir <- renderText({
                         global$datapath
                     })
                 })

    # Save the analysis
    observeEvent(input$saveAnalysis, {
        #browser()
        saveReVals <- reVals$analysisOb
        
        if (input$saveFn == ""){save(saveReVals, file = paste0(global$datapath, "/", "RNASeqAnalysis.rData"))}
        else {save(saveReVals, file = paste0(global$datapath, "/",input$saveFn, ".rData"))}

    })
    
    # load analysis
    observeEvent(input$loadAnalysis,
                 {
                     
                     #browser()
                     if ( is.null(input$loadAnalysis)) return(NULL)
                     inFile <- input$loadAnalysis
                     file <- inFile$datapath
                     # load the file into new environment and get it from there
                     e = new.env()
                     name <- load(file, envir = e)
                     reVals$analysisOb <- e[[name]]
                     
                     #browser()
                     # Add the directory and file name to the save directory
                     updateTextInput(session, "saveFn", value = gsub(".rData", "", inFile$name)) 
                     
                     # Recreate the tables etc
                     # Gene counts table table
                     output$genCntTabTab = DT::renderDataTable(reVals$analysisOb@GeneMeta, server = FALSE, options = list(dom = 't'))
                     
                     # Factors tables and drp downs
                     #browser()
                     if (length(reVals$analysisOb@factorsTab) > 0){
                     updateSelectInput(session, "selectFacTab", choices = as.character(1:length(reVals$analysisOb@factorsTab)), selected =
                                           length(reVals$analysisOb@factorsTab) )
                        }
                     # Sample experiment table
                     updateSelectInput(session, "selectExpSmp", choices = names(reVals$analysisOb@ExpSmpTab), selected = names(reVals$analysisOb@ExpSmpTab)[[1]]) # update drop down
                     # read in values
                     
                     # Normalized reads table
                     #browser()
                     if (length(reVals$analysisOb@NrmCnts) > 0){
                       
                       
                     reVals$nrmedCntsTab <- tibble(Name = names(reVals$analysisOb@NrmCnts), `Exp Smp List` = unlist(reVals$analysisOb@NrmCntsExpSmp), 
                                                   `Norm Method` = unlist(lapply(1:length(reVals$analysisOb@NrmCnts), function(x) reVals$analysisOb@NrmCnts[[x]]@NrmMethod)), 
                                                   `Design Factors` = unlist(lapply(1:length(reVals$analysisOb@NrmCnts), function(x) reVals$analysisOb@NrmCnts[[x]]@design)))
                     }
                     
                     # DifferentiaL analysis tab
                     # Populate drop downs
                     if (length(reVals$analysisOb@NrmCnts) > 0){
                     updateSelectInput(session, "selectNrmCntsDiff", choices = names(reVals$analysisOb@NrmCnts), selected = names(reVals$analysisOb@NrmCnts)[[1]])
                       # p = reVals$analysisOb@factorsTab[[reVals$analysisOb@ExpSmpFactorsTab[[reVals$analysisOb@NrmCnts[[1]]@ExpSampNm]]]]  # this is evil
                       # # 
                       # updateSelectInput(session, "selectNrmCntsCntrst",
                       #                   choices = p$Factors,  # this is evil
                       #                   selected = p$Factors[1])
                       # 
                       # 
                       # 
                       # updateSelectInput(session, "selectNrmCntsCnd",
                       #                   choices = as.character(p[p$Factors == p$Factors[1],][2:ncol(p)]),
                       #                   select = as.character(p[p$Factors == p$Factors[1],][2:ncol(p)])[1])
                     }
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
         # check if the data has 4 columns
        if (ncol(reVals$geneCntIn) != 4){
          
          showModal(modalDialog(title = "Warning", "Data does not have 4 columns, cannot load"))
          
        }
        else{
        colnames(reVals$geneCntIn) <- c("ENSEMBL_ID", "+", "-", "All")
        #print(head(reVals$geneCntIn))
        #print(analysisOb)
        
        # show description dialog
        showModal(dataModal())
        }
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
        #browser()
        if (!(input$selectFacTab == "")){
        reVals$curFacTabDx <- as.numeric(input$selectFacTab)
        reVals$factorsTab <- reVals$analysisOb@factorsTab[[reVals$curFacTabDx]]
        reVals$analysisOb@State["factorsTab", "selected"] <- input$selectFacTab
        }
        #print(reVals$factorsTab)

        

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
    
    # Helper function for making dropdown multi row
    shinyInput = function(FUN, len, id,...) {
        inputs = character(len)
        for (i in seq_len(len)) {
            inputs[i] = as.character(FUN(paste0(id, i), label = NULL, ...))
        }
        inputs
    }
    
    # Helper function for making dropdown single row
    shinyInput1 = function(FUN, idx, id,...) {
      inputs = character(idx)
      #for (i in seq_len(len)) {
      inputs = as.character(FUN(paste0(id, idx), label = NULL, ...))
      #}
      inputs
    }
    
    # helper function for reading selections
    shinyValue = function(id, len) {
        unlist(lapply(seq_len(len), function(i) {
            value = input[[paste0(id, i)]]
            if (is.null(value)){
                NA
            } else {
                value
            }
        }))
    }
    
    # helper function to read drop down menu tabs from ExpSamp DT into tibble (with correct options)
    readDropDownDTtoTib = function(dt){
      
      for (dx in 2:ncol(dt)){
        
        dt[colnames(dt[dx])] <- shinyValue(colnames(dt[dx]),nrow(dt))
        
      }
   
      return(dt)
      
    }
    
    # helper function to read drop down menu tabs from tibble DT into ExpSamp (with correct options)
    readTibToDropDown = function(tib, facs){
      #browser()
      # tib = the selected factors
      # facs = the associated factors table
       
    for (cx in 1:nrow(tib)){

      for (dx in 1:nrow(facs[,1])){
      sel <- tib[[cx, facs[[dx, 1]]]] # get the selected option
        tib[cx, facs[[dx, 1]]] <-
          shinyInput1(selectInput,
                     cx,
                     facs[[dx,1]],
                     choices=as.character(facs[dx, 2:ncol(facs)]), selected = sel) # overwrite the assignFactorsTab

      }
    }

      reVals$assignFactorsTab <- tib
   # browser()
      # render
      output$assignfactors <- DT::renderDataTable(isolate(tib), server = FALSE, escape = FALSE, selection = 'none', options = list(
        dom = 't', paging = FALSE, ordering = FALSE,
        preDrawCallback = JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
        drawCallback = JS('function() { Shiny.bindAll(this.api().table().node()); } ')
      ),rownames=FALSE)

      
    }
 
    # Update the assign factors table when update assign factors button is pushed
    observeEvent(input$newExpSmpTab,{
        #browser()
        
      ## You need to remove the ui drop downs from the currently displayed table. If you don't and the new table has factor
      ## names which are the same as the old one, it wil use the selected values from the old table and not update
      
      # Get the factor names from the currently displayed table
      if (length(reVals$analysisOb@ExpSmpTab) > 0){
      curTb <- reVals$analysisOb@ExpSmpTab[[input$selectExpSmp]] 
      curFacN <- ncol(curTb) - 1
      curFacs <- colnames(curTb[, 2:(curFacN + 1)])
      curGnCntN <- nrow(curTb)
      
      # Loop through factors and remove ui
      
      for (dx in 1:curFacN){
        
        for(ex in 1:curGnCntN){
          
          removeUI(selector = paste0("#", curFacs[dx], ex), immediate = T)
          
        }
   
      }
      }
      # Create a new empty table
        reVals$assignFactorsTab <- tibble(`Gene Count tab` = reVals$analysisOb@GeneMeta[input$genCntTabTab_rows_selected, 1])
        
        
        
        for (dx in 1:nrow(reVals$factorsTab[,1])){ # add factor columns
            
            reVals$assignFactorsTab[, reVals$factorsTab[[dx, 1]]] <-
                                          shinyInput(selectInput,
                                                               length(input$genCntTabTab_rows_selected),
                                                               reVals$factorsTab[[dx,1]],
                                                               choices=as.character(reVals$factorsTab[dx, 2:ncol(reVals$factorsTab)])) # overwrite the assignFactorsTab

        }
        
         
    output$assignfactors <- DT::renderDataTable(isolate(reVals$assignFactorsTab), server = FALSE, escape = FALSE, selection = 'none', options = list(
        dom = 't', paging = FALSE, ordering = FALSE,
        preDrawCallback = JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
        drawCallback = JS('function() { Shiny.bindAll(this.api().table().node()); } ')
    ),rownames=FALSE)
    
    })
    
    # save an example table
    observeEvent(input$saveExpSmpTab, {
      #browser()
      
      # Loop through and remove the UI objects
      
      ## I'm not sure why, but upon exiting the save the reVals$assignFacTab object stops being a 
      ## tibble and becomes a 'shiny.render.function'. In this case we can't save the table. This
      ## bit of code converts it back so you can change an already saved table
      
      # At the moment you can't change an exp factors table
      if (any(class(reVals$assignFactorsTab) == "shiny.render.function")){
        reVals$assignFactorsTab <- tibble(`Gene Count tab` = reVals$analysisOb@ExpSmpTab[[input$selectExpSmp]][1])
      
      reVals$factorsTab <- reVals$analysisOb@factorsTab[[reVals$analysisOb@ExpSmpFactorsTab[[input$selectExpSmp]]]]
      
      for (dx in 1:nrow(reVals$factorsTab[,1])){ # add factor columns
        
        reVals$assignFactorsTab[, reVals$factorsTab[[dx, 1]]] <-
          shinyInput(selectInput,
                     length(reVals$analysisOb@ExpSmpTab[[input$selectExpSmp]][1]),
                     reVals$factorsTab[[dx,1]],
                     choices=as.character(reVals$factorsTab[dx, 2:ncol(reVals$factorsTab)])) # overwrite the assignFactorsTab
        
      }
      }
      
     
      if (ncol(reVals$assignFactorsTab) > 1){
      # create new dataframe
      
      curExpSmpVals <- readDropDownDTtoTib(reVals$assignFactorsTab)
      
      # write name to analysis object (from the text input)
      if (input$nameExpSamTab == ""){ # check if the object has a name
        
        showModal(modalDialog(title = "Warning", "Please give your experiment sample table a name"))
        
      } else if (input$nameExpSamTab %in% names(reVals$analysisOb@ExpSmpFactorsTab)){
        showModal(modalDialog(title = "Warning", "Experiment sample table name already exists"))
      } 
            else{
        
        # write table to analysis object
        reVals$analysisOb <- addExpSmpTab(reVals$analysisOb, curExpSmpVals, input$nameExpSamTab, as.numeric(input$selectFacTab))
        
        # Update dropdown
        updateSelectInput(session, "selectExpSmp", choices = names(reVals$analysisOb@ExpSmpTab), selected = input$nameExpSamTab)
        
      }
      

    }
    }
    )
    
    # re render table with sample drop down
    observeEvent(input$selectExpSmp, {
      
      #browser()
      if(!(input$selectExpSmp == "")){
      reVals$assignFactorsTab <- readTibToDropDown(reVals$analysisOb@ExpSmpTab[[input$selectExpSmp]], 
                                                   reVals$analysisOb@factorsTab[[reVals$analysisOb@ExpSmpFactorsTab[[input$selectExpSmp]]]])
      }
      
      
    })
    
    # remove experimental sample
    observeEvent(input$rmExpSmpTab, {
      #browser()
      reVals$analysisOb <- rmExpSmpTab(reVals$analysisOb, input$selectExpSmp)
      
      if (length(reVals$analysisOb@ExpSmpTab) < 1){
        showModal(modalDialog(title = "Warning", "No I like this one, I don't want to delete it. (Ok so it makes the code so much easier if I don't have at least one table)"))
      }else{
      updateSelectInput(session, "selectExpSmp", choices = names(reVals$analysisOb@ExpSmpTab), selected = names(reVals$analysisOb@ExpSmpFactorsTab)[[1]])
      }
      
    })
    
    
## Normalization tab ============================================================================================== ###
     # Set the select experimental sample table, set to be same as the one on the factors page
   # browser()
   output$nrmCntsDT = DT::renderDataTable(reVals$nrmedCntsTab, server = FALSE, options = list(dom = 't'), rownames = F, class = 'cell-border stripe')
    observeEvent(input$selectExpSmp, {
      
      #browser()
      updateSelectInput(session, "selectExpSmpNrm", choices = names(reVals$analysisOb@ExpSmpTab), selected = input$selectExpSmp)
      
    })
    
    # When select Exp sample for normalization is selected update selectDesignFactors table
    observeEvent(input$selectExpSmpNrm, {
      #browser()
      if (!input$selectExpSmpNrm == ""){
      output$selectDesignFactors = DT::renderDataTable((reVals$analysisOb@factorsTab[[reVals$analysisOb@ExpSmpFactorsTab[[input$selectExpSmpNrm]]]])["Factors"], server = FALSE, options = list(dom = 't'))
     }
    })
    
    # Create a new normed counts deseq2dataset
    observeEvent(input$newDESeqNrmCnts, {
      #browser()
      
      if (input$nrmedCntsName == ""){showModal(modalDialog(title = "Warning", "Please give the normed counts DESeq dataset a name"))}
      else if(input$nrmedCntsName %in% names(reVals$analysisOb@NrmCnts)){
        showModal(modalDialog(title = "Warning", "Normed counts name already used"))
        
      }
      else{
        
        #browser()
        desFac <- (reVals$analysisOb@factorsTab[[reVals$analysisOb@ExpSmpFactorsTab[[input$selectExpSmpNrm]]]])["Factors"][input$selectDesignFactors_rows_selected,]
        if (isEmpty(desFac)){
          showModal(modalDialog(title = "Warning", "No Design factors selected"))
        }else {
          
          reVals$analysisOb <- deseq2CntNrm(reVals$analysisOb, input$selectExpSmpNrm, desFac, input$nrmedCntsName, input$rmLowCnts)
          
          #browser()
          
          if (reVals$nrmedCntsTab$Name[1] == ""){
            
            reVals$nrmedCntsTab <- tibble(Name = names(reVals$analysisOb@NrmCnts), `Exp Smp List` = reVals$analysisOb@NrmCntsExpSmp[[1]], `Norm Method` = reVals$analysisOb@NrmCnts[[1]]@NrmMethod,
                                          `Design Factors` = reVals$analysisOb@NrmCnts[[1]]@design)
            
          } else {
            
            reVals$nrmedCntsTab <- rbind(reVals$nrmedCntsTab , tibble(Name =last(names(reVals$analysisOb@NrmCnts)), `Exp Smp List` = last(reVals$analysisOb@NrmCntsExpSmp), 
                                                                      `Norm Method` = last(reVals$analysisOb@NrmCnts)@NrmMethod, `Design Factors` = last(reVals$analysisOb@NrmCnts)@design))
            
          }         
          
          #update normalised counts to select drop down
          updateSelectInput(session, "selectNrmCntsDiff", choices = names(reVals$analysisOb@NrmCnts), selected = names(reVals$analysisOb@NrmCnts)[[1]])
          #browser()
          # p = reVals$analysisOb@factorsTab[[reVals$analysisOb@ExpSmpFactorsTab[[reVals$analysisOb@NrmCnts[[1]]@ExpSampNm]]]]  # this is evil
          # # 
          # updateSelectInput(session, "selectNrmCntsCntrst",
          #                   choices = p$Factors,  # this is evil
          #                   selected = p$Factors[1])
          # 
          # 
          # 
          # updateSelectInput(session, "selectNrmCntsCnd",
          #                   choices = as.character(p[p$Factors == p$Factors[1],][2:ncol(p)]),
          #                   select = as.character(p[p$Factors == p$Factors[1],][2:ncol(p)])[1])
          # 
          # 
        }
      }

      
      #output$nrmCntsDT <- DT::renderDataTable(isolate(reVals$nrmedCntsTab), server = FALSE, options = list(dom = 't'))
      
      
    })
    
    ## Differential tab ============================================================================================== ###
     # update normalised counts to select drop down
    #updateSelectInput(session, "selectNrmCntsDiff", choices = names(reVals$analysisOb@nrmCnts))
    
    # need to update the contrast drop down if the nrmed counts change changes
    observeEvent(input$selectNrmCntsDiff, {
      reVals$analysisOb@State["nrmCnts", "selected"] <- input$selectNrmCntsDiff
      #browser()
      if (length(reVals$analysisOb@NrmCnts) > 0){
      p = unlist(strsplit(reVals$analysisOb@NrmCnts[[input$selectNrmCntsDiff]]@design, '\\+'))  # this is evil
      #p <- p[!(p == "+")]
      # 
      updateSelectInput(session, "selectNrmCntsCntrst",
                        choices = p,  # this is evil
                        selected = p[[1]])
      }
      
      if(length(reVals$analysisOb@NrmCntsExpSmp) > 1){
        tNrmCntsSel <- reVals$analysisOb@NrmCntsExpSmp[[input$selectNrmCntsDiff]] # get the expSmpTab in the currently selected nrmcounts tbale
      tFacTabSel <- reVals$analysisOb@ExpSmpFactorsTab[[tNrmCntsSel]] # get index of factors table
      tFacTab <- reVals$analysisOb@factorsTab[[tFacTabSel]] # get the factors tab
      tLev <- as.character(tFacTab[tFacTab[,1] == input$selectNrmCntsCntrst, ][2:ncol(tFacTab)])
      updateSelectInput(session, "selectNrmCntsCnd",
                        choices = tLev,
                        selected = tLev[[1]])
      }
      
    })
    
    # need to update the condition drop down if the contrast changes
    
    observeEvent(input$selectNrmCntsCntrst, {
      reVals$analysisOb@State["contrastFac", "selected"] <- input$selectNrmCntsCntrst

      #browser()

    if(length(reVals$analysisOb@NrmCntsExpSmp) > 1){
      tNrmCntsSel <- reVals$analysisOb@NrmCntsExpSmp[[input$selectNrmCntsDiff]] # get the expSmpTab in the currently selected nrmcounts tbale
    tFacTabSel <- reVals$analysisOb@ExpSmpFactorsTab[[tNrmCntsSel]] # get index of factors table
    tFacTab <- reVals$analysisOb@factorsTab[[tFacTabSel]] # get the factors tab
    tLev <- as.character(tFacTab[tFacTab[,1] == input$selectNrmCntsCntrst, ][2:ncol(tFacTab)])
    updateSelectInput(session, "selectNrmCntsCnd",
                      choices = tLev,
                      selected = tLev[[1]])
      }

    })
    
    observeEvent(input$selectNrmCntsCnd, {
      #browser()
      reVals$analysisOb@State["contrastCond", "selected"] <- input$selectNrmCntsCnd
      
    })
    
    # Run DESeq
    observeEvent(input$runDESeq2, {
      
      if(length(reVals$analysisOb@NrmCntsExpSmp) > 1){
        nt <- showNotification("DESeq2 Running", duration = NULL)
      reVals$analysisOb <- deseq2DA(reVals$analysisOb, input$selectNrmCntsDiff)
        removeNotification(nt)
      } else {
        
        showModal(modalDialog(title = "Warning", "No Normalized counts found"))
      }
    })
    
    # get DESeq results
    observeEvent(input$deseqResAct, 
                 {
                   #browser()
                   if (length(reVals$analysisOb@NrmCntsExpSmp) > 1) # 
                   {x <- deseq2Res(reVals$analysisOb, input$selectNrmCntsDiff, input$selectNrmCntsCntrst, input$selectNrmCntsCnd, input$LFCSelect)
                    if(any(class(x) == "error")){
                      showModal(modalDialog(title = "Warning", x$message))
                    }
                    else {
                      reVals$analysisOb <- x
                      browser()
                      # Make MA plot
                      
                      fig <- maplot(reVals$analysisOb@NrmCnts[[input$selectNrmCntsDiff]]@res)
                       output$MAPlot <-renderPlotly(fig)
                       
                       fig2 <- volplot(reVals$analysisOb@NrmCnts[[input$selectNrmCntsDiff]]@res)
                       output$volPlot <- renderPlotly(fig2)
                      
                    }
                   
                   }
                   else {
                     
                     showModal(modalDialog(title = "Warning", "No Normalized counts found"))
                   }
                  # browser()
                   
                   
                 })
    
    
    observeEvent(input$plotMAAct, {
      if (length(reVals$analysisOb@NrmCnts[[input$selectNrmCntsDiff]]@res) > 0 )
     { fig <- maplot(reVals$analysisOb@NrmCnts[[input$selectNrmCntsDiff]]@res)
      output$MAPlot <-renderPlotly(fig)}
      else
      {
        showModal(modalDialog(title = "Warning", "DESeq2 results not found, have you computed them yet?"))
        
      }
    })
    
    observeEvent(input$lfcplotMAAct, {
      if (length(reVals$analysisOb@NrmCnts[[input$selectNrmCntsDiff]]@resLFC) > 0 )
      { fig <- maplot(reVals$analysisOb@NrmCnts[[input$selectNrmCntsDiff]]@resLFC)
      output$MAPlot <-renderPlotly(fig)}
      else
      {
        showModal(modalDialog(title = "Warning", "DESeq2 LFC results not found, have you computed them?"))
        
      }
    })
    
### ================================================================================================================ ###
### Output Objects ================================================================================================= ###
    
    
    
 
}

# Run the application 
shinyApp(ui = ui, server = server)
