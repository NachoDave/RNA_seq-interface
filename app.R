#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(BiocManager)
options(repos = BiocManager::repositories())
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
source("geneTable.R")
source("pathWaySave.R")
options(shiny.maxRequestSize = 50*1024^2)
library(EnsDb.Hsapiens.v86)
library(WebGestaltR)
library(apeglm)
library(ashr)
#library()
# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # use css to define boxes around groups of buttons
  tags$head(tags$style(
    HTML('
         #DESeqRun, #DESeqGetRes {
            border: 1px solid black;
            background-color:LightCyan;
        }
        body, label, input, button, select { 
          font-family: "Arial";
        }')
  )),
  
    
    fluidRow(
        # Left menu shwing analysis and loaded files
        column(3, style = "background-color:lightblue;", 
              # fluidRow(h4("Data and Workspace")),
               
               sidebarLayout(
                 sidebarPanel(width = 12, id = "DESeqRun",
                              fluidRow(h3('Load and Save Workspace')),
                              #fluidRow(actionButton("newAnalysisactBut","Start New Analysis")),
                              fluidRow(fileInput(label = "Load WorkSpace", inputId = "loadAnalysis")),
                              fluidRow(h4('Download Workspace')),
                              # fluidRow(shinyDirButton("dir", "Select Save Directory", "Upload"),
                              #          verbatimTextOutput("dir", placeholder = TRUE)  
                              # ),
                              # 
                              fluidRow(textInput(label = "Workspace Filename", inputId = "saveFn")),
                              #fluidRow(actionButton("saveAnalysis", "Save Work space")),
                              downloadButton(label = "Save Workspace", outputId = "dwnLdWrkSpc"),
                              #fluidRow(tableOutput("cntTables"))
       
                              
                 ),
                 mainPanel(width = 0)
               ),
               
               
               sidebarLayout(
                 sidebarPanel(width = 12, id = "DESeqRun",
                              fluidRow(h3('Input Data')),
                              fluidRow(fileInput(label = "Load Gene Count table", inputId = "loadCntTab")),
                              
                              # table showing count data
                              fluidRow(h4("Loaded Gene Count tables")),
                              #fluidRow(DTOutput("genCntTbl")),
                              fluidRow(
                                
                                
                                DT::dataTableOutput('genCntTabTab'),
                                #verbatimTextOutput('y12')
                              ),
                              fluidRow(actionButton("removeGnCnt", "Remove Gene Counts")),
                              #fluidRow(h4('Factors tables'), DT::dataTableOutput('ldedGnFacsTab')),
                              
                 ),
                 mainPanel(width = 0)
               ),
               
    
              

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
                            #column(1),
                            #column(11,
                           # fluidRow(
                           sidebarLayout(
                             sidebarPanel(width = 12, id = "DESeqRun",
                                          fluidRow(h3('Create Factors Table')),
                                          fluidRow(
                                            actionButton(label = "New factors table", inputId = "newFacTab", icon = icon("plus")),
                                            actionButton(label = "Remove factors table", inputId = "rmFacTab", icon = icon("minus"))
                                            
                                          ), 
                                          fluidRow(
                                            h4('Factors (double click on table cell to enter factors and levels)'),
                                            DT::dataTableOutput('factors'),
                                          ),
                                          
                                          
                                          fluidRow(
                                            actionButton(label = "Add Factor", inputId = "addFactor", icon = icon("plus")),
                                            actionButton(label = "Add level", inputId = "addLevel", icon = icon("plus")),
                                            actionButton(label = "Remove Factor", inputId = "rmFactor", icon = icon("minus")),
                                            actionButton(label = "Remove level", inputId = "rmLevel", icon = icon("minus")),
                                          ),
                                          fluidRow(actionButton(label = "Add Factors to Workspace", inputId = "updateFactors")
                                          ),
                                          
                             ),
                             mainPanel(width = 0)
                           ),

                           sidebarLayout(
                             sidebarPanel(width = 12, id = "DESeqRun",
                                          fluidRow(h3('Create Experiment Sample Table'),
                                          selectInput(label = "1. Select Factors Table", inputId = "selectFacTab", choices = NULL),
                                          ),
                                           
                                          fluidRow(h5(strong('2. Select Samples you want to analyse from the loaded Gene Count tables (in left panel)'))),
                                          fluidRow(h5(strong('3. Create a new Experiment sample table'))),
                                          fluidRow(actionButton("newExpSmpTab", "New Sample Exp table")), # start a new assign factors table
                                          fluidRow(
                                            
                                            textInput("nameExpSamTab", "4. Give your experiment sample table a name", placeholder = "Experiment name"),
                                            
                                            DT::dataTableOutput('assignfactors'),
                                            
                                          ),
                                          
                                            
                                            fluidRow(h5(strong('5. Add Experiment sample to workspace'))),
                                          fluidRow(actionButton("saveExpSmpTab", "Add Sample table to workspace"),
                                            
                                            
                                          ),
                                          
                                          fluidRow(selectInput(inputId = "selectExpSmp", label = "Select Experiment sample table (view tables created)", choices = NULL),
                                                   actionButton("rmExpSmpTab", "Remove Sample table from workspace"))
                                          
                             ),
                             mainPanel(width = 0)
                           ),
                            
                           
                            #),
                   ),
                   ## Normalization tab panel ---------------------------------------------------------------------------## 
                   tabPanel("Normalization", #column(1),
                            #column(11,
                            
                            sidebarLayout(
                              sidebarPanel(width = 12, id = "DESeqRun",
                                 
                                           h3("Normalization using DESeq2 median of ratios (other methods can be added if required)"),
                                           fluidRow(selectInput(inputId = "selectExpSmpNrm", label = "1. Select Experiment sample table for normalization", choices = NULL),),
                                           fluidRow(h5(strong("2. Select Design Factors for differential Analysis")),),
                                           fluidRow(
                                             DT::dataTableOutput("selectDesignFactors"),
                                           ),
                                           fluidRow(
                                             
                                             numericInput(inputId = "rmLowCnts", label = "3. Remove genes with total counts less than:", value = 10),
                                           ),
                                           fluidRow(textInput("nrmedCntsName", "4. Give your normed counts a name"),),
                                           fluidRow(h5(strong("5. Add normed counts to the workspace")),),
                                           
                                           fluidRow(actionButton("newDESeqNrmCnts", "New DESeq2 Norm"),)  
                                                    
                                           ),         
                                           
                              
                              mainPanel(width = 0)
                            
                            ),
                            
                            sidebarLayout(
                              sidebarPanel(width = 12, id = "DESeqRun",
                                           
                                           fluidRow(h3("Normed Count Matrices"),
                                                    DT::dataTableOutput("nrmCntsDT")),
                                           fluidRow(actionButton("rmNrmCnts", "Remove normalized counts"),)
                                           
                              ),
                              mainPanel(width = 0)
                            ),
                            
                            #fluidRow(,
                            
                            
                            
                            
                            ),                        
                   
                   ## Data Exploration tab panel ---------------------------------------------------------------------------## 
                   tabPanel("Data Exploration", column(2,
                                                       sidebarLayout(
                                                         sidebarPanel(width = 12, id = "DESeqRun",
                                                                      fluidRow(selectInput(label = "Select Normalized Counts", inputId = "selectNrmCntsExp", choices = NULL)),
                                                                      fluidRow(radioButtons(inputId = "expDtTrn", label = "Data Transform", choices = c("None", "vst", "rlog"))),
                                                                      fluidRow(actionButton(label = "PCA",inputId =  "pca")),
                                                                      fluidRow(actionButton(label = "Sample Distance",inputId =  "smpDist")),
                                                                      fluidRow(actionButton(label = "Heat Map",inputId =  "htMp")),
                                                                      
                                                         ),
                                                         mainPanel(width = 0)
                                                         
                                                         
                                                       ),
                   ),
                   column(5,
                          h3("PCA"),
                          fluidRow(plotlyOutput("PCAPlt")),
                          h3("Sample Distance"),
                          fluidRow(plotOutput("smpDstPlt")),
                          
                   ),
                   
                   column(5, 
                          
                          h3("Heat Map of Count Matrix"),
                          fluidRow(plotOutput("GnCntPlt")),
                          
                   )
                   ),
                   
                   ## Differential Analysis tab panel ---------------------------------------------------------------------------## 
                   tabPanel("Differential Analysis", # Check overlaps
                   
                   column(3,
                          sidebarLayout(
                              sidebarPanel(width = 12, id = "DESeqRun",
                                fluidRow(selectInput(label = "Select Normalized Counts", inputId = "selectNrmCntsDiff", choices = NULL), 
                                   selectInput(label = "Select Contrast Factor", inputId = "selectNrmCntsCntrst", choices = NULL),
                                   selectInput(label = "Select constrast condition", inputId = "selectNrmCntsCnd", choices = NULL)),
                                fluidRow(actionButton(label = "Run DESeq2", inputId = "runDESeq2")),
                                
                              ),
                                      mainPanel(width = 0)
                                      
                                
                          ),
                          
                          
                          sidebarLayout(
                            sidebarPanel(width = 12, id = "DESeqGetRes",
                              fluidRow(actionButton(label = "Get DESeq2 results", inputId = "deseqResAct"), 
                                       checkboxInput(label = "Calculate LFC", inputId = "LFCSelect"),
                                       selectInput(label = "Select LFC Method", inputId = "LFCMethod", choices = c("apeglm", "ashr", "normal"))
                                       ),
                          
                            ),
                            mainPanel(width = 0)
                          ),
                          
                          
                          sidebarLayout(
                              sidebarPanel(width = 12, id = "DESeqGetRes", 
                                fluidRow(actionButton(label = "MA Plot", inputId = "plotMAAct")), 
                                fluidRow(actionButton(label = "Volcano Plot", inputId = "plotVolcanoAct")),
                                checkboxInput(label = "Plot LFC Result?", inputId = "plotLFCChk")
                              ), 
                              mainPanel(width = 0)
                          )
                   ),

                          column(7, fluidRow(column(11, offset = 1, h3("Differential Analysis MA Plot"))),
                          (fluidRow(column(9,plotlyOutput("MAPlot")), column(2,
                                  numericInput(label = "Lower Limit Y", inputId = "MAlwLimY", value = ""), 
                                   numericInput(label = "Upper Limit Y", inputId = "MAupLimY", value = ""), 
                                   numericInput(label = "Lower Limit X", inputId = "MAlwLimX", value = ""), 
                                   numericInput(label = "Upper Limit X", inputId = "MAupLimX", value = ""),
                                   actionButton(label = "Replot", inputId = "replotMA")))),
                          fluidRow(column(11, offset = 1, h3("Differential Analysis Volcano Plot"))),
                          (fluidRow(column(9,plotlyOutput("volPlot")), column(2,
                                                                              numericInput(label = "Lower LimitY", inputId = "vollwLimY", value = ""), 
                                                                             numericInput(label = "Upper LimitY", inputId = "volupLimY", value = ""), 
                                                                             numericInput(label = "Lower Limit X", inputId = "vollwLimX", value = ""), 
                                                                             numericInput(label = "Upper Limit X", inputId = "volupLimX", value = ""),
                                                                             actionButton(label = "Replot", inputId = "replotVol")))),
                          ),
                          
                          
                   
                   ),
                   tabPanel("Genelists/Sequences",
                      column(3,
                        
                             
                             
                        # Create New gene list panel
                       sidebarLayout(
                        sidebarPanel(width = 12, id = "DESeqRun",
                                     fluidRow(h3("Create Gene list from Differential Analysis")),
                          fluidRow(selectInput(label = "1. Select Differential Analysis", inputId = "selectDiffAnl", choices = NULL)),
                          fluidRow(checkboxInput(label = "Use LFC", inputId = "LFCGeneList")),
                          fluidRow(radioButtons(inputId = "convertGeneIDtoSym", label = "2. Convert Gene IDs to symbols?", choices = c("No", "Ensembl", "Entrez"))),
                          fluidRow(h5(strong("3. Create New Gene List"))),
                          fluidRow(actionButton(label = "Create Gene List",inputId =  "createGeneList")), 
                          fluidRow(h5(strong("4. Filter Gene List (using the table filter)"))),
                                   fluidRow(textInput(inputId = "nameGenelist", label = "5. Name Gene List"), ),
                                  fluidRow(h5(strong("6. Save Gene list to file"))),
                                   fluidRow(actionButton(label = "Save Gene List", "saveGeneList"), 
                                   ),
                          #fluidRow(actionButton(label = "Write Gene List to file", inputId = "writeGeneList"))
                         ),
                                     mainPanel(width = 0)
                       ),
                       
                       # Select Gene List panel
                       sidebarLayout(
                         sidebarPanel(width = 12, id = "DESeqRun",
                                      selectInput(label = "View Gene Lists in Workspace", inputId = "selectGeneTableGnTab", choices = NULL)
                                      
                         ),
                         mainPanel(width = 0)
                       ),
                       
                       # Select Gene List panel
                       sidebarLayout(
                         sidebarPanel(width = 12, id = "DESeqRun",
                                      fluidRow(h3("Download Current Gene list to csv")),
                                      fluidRow(downloadButton(label = "Download Gene List", outputId = "writeGeneList")),
                                      
                         ),
                         mainPanel(width = 0)
                       ),
                       
                       
                       
                       
                      ),
                      column(9,
                             sidebarLayout(
                               sidebarPanel(width = 12, id = "DESeqRun",
                                fluidRow(h3("Gene List"),
                                DT::dataTableOutput("geneListDT")),
                                helpText("To filter the based on values type string 'lower Value ... Upper Value' into 
                                         filter boxes after table, i.e to filter p value between 0 and 0.05 type 0 ... 0.5 into the p value filter")
                              ),
                              mainPanel(width = 0)
                               
                             )
                             
                       )
      
                   ),
                   
                   
                   tabPanel("Pathway Analysis", 
                            column(2,
                                   # Select Gene List panel
                                   sidebarLayout(
                                     sidebarPanel(width = 12, id = "DESeqRun",
                                                  fluidRow(h3("Pathway Analysis (you can perform as many as you like)")),
                                                  selectInput(label = "1. Select Gene List", inputId = "selectGOGnTb", choices = NULL),
                                                  
                                   #   ),
                                   #   mainPanel(width = 0)
                                   # ),
                                   
                                   # Select GO method
                                   # sidebarLayout(
                                   #   sidebarPanel(width = 12, id = "DESeqRun",
                                        
                                                  fluidRow(checkboxGroupInput(label = "2. Select GO Method(s)", inputId = "selectGOMethod", choices = c("GSEA", "ORA"))),
                                                  fluidRow(selectInput(label = "3. Select Enrichment Metric for GSEA", inputId = "selectGSEAEnrich", choices = NULL)),
                                                  
                                   #   ),
                                   #   mainPanel(width = 0)
                                   # ),
                                   

                                   # Select GO process
                                   # sidebarLayout(
                                   #   sidebarPanel(width = 12, id = "DESeqRun",
                                                  
                                                  checkboxGroupInput(label = "4. Select Analysis", inputId = "selectGODB", 
                                                                     
                                                                     choices = c("GO Molecular Function", "GO Cellular Component", "GO Biological Process", # GO 
                                                                                 "KEGG Pathway", "Panther Pathway", 
                                                                                 "Reactome Pathway")), # pathway
                                                  
                                   #   ),
                                   #   mainPanel(width = 0)
                                   # ),
                                   
                                   # Select GO database
                                   # sidebarLayout(
                                   #   sidebarPanel(width = 12, id = "DESeqRun",
                                   #                
                                   #                checkboxGroupInput(label = "Select GO Database(s)", inputId = "selectGODB", choices = c("GO", "KEGG", "PANTHER")),
                                   #                
                                   #   ),
                                   #   mainPanel(width = 0)
                                   # ),
                                   
                                   # Run GO
                                   # sidebarLayout(
                                   #   sidebarPanel(width = 12, id = "DESeqRun",
                                                  
                                                  # checkboxInput(inputId = "saveGOChk", label = "Save Analysis Results"),
                                                  # fluidRow(shinyDirButton("dirGO", "Save Directory", "Upload"),
                                                  #          verbatimTextOutput("dirGO", placeholder = TRUE)  
                                                  # ),
                                                  # textInput(inputId = "textSaveGO", label = "Name GO Analysis"),
                                                  fluidRow(h5(strong("5. Run Pathway Analysis"))), 
                                                  actionButton(inputId = "runGO", label = "Run Analysis"),
                                                  
                                      ),
                                     mainPanel(width = 0)
                                   ),
                            ),
                            
                            column(7, fluidRow(h3("Pathway Analysis"),
                                               plotlyOutput("pthWyPlt")),
                                   fluidRow( uiOutput("pthWyUrl"), #textOutput("pthWyTxt"), 
                                             DT::dataTableOutput("pthWyTb")),),
                            column(3,  sidebarLayout(
                              sidebarPanel(width = 12, id = "DESeqRun",
                                           selectInput(label = "Select Pathway Analysis", inputId = "selectPthWy", choices = NULL),
                                           actionButton(inputId = "pthWyBarPlot", label = "Bar Plot")
                                           
                              ),
                              mainPanel(width = 0)
                            ),
                                   
                                   
                                   
                            )
                            
                   ),
                   
                   
                   
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
                             nrmedCntsTab = tibble(Name = "", `Exp Smp List` = "", `Norm Method` = "", `Design Factors` = ""),
                             maFig = ggplotly(), volFig = ggplotly(), 
                             curGeneTab = data.frame()
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
    
    # Save analysis = for saving to server
    # shinyDirChoose(
    # 
    #     input,
    #     'dir',
    #     roots = c(home = "/data/RNASeqAnalysis"),
    #     #filetypes = c('', 'txt', 'bigWig', "tsv", "csv", "bw")
    # )
    # #browser()
    # global <- reactiveValues(datapath = getwd(), datapathGO = getwd())
    # 
    # dir <- reactive(input$dir)
    # 
    # output$dir <- renderText({
    #     "Select Save Path"
    # })
    # 
    # observeEvent(ignoreNULL = TRUE,
    #              eventExpr = {
    #                  input$dir
    #              },
    #              handlerExpr = {
    #                  if (!"path" %in% names(dir())) return()
    #                  home <- normalizePath("/data/RNASeqAnalysis")
    #                  global$datapath <-
    #                      file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
    #                  #print(global$datapath)
    #                  output$dir <- renderText({
    #                      global$datapath
    #                  })
    #              })

    # Save the analysis
    # observeEvent(input$saveAnalysis, {
    #     #browser()
    #     # saveReVals <- reVals$analysisOb
    #     # 
    #     # if (input$saveFn == ""){save(saveReVals, file = paste0(global$datapath, "/", "RNASeqAnalysis.rData"))}
    #     # else {save(saveReVals, file = paste0(global$datapath, "/",input$saveFn, ".rData"))}
    # 
    # })
    
    output$dwnLdWrkSpc <- downloadHandler(
      
      filename = function() {
        if (input$saveFn == ""){
          fn <- "RNASeqAnalysis"
        }
          else
          {
            fn <- input$saveFn
          }
        
        paste(fn, ".rData", sep = "")
      },
      content = function(file) {
        #browser()
        saveReVals <- reVals$analysisOb
        save(saveReVals, file = file)
      }
    )
    
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
                       updateSelectInput(session, "selectNrmCntsExp", choices = names(reVals$analysisOb@NrmCnts), selected = names(reVals$analysisOb@NrmCnts)[[1]])
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
                     
                     # Populate the Gene table 
                     if (length(reVals$analysisOb@GeneTables) >0){
                       #browser()
                       updateSelectInput(session, "selectGeneTableGnTab", choices = names(reVals$analysisOb@GeneTables), selected =  names(reVals$analysisOb@GeneTables)[[1]])
                       updateSelectInput(session, "selectGOGnTb", choices = names(reVals$analysisOb@GeneTables), selected =  names(reVals$analysisOb@GeneTables)[[1]]) # update gene table on GO tab
                       reVals$curGeneTab <- reVals$analysisOb@GeneTables[[1]]@gnTbl
                       gnTbDT <- reVals$curGeneTab
                       gnTbDT[sapply(reVals$curGeneTab, is.numeric)] <- round(reVals$curGeneTab[sapply(reVals$curGeneTab, is.numeric)], 5)

                       output$geneListDT <- DT::renderDataTable(gnTbDT
                                                                , server = TRUE,  rownames = T)


                     }
                     #browser()
                     if (length(reVals$analysisOb@PthWyAnl) > 0 )
                     {
                        # browser()
                       tMet <- names(reVals$analysisOb@PthWyAnl)
                     updateSelectInput(session, inputId = "selectPthWy", 
                                       choices = tMet)}
                       
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
          updateSelectInput(session, "selectNrmCntsExp", choices = names(reVals$analysisOb@NrmCnts), selected = names(reVals$analysisOb@NrmCnts)[[1]])
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
    
    ## Data Exploration tab ========================================================================================== ###
    # plot distance heat map
    observeEvent(input$smpDist, {
      #browser()
      nt <- showNotification("Making Distance heatmap", duration = NULL)
      fig <- pltSmpDist(reVals$analysisOb@NrmCnts[[input$selectNrmCntsExp]], input$expDtTrn)
      output$smpDstPlt <- renderPlot(fig)
      removeNotification(nt)
    })
    
    # plot count heat map
    observeEvent(input$htMp, {
      nt <- showNotification("Making Count Matrix Heatmap", duration = NULL)
      fig <- pltHtMp(reVals$analysisOb@NrmCnts[[input$selectNrmCntsExp]], input$expDtTrn)
      output$GnCntPlt <- renderPlot(fig)
      removeNotification(nt)
    })
    
    # PCA plot
    observeEvent(input$pca, {
      nt <- showNotification("Making PCA plot", duration = NULL)
      fig <- pltPCA(reVals$analysisOb@NrmCnts[[input$selectNrmCntsExp]], input$expDtTrn)
      output$PCAPlt <- renderPlotly(fig)
      removeNotification(nt)      
      
      
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
      
      if(length(reVals$analysisOb@NrmCntsExpSmp) > 0){
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

    if(length(reVals$analysisOb@NrmCntsExpSmp) > 0){
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
    
      if(length(reVals$analysisOb@NrmCntsExpSmp) > 0){
        #browser()
        nt <- showNotification("DESeq2 Running", duration = NULL)
      reVals$analysisOb <- deseq2DA(reVals$analysisOb, input$selectNrmCntsDiff, input$selectNrmCntsCntrst, input$selectNrmCntsCnd)
        removeNotification(nt)
      } else {
        
        showModal(modalDialog(title = "Warning", "No Normalized counts found"))
      }
    })
    
    # get DESeq results
    observeEvent(input$deseqResAct, 
                 {
                   #browser()
                   if (length(reVals$analysisOb@NrmCntsExpSmp) > 0) # 
                   {x <- deseq2Res(reVals$analysisOb, input$selectNrmCntsDiff, input$selectNrmCntsCntrst, input$selectNrmCntsCnd, input$LFCSelect, input$LFCMethod)
                    if(any(class(x) == "error")){
                      showModal(modalDialog(title = "Warning", x$message))
                    }
                    else {
                      nt <- showNotification("Getting DESeq2 Results", duration = NULL)
                      reVals$analysisOb <- x
                      #browser()
                      # Make MA plot
                      
                      reVals$maFig <- maplot(reVals$analysisOb@NrmCnts[[input$selectNrmCntsDiff]]@res)
                       output$MAPlot <-renderPlotly(reVals$maFig)
                       
                       updateNumericInput(session, "MAlwLimX", value = reVals$maFig$x$layout$xaxis$range[1])
                       updateNumericInput(session, "MAupLimX", value = reVals$maFig$x$layout$xaxis$range[2])
                       
                       updateNumericInput(session, "MAlwLimY", value = reVals$maFig$x$layout$yaxis$range[1])
                       updateNumericInput(session, "MAupLimY", value = reVals$maFig$x$layout$yaxis$range[2])
                       
                       #browser()
                       reVals$volFig <- volplot(reVals$analysisOb@NrmCnts[[input$selectNrmCntsDiff]]@res)
                       output$volPlot <- renderPlotly(reVals$volFig)
                       
                       updateNumericInput(session, "vollwLimX", value = reVals$volFig$x$layout$xaxis$range[1])
                       updateNumericInput(session, "volupLimX", value = reVals$volFig$x$layout$xaxis$range[2])
                       
                       updateNumericInput(session, "vollwLimY", value = reVals$volFig$x$layout$yaxis$range[1])
                       updateNumericInput(session, "volupLimY", value = reVals$volFig$x$layout$yaxis$range[2])                       
                       
                       removeNotification(nt)

                    }
                   
                   }
                   else {
                     
                     showModal(modalDialog(title = "Warning", "No Normalized counts found"))
                   }
                  # browser()
                   
                   
                 })
    
    
    observeEvent(input$plotMAAct, {
      if (length(reVals$analysisOb@NrmCnts[[input$selectNrmCntsDiff]]@res) > 0 )
     { 
        
        if (input$plotLFCChk) # check if plot the lfc checkbox is clicked
        {
          if (length(reVals$analysisOb@NrmCnts[[input$selectNrmCntsDiff]]@resLFC) > 0 ) # check if the lfc results have been computed
          {
            reVals$maFig <- maplot(reVals$analysisOb@NrmCnts[[input$selectNrmCntsDiff]]@resLFC, tit = " (LFC)")
            output$MAPlot <-renderPlotly(reVals$maFig)
            
            updateNumericInput(session, "MAlwLimX", value = reVals$maFig$x$layout$xaxis$range[1])
            updateNumericInput(session, "MAupLimX", value = reVals$maFig$x$layout$xaxis$range[2])
            
            updateNumericInput(session, "MAlwLimY", value = reVals$maFig$x$layout$yaxis$range[1])
            updateNumericInput(session, "MAupLimY", value = reVals$maFig$x$layout$yaxis$range[2])
            
          } 
          else
          {
            showModal(modalDialog(title = "Warning", "DESeq2 LFC results not found, have you computed them?"))
          }
        } 
        else
        {
          reVals$maFig <- maplot(reVals$analysisOb@NrmCnts[[input$selectNrmCntsDiff]]@res)
          output$MAPlot <-renderPlotly(reVals$maFig)
          
          updateNumericInput(session, "MAlwLimX", value = reVals$maFig$x$layout$xaxis$range[1])
          updateNumericInput(session, "MAupLimX", value = reVals$maFig$x$layout$xaxis$range[2])
          
          updateNumericInput(session, "MAlwLimY", value = reVals$maFig$x$layout$yaxis$range[1])
          updateNumericInput(session, "MAupLimY", value = reVals$maFig$x$layout$yaxis$range[2])

        }
            
      }
      else
      {
        showModal(modalDialog(title = "Warning", "DESeq2 results not found, have you computed them yet?"))
        
      }
    })
    
    # observeEvent(input$lfcplotMAAct, {
    #   if (length(reVals$analysisOb@NrmCnts[[input$selectNrmCntsDiff]]@resLFC) > 0 )
    #   { fig <- maplot(reVals$analysisOb@NrmCnts[[input$selectNrmCntsDiff]]@resLFC, tit = " (LFC)")
    #   output$MAPlot <-renderPlotly(fig)}
    #   else
    #   {
    #     showModal(modalDialog(title = "Warning", "DESeq2 LFC results not found, have you computed them?"))
    #     
    #   }
    # })
    
    observeEvent(input$plotVolcanoAct, {
      if (length(reVals$analysisOb@NrmCnts[[input$selectNrmCntsDiff]]@res) > 0 )
      { 
        
        if (input$plotLFCChk) # check if plot the lfc checkbox is clicked
        {
          if (length(reVals$analysisOb@NrmCnts[[input$selectNrmCntsDiff]]@resLFC) > 0 ) # check if the lfc results have been computed
          {
            reVals$volFig <- volplot(reVals$analysisOb@NrmCnts[[input$selectNrmCntsDiff]]@resLFC, tit = " (LFC)")
            output$volPlot <-renderPlotly(reVals$volFig)
            
            updateNumericInput(session, "vollwLimX", value = reVals$volFig$x$layout$xaxis$range[1])
            updateNumericInput(session, "volupLimX", value = reVals$volFig$x$layout$xaxis$range[2])
            
            updateNumericInput(session, "vollwLimY", value = reVals$volFig$x$layout$yaxis$range[1])
            updateNumericInput(session, "volupLimY", value = reVals$volFig$x$layout$yaxis$range[2]) 
          } 
          else
          {
            showModal(modalDialog(title = "Warning", "DESeq2 LFC results not found, have you computed them?"))
          }
        } 
        else
        {
          reVals$volFig <- volplot(reVals$analysisOb@NrmCnts[[input$selectNrmCntsDiff]]@res)
          output$volPlot <-renderPlotly(reVals$volFig)

          updateNumericInput(session, "vollwLimX", value = reVals$volFig$x$layout$xaxis$range[1])
          updateNumericInput(session, "volupLimX", value = reVals$volFig$x$layout$xaxis$range[2])
          
          updateNumericInput(session, "vollwLimY", value = reVals$volFig$x$layout$yaxis$range[1])
          updateNumericInput(session, "volupLimY", value = reVals$volFig$x$layout$yaxis$range[2]) 
          
        }
        
      }
      else
      {
        showModal(modalDialog(title = "Warning", "DESeq2 results not found, have you computed them yet?"))
        
      }
    })
    
    # observeEvent(input$lfcplotVolcanoAct, {
    #   if (length(reVals$analysisOb@NrmCnts[[input$selectNrmCntsDiff]]@resLFC) > 0 )
    #   { fig <- volplot(reVals$analysisOb@NrmCnts[[input$selectNrmCntsDiff]]@resLFC, tit = " (LFC)")
    #   output$volPlot <-renderPlotly(fig)}
    #   else
    #   {
    #     showModal(modalDialog(title = "Warning", "DESeq2 LFC results not found, have you computed them?"))
    #     
    #   }
    # })
    
    # Replot MA with limits
    observeEvent(input$replotMA, {
      if (input$plotLFCChk)
      {
        reVals$maFig <- maplot(reVals$analysisOb@NrmCnts[[input$selectNrmCntsDiff]]@resLFC, 
                               xlim = c(input$MAlwLimX, input$MAupLimX), 
                               ylims = c(input$MAlwLimY, input$MAupLimY), tit = " (LFC)") # need to say if its LFC or not??
      }
      else
      {
        reVals$maFig <- maplot(reVals$analysisOb@NrmCnts[[input$selectNrmCntsDiff]]@res, 
                             xlim = c(input$MAlwLimX, input$MAupLimX), 
                             ylims = c(input$MAlwLimY, input$MAupLimY)) # need to say if its LFC or not??
      }
      
      output$MAPlot <-renderPlotly(reVals$maFig)
      
      updateNumericInput(session, "MAlwLimX", value = reVals$maFig$x$layout$xaxis$range[1])
      updateNumericInput(session, "MAupLimX", value = reVals$maFig$x$layout$xaxis$range[2])
      
      updateNumericInput(session, "MAlwLimY", value = reVals$maFig$x$layout$yaxis$range[1])
      updateNumericInput(session, "MAupLimY", value = reVals$maFig$x$layout$yaxis$range[2])
    }
    )
    
    
    observeEvent(input$replotVol, {
      #browser()
      if (input$plotLFCChk)
      {
        reVals$volFig <- volplot(reVals$analysisOb@NrmCnts[[input$selectNrmCntsDiff]]@resLFC, 
                               xlim = c(input$vollwLimX, input$volupLimX), 
                               ylims = c(input$vollwLimY, input$volupLimY), tit = " (LFC)") # need to say if its LFC or not??
      }
      else
      {
        reVals$volFig <- volplot(reVals$analysisOb@NrmCnts[[input$selectNrmCntsDiff]]@res, 
                               xlim = c(input$vollwLimX, input$volupLimX), 
                               ylims = c(input$vollwLimY, input$volupLimY)) # need to say if its LFC or not??
      }
      
      output$volPlot <-renderPlotly(reVals$volFig)
      
      updateNumericInput(session, "vollwLimX", value = reVals$volFig$x$layout$xaxis$range[1])
      updateNumericInput(session, "volupLimX", value = reVals$volFig$x$layout$xaxis$range[2])
      
      updateNumericInput(session, "vollwLimY", value = reVals$volFig$x$layout$yaxis$range[1])
      updateNumericInput(session, "volupLimY", value = reVals$volFig$x$layout$yaxis$range[2])
    }
    )
    
### Gene list tab ================================================================================================== ###
### ================================================================================================================ ###
### Output Objects ================================================================================================= ###
    
    # Update the select diff analysis drop down
    observeEvent(reVals$analysisOb@NrmCnts, 
                 {
                   updateSelectInput(session, "selectDiffAnl", choices = names(reVals$analysisOb@NrmCnts))
                   
                 })
    
    # create new genelist
    observeEvent(input$createGeneList, {
      # I apoligze for the state of this code!
      
      # Check if using LFC result and there are LFC results present
      if (input$LFCGeneList && length(reVals$analysisOb@NrmCnts[[input$selectDiffAnl]]@resLFC) > 0){
        #browser()
        
        if (input$convertGeneIDtoSym == "Ensembl"){ 
          reVals$curGeneTab <- as.data.frame(reVals$analysisOb@NrmCnts[[input$selectDiffAnl]]@resLFC)
          # Get the gene IDs
          
          ensId <- gsub("\\.\\d+", "",rownames(reVals$curGeneTab))
          rownames(reVals$curGeneTab) <- ensId
          gnSym <- ensembldb::select(EnsDb.Hsapiens.v86, keys = ensId, keytype = "GENEID", columns = c("SYMBOL"))
          reVals$curGeneTab$`Gene Symbol` <- ""
          reVals$curGeneTab <- reVals$curGeneTab[ , c(ncol(reVals$curGeneTab), 1:ncol(reVals$curGeneTab) - 1)]
          reVals$curGeneTab$`Gene Symbol`[ensId %in% gnSym$GENEID] <- gnSym$SYMBOL
          
          # Make a datatable object for rendering and reducing column width
          # gnTbDT <- formatRound(datatable(reVals$curGeneTab), digits = 5, 
          #                       columns = colnames(reVals$curGeneTab)[sapply(reVals$curGeneTab, is.numeric)])

        }
        else
        {
          reVals$curGeneTab <- as.data.frame(reVals$analysisOb@NrmCnts[[input$selectDiffAnl]]@resLFC)
          # gnTbDT <- formatRound(datatable(reVals$curGeneTab), digits = 5, 
          #                       columns = colnames(reVals$curGeneTab)[sapply(reVals$curGeneTab, is.numeric)])
        }
      }
      else if (length(reVals$analysisOb@NrmCnts[[input$selectDiffAnl]]@res) > 0) {# don't use LFC results'
        
        if (input$convertGeneIDtoSym == "Ensembl"){ # convert ensembl id to gene symbol
          reVals$curGeneTab <- as.data.frame(reVals$analysisOb@NrmCnts[[input$selectDiffAnl]]@res)
          # Get the gene IDs
          
          ensId <- gsub("\\.\\d+", "",rownames(reVals$curGeneTab))
          rownames(reVals$curGeneTab) <- ensId
          gnSym <- ensembldb::select(EnsDb.Hsapiens.v86, keys = ensId, keytype = "GENEID", columns = c("SYMBOL"))
          reVals$curGeneTab$`Gene Symbol` <- ""
          reVals$curGeneTab <- reVals$curGeneTab[ , c(ncol(reVals$curGeneTab), 1:ncol(reVals$curGeneTab) - 1)]
          reVals$curGeneTab$`Gene Symbol`[ensId %in% gnSym$GENEID] <- gnSym$SYMBOL

          # Make a datatable object for rendering and reducing column width
          # gnTbDT <- formatRound(datatable(reVals$curGeneTab), digits = 5, 
          #                       columns = colnames(reVals$curGeneTab)[sapply(reVals$curGeneTab, is.numeric)])
          
        }
        else
        {
          reVals$curGeneTab <- as.data.frame(reVals$analysisOb@NrmCnts[[input$selectDiffAnl]]@res)
          
        }
        
      }
      else {
        

        showModal(modalDialog(title = "warning", "No DESeq results found. Have you run 'Get DESeq2 Results' in the differential Analysis tab?"))
        reVals$curGeneTab <- data.frame("")
        output$geneListDT <- DT::renderDataTable(reVals$curGeneTab
                                                 , server = TRUE,  rownames = T, filter = 'top')
        return()
      }
      #reVals$curGeneTab$Abs
      
      #gnTbDT <- #formatRound(datatable(reVals$curGeneTab), digits = 5, 
                 #           columns = colnames(reVals$curGeneTab)[sapply(reVals$curGeneTab, is.numeric)])      
      
      # absFc <- abs(reVals$curGeneTab$log2FoldChange)
      # reVals$curGeneTab$'abs log2FoldChange' <- absFc
      #browser()
      #reVals$curGeneTab <- reVals$curGeneTab[, c(1, 2,3, 5, 6)]
      gnTbDT <- reVals$curGeneTab
      gnTbDT[sapply(reVals$curGeneTab, is.numeric)] <- round(reVals$curGeneTab[sapply(reVals$curGeneTab, is.numeric)], 5)
      
      output$geneListDT <- DT::renderDataTable(gnTbDT
                                               , server = TRUE,  rownames = T, filter = 'top')
    })
    
    
    
    # Save a genelist to the workspace
    
    observeEvent(input$saveGeneList, {
      
      #browser()
       # this is the command reVals$curGeneTab[input$geneListDT_rows_all,]
     #tDt
      # Check there's a table there, that there is a name and the name hasn't been used before
      if (ncol(reVals$curGeneTab) > 1 && !(input$nameGenelist == "" && !(input$nameGenelist %in% names(reVals$analysisOb@GeneTables)))){
        
        tGnTbOb <- new("geneTableDA")
        tGnTbOb <- addGeneTableDA(tGnTbOb, reVals$curGeneTab[input$geneListDT_rows_all,], input$selectDiffAnl, input$geneListDT_search_columns)
        reVals$analysisOb <- addGeneTable(reVals$analysisOb, tGnTbOb, input$nameGenelist)
        #browser()
        updateSelectInput(session, "selectGeneTableGnTab", choices = names(reVals$analysisOb@GeneTables), selected =  input$nameGenelist)
        updateSelectInput(session, "selectGOGnTb", choices = names(reVals$analysisOb@GeneTables), selected =  names(reVals$analysisOb@GeneTables)[[1]]) # update gene table on GO tab
      } else if(ncol(reVals$curGeneTab) < 1)
      {
        
        showModal(modalDialog(title = "Warning", "No Gene Table"))
        
      } else if (input$nameGenelist == "") {
        
        showModal(modalDialog(title = "Warning", "Please Give the gene list a name"))
        
      } else if (input$nameGenelist %in% names(reVals$analysisOb@GeneTables)) {
        
        showModal(modalDialog(title = "Warning", "Gene List Name has already been used")) 
        
      }
      
    })
    
    # Select Gene List using drop down
    observeEvent(input$selectGeneTableGnTab, {
      if (length(reVals$curGeneTab) > 0)
      {
        reVals$curGeneTab <- reVals$analysisOb@GeneTables[[input$selectGeneTableGnTab]]@gnTbl
      gnTbDT <- reVals$curGeneTab
      gnTbDT[sapply(reVals$curGeneTab, is.numeric)] <- round(reVals$curGeneTab[sapply(reVals$curGeneTab, is.numeric)], 5)
      
      output$geneListDT <- DT::renderDataTable(gnTbDT, server = TRUE,  rownames = T)
      }
      
    })

       output$writeGeneList <- downloadHandler(
        
         filename = function() {
           paste(input$selectGeneTableGnTab, ".csv", sep = "")
         },
         content = function(file) {
           write.csv(reVals$curGeneTab, file, row.names = TRUE)
         }
       )

    ###  GO Analysis tab =============================================================================================== ###
       
       # Save GO analysis ---------------------------------------------------------------------------------------------- ###
       shinyDirChoose(
         input,
         'dirGO',
         roots = c(home = "/data/RNASeqAnalysis"),
         #filetypes = c('', 'txt', 'bigWig', "tsv", "csv", "bw")
       )
       
       # global <- reactiveValues(datapathGO = getwd())
       # 
       dirGO <- reactive(input$dirGO)

       output$dirGO <- renderText({
         "Select Save Path"
       })

       observeEvent(ignoreNULL = TRUE,
                    eventExpr = {
                      input$dirGO
                    },
                    handlerExpr = {
                      if (!"path" %in% names(dirGO())) return()
                      home <- normalizePath("/data/RNASeqAnalysis")
                      global$datapathGO <-
                        file.path(home, paste(unlist(dirGO()$path[-1]), collapse = .Platform$file.sep))
                      #print(global$datapath)
                      output$dirGO <- renderText({
                        global$datapathGO
                      })
                    })
       # 
       # Run GO
       observeEvent(input$runGO, {
         
         
         if (is.null(input$selectGODB)){
           
           showModal(modalDialog(title = "Warning", "No Database selected"))
          return()  
         }
         if (is.null(input$selectGOMethod)){
           
           showModal(modalDialog(title = "Warning", "No Analysis Method selected"))
           return()  
         }         
         # if (input$textSaveGO == ""){
         #   
         #   showModal(modalDialog(title = "Warning", "Please name your GO analysis"))
         #   return() 
         #   
         # }
         # if (input$dirGO == "" ){
         #   
         #   showModal(modalDialog(title = "Warning", "Please set path for output directory"))
         #   return()            
         #   
         # }
         # 
         # Find ID type from gene table
         
         nt <- showNotification("Pathway Analysis Running", duration = NULL)
         gns <- data.frame(rownames(reVals$analysisOb@GeneTables[[input$selectGOGnTb]]@gnTbl),
                  reVals$analysisOb@GeneTables[[input$selectGeneTableGnTab]]@gnTbl[[input$selectGSEAEnrich]])
         colnames(gns) <- c("ID", "input$selectGSEAEnrich")
         
         # x <- runWebGestaltR(reVals$analysisOb, input$selectGODB, input$selectGOMethod,
         #                     global$datapathGO, input$textSaveGO, gns, "ensembl_gene_id", input$saveGOChk, 
         #                     input$selectGOGnTb)
         x <- runWebGestaltR(reVals$analysisOb, input$selectGODB, input$selectGOMethod,
                             "None", "", gns, "ensembl_gene_id", FALSE, 
                             input$selectGOGnTb)
         # x <- runWebGestaltR(reVals$analysisOb, input$selectGODB, input$selectGOMethod,
         #                     "None", input$textSaveGO, gns, "ensembl_gene_id", FALSE, 
         #                     input$selectGOGnTb)
         removeNotification(nt)
         # Check if there was an error
         #browser()
         if (any(class(x) == "error")){
           
           showModal(modalDialog(title = "WebGStalt Error", x$message))
           
         } else
         {
           
           reVals$analysisOb <- x
           updateSelectInput(session, "selectPthWy", choices = names(reVals$analysisOb@PthWyAnl))
         }
         
       }
       )
       
       observeEvent(input$selectGOGnTb, {

         if (length(reVals$analysisOb@GeneTables) > 0)
         {
           #browser()
           tMet <- colnames(reVals$analysisOb@GeneTables[[input$selectGOGnTb]]@gnTbl)
           updateSelectInput(session, inputId = "selectGSEAEnrich", 
                           choices = tMet, selected = tMet[length(tMet)])
           }
       })
    
       # Plot a bar chart of pathway analysis
       observeEvent(input$pthWyBarPlot, {
         
        #browser()
         plt <- plotPthWybar(reVals$analysisOb@PthWyAnl[[input$selectPthWy]])
         #browser()
         output$pthWyPlt <- renderPlotly(plt)
       })
       
       # observeEvent(event_data("plotly_click", source = "A"), 
       #              {
       #                browser()
       #              })

       observeEvent(event_data("plotly_click", source = "C"), 
                    {
                      # output$pthWyPlt
                      #
                      #browser()
                      dx <- event_data("plotly_click", source = "C")$pointNumber + 1
                      g <- data.frame(strsplit(reVals$analysisOb@PthWyAnl[[input$selectPthWy]]@pthWyAnl$userId[[dx]], split = ";", fixed = T)[[1]])
                      colnames(g) <- c("Gene ID")
                      
                      if (grep("ENSG", g$`Gene ID`[[1]])){
                        
                        
                        gnSym <- ensembldb::select(EnsDb.Hsapiens.v86, keys = as.character(g$`Gene ID`), keytype = "GENEID", columns = c("SYMBOL"))
                        g$`Gene Symbol` <- gnSym$SYMBOL
                      }
                      
                      output$pthWyTb <- DT::renderDataTable(g, server = TRUE,  rownames = T)
                      url1 <- a(reVals$analysisOb@PthWyAnl[[input$selectPthWy]]@pthWyAnl$geneSet[[dx]] , 
                                        href = reVals$analysisOb@PthWyAnl[[input$selectPthWy]]@pthWyAnl$link[[dx]])
                      output$pthWyUrl <- renderUI({
                        tagList("URL link:", url1)
                      })
                      # output$pthWyTxt <- renderText(paste(
                      #                          reVals$analysisOb@PthWyAnl[[input$selectPthWy]]@pthWyAnl$link[[dx]]))
                    })
}
# Run the application 
shinyApp(ui = ui, server = server)
