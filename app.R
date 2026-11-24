# ___________________ ----
# App Startup ----

## 1.0 Load Libraries ----
library(EZRShiny)
library(shiny)
library(bslib)
library(shinyjs)
library(bsicons)
library(plotly)
library(DT)
library(readr)
library(tidytable)
library(colourpicker)
library(pheatmap)
library(grid)
library(ggnewscale)
library(stringr)
library(viridis)
library(svglite)
library(tibble)
library(splitstackshape)
library(ggh4x)
library(ggplot2)
library(ggsci)
library(ggbeeswarm)
library(ggprism)
library(lemon)
library(dplyr)
library(MASS)
library(emmeans)
library(sgpv)
library(data.table)
library(lme4)

## 2.0 Load Basics ----
options(shiny.maxRequestSize = 300 * 1024^2)
sourceFunctions("Functions")


## 3.0 Universal Vars ----
# App Name here:
appName <- "M2U"
# Necessary Files here:
tree <<- readRDS("Necessary Files/flippedTreeMouse.rds")
template <<- read.csv("Necessary Files/Non-Nutil Template.csv")
recCols <<- c("mouse","sex","treatment","mpi","genotype","marker","batch")
baseCols <<- c("fileName", "region", "ABAID", "hemi")
thePalette <<- readRDS("Necessary Files/colorPalette.rds")


# ___________________ ----
# UI ----
# Ensure use of document outline and preloading as many inputs as possible
ui <- UINav(

  # Logo files must be placed in www folder and must include ".{type}"
  logoFile = c("VAI 2 Line White.png", "Updated M2U Logo.jpg"),
  logoHeight = c("35vh", "40vh"),
  appName = appName,
  barColor = "#005596",

  ## 1.0 File Import ----
  biLevelTab( "File Import",

    ### 1.1 Initial Import ----
    subTab( "Initial Import",
      navDownload("nnTemplate", "Download the Template for Non-Nutil Data"),
      navUpload("fileInput", "Select All Files to Import", multiple = TRUE),
      navDownload("downloadRawData", "Download the Raw Data Save")
    ),

    ### 1.2 Annotation Creation ----
    subTab( "Annotation Creation",
      navUpload("rawAnnoInput", "Input Raw Data if needed", multiple = TRUE),
      navCheckbox("singleHemi", "Single Hemisphere (Defaults to right)", value = FALSE),
      navSelect("annoRecOptions", "Recommended Options",
                choices = recCols, selected = recCols, multiple = TRUE),
      navSelect("annoCustomOptions", "Custom Options", multiple = TRUE, create = TRUE),
      navButton("createAnnoFile", "Create Annotation File"),
      navDownload("downloadAnnoFile", "Download the Blank Annotation File")
    ),

    ### 1.3 Merge Data and Annotations ----
    subTab( "Merge Data and Annotations",
      navUpload("rawDataInput", "Select your Raw Data File"),
      navUpload("annoFileInput", "Select your completed Annotation File"),
      navDownload("downloadCheckpointData", "Download the Checkpoint for the Data")
    ),

    ### 1.4 Checkpoint Start ----
    subTab( "Checkpoint Start",
      navUpload("checkpointInput", "Select your Checkpoint File"),
      navUpload("varInput", "Select your Variable File")
    )
  ),

  ## 2.0 Variable Selection ----
  biLevelTab( "Variable Selection",

    ### 2.1 Layout Selection ----
    subTab( "Layout Selection",
      navSelect("AoI", "Select Annotation(s) of Interest (AoI)", multiple = TRUE),
      navSelect("regionLevel", "Select Region Level",
                choices = c("daughter", "parent", "minor", "major"), selected = "daughter"),
      navSelect("legendVar", "Select Legend Variables (based on selected region level)",
                choices = c("current level", "parent", "major"),
                selected = c("current level","major"), multiple = TRUE),
      navSelect("rowVars", "Select Row Names", multiple = TRUE),
      navSelect("colVars", "Select Column Names (based on selected region level)",
                choices = c("none","current level", "parent", "major"), selected = "none")
    ),

    ### 2.2 Data and Transformation ----
    subTab( "Data and Transformation",
      navSelect("valueVar", "Select Column to get Data from"),
      navCheckbox("logged", "Log Transform", value = TRUE),
      navCheckbox("multiPercent", "Multiply by 100", value = FALSE),
      navCheckbox("deviPercent", "Divide by 100", value = FALSE),
      navCheckbox("trim", "Trim", value = TRUE),
      navSelect("regionsToRemove", "Regions to ignore",
                choices = c("none",tree$region), selected = "none", multiple = TRUE)
    ),

    ### 2.3 Color Scheme ----
    subTab( "Color Scheme",
      navCheckbox("invert", "Invert the Color Scheme", value = FALSE),
      navOutputText("minMaxText"),
      navNumeric("minVal", "Set the Minimum Value of Scale", value = 0),
      navNumeric("maxVal", "Set the Maximum Value of Scale", value = 1),
      navSelect("colorPalette", "Select Color Palette",
                choices = c("Viridis", "2 Color", "3 Color")),
      uiOutput("colorOptions")
    ),

    ### 2.4 Variable Set and Save ----
    subTab( "Variable Set and Save",
      navButton("varSetButton", "Finalize the Variables"),
      navDownload("saveVar", "Download the Variable Save"),
      navDownload("normDataDownload", "Download Normal Plot Data")
    )
  ),

  ## 3.0 Data Preview ----
  biLevelTab( "Data Preview",

    ### 3.1 Normal Heatmap ----
    subSidebarTab( "Normal Heatmap",
      sidebarElements = list(
        navButton("normHeatmap", "Create normal heatmap"),
        navDownload("normHeatmapDownload", "Download the data preview heatmap"),
        navDownload("normJustPlotDownload", "Download the heatmap without the legend")
      ),
      navOutputPic("normPlot")
    ),

    ### 3.2 Anatomical Heatmap ----
    subSidebarTab( "Anatomical Heatmap",
      sidebarElements = list(
        navSelect("slices", "Select Allen Brain Atlas slices to plot",
                  choices = c(1:132), selected = c(32,46,67,75,82,96), multiple = TRUE),
        navButton("anatomicalHeatmap", "Create the anatomical heatmap"),
        navDownload("anatomicalHeatmapDownload", "Download the anatomical heatmap (click once and then wait)")
      ),
      navOutputPic("anatomicalPlot")
    )
  ),

  ## 4.0 Statistical Analysis ----
  biLevelTab( "Statistical Analysis",

    ### 4.1 Model Selection ----
    subTab( "Model Selection",
      navSelect("varInterest", "Select Variable of Interest",
                choices = c("notYet"), selected = "notYet"),
      navSelect("strata", "Select Strata",
                choices = c("notYet"), selected = "notYet", multiple = TRUE),
      navSelect("statMethod", "Select statistical method",
                choices = c("2nd Gen", "FDR"), selected = "2nd Gen"),
      navNumeric("sigLevel", "Select significance level for method", value = 0.05),
      navCheckbox("flipComp", "Flip Comparison", value = FALSE),
      navButton("runStats", "Run statistical analysis"),
      navDownload("statsDataDownload", "Pre Stats Download"),
      navDownload("statsDownload", "Stats Table Download"),
      navOutputText("statsError")
    ),

    ### 4.2 Statistical Anatomical Heatmap ----
    subSidebarTab( "Statistical Anatomical Heatmap",
      sidebarElements = list(
        navSelect("statSlices", "Select Allen Brain Atlas slices to plot",
                  choices = c(1:132), selected = c(32,46,67,75,82,96), multiple = TRUE),
        navSelect("statValueVar", "Select the value variable",
                  choices = c("Est. Difference","Adj. P-Value","Significance","Mean Difference"),
                  selected = "Est. Difference"),
        navNumeric("compTextSize", "Select comparison text size", value = 5),
        navButton("statAnatomicalHeatmap", "Create the anatomical heatmap"),
        navDownload("statAnatomicalHeatmapDownload", "Download the anatomical heatmap (click once and then wait)")
      ),
      navOutputPic("statAnatomicalPlot")
    ),

    ### 4.3 Per Region Boxplot ----
    subSidebarTab( "Per Region Boxplot",
      sidebarElements = list(
        navSelect("statBoxRegion", "Select a region to look at",
                  choices = tree$region, selected = tree$region[1]),
        navButton("runBoxPlot", "Create the boxplot"),
        navDownload("statBoxplotDownload", "Download the boxplot")
      ),
      navOutputPic("boxPlot")
    )

  )
)


# ___________________ ----
# Server ----
server <- function(session, input, output) {

  ## 1.0 Global Vars ----
  # Any variables here that need to carry across app, but should be local to the user
  # Try to keep all variables inside a list to help with debugging later
  global <- reactiveValues(
    dataFrames = list(
      rawData = NULL, annoFile = NULL, preVar = c(),
      fullData = c(), preStatsData = c(), 
      statsData = c(), MBH_data = c()
    ),
    info = list(
      numericCols = c(), recAnnoCols = c(), customAnnoCols = c(),
      annoCheck1 = F, annoCheck2 = F
    ),
    variables = list(
      #layout
      AoI = "", regionLevel = "", legendVar = c(), 
      rowVars = c(), colVars = c(),
      #dataTrans
      valueVar = c(), logged = F, multiPercent = F, deviPercent = F,
      trim = T, regionsNoData = c("fiber tracts"),
      #scaleOptions
      invert = F, minVal = 0, maxVal = 1,
      colorPalette = "",
      #colorScheme
      viridisOptions = "infirno",
      twoCol1 = "white", twoCol2 = "red",
      threeCol1 = "white", threeCol2 = "orange", threeCol3 = "red",
      #stats
      varInterest = c(), strata = c(), statMethod = "2nd Gen", sigLevel = 0.05,
      customStrata = c(), statValueVar = "Est. Difference"
    )
  )

  ## 2.0 Run on Start ----
  # Hide pages and deactivate buttons that should not be used yet
  observe({
    startSection("Run on Start")
    
    # Hide sub Tabs for tabs we don't want shown at start
    hideNavTabs(rootID = "Variable Selection", 
                tabIDs = c("Layout Selection",
                           "Data and Transformation",
                           "Color Scheme",
                           "Variable Set and Save"))
    
    hideNavTabs(rootID = "Data Preview",
                tabIDs = c("Normal Heatmap",
                           "Anatomical Heatmap"))
    
    hideNavTabs(rootID = "Statistical Analysis", 
                tabIDs = c("Model Selection",
                           "Statistical Anatomical Heatmap",
                           "Per Region Boxplot"))
    
    # Deactivate buttons that need other things to function
    deactivateItems(c("createAnnoFile",
                      "twoCol1","twoCol2","threeCol1","threeCol2","threeCol3",
                      "statAnatomicalHeatmap","runBoxPlot"))

    # Deactivate all downloads since there is nothing to download at start
    deactivateItems(c(
      "downloadRawData","downloadAnnoFile","downloadCheckpointData",
      "saveVar","normDataDownload","statsDataDownload","normHeatmapDownload",
      "normJustPlotDownload","anatomicalHeatmapDownload","statsDownload",
      "statAnatomicalHeatmapDownload","statBoxplotDownload"
    ))
    
    endSection("Run on Start")
  })

  ## 3.0 All File Imports ----
  # Make sure the initial download buttons that are needed have the proper files in them
  output$nnTemplate <- downloadHandler(
    filename = function() {
    paste("Non-Nutil Template.csv")
    },
    content = function(file) {
    write.csv(template, file, row.names = FALSE)
    }
  )

  # Make sure all places where file imports happen are reading in files appropriately

  ### 3.1 Basic File Import ----
  observeEvent(input$fileInput, {
    startSection("Basic File Import")

    # Inputs
    fileNames <- input$fileInput$name
    paths <- input$fileInput$datapath

    # Import all files
    importItems <- fileImport(fileNames,paths)
    
    # Pull out items
    rawData <- importItems[[1]]
    numericCols <- importItems[[2]]

    # App interactions
    enable("createAnnoFile")
    enable("downloadRawData")

    output$downloadRawData <- downloadHandler(
        filename = function() {
            paste("RawDataSave.csv")
        },
        content = function(file) {
          disable("downloadRawData")
          write.csv(rawData, file, row.names = FALSE)
        }
    )
    
    # Save globals
    global$info$numericCols <- numericCols
    global$dataFrames$rawData <- rawData
    global$info$annoCheck1 <- T

    endSection("Basic File Import")
  })

  ### 3.2 Raw Data for Annotation File Import ----
  observeEvent(input$rawAnnoInput, {
    startSection("Annotation Raw Data Save Import")
    
    # Read in raw data from the annotation page
    rawData <- read.csv(input$rawAnnoInput$datapath)
    numericCols <- colnames(rawData)[!colnames(rawData) %in% baseCols]

    # App interactions
    enable("createAnnoFile")

    #saved globals
    global$dataFrames$rawData <- rawData
    global$info$annoCheck1 <- T
    global$info$numericCols <- numericCols
    
    endSection("Annotation Raw Data Save Import")
  })

  ### 3.3 Annotation File Creation ----
  observeEvent(input$createAnnoFile, {
    startSection("Annotation File Creation")
    
    # Load globals
    rawData <- global$dataFrames$rawData

    # Inputs
    recOptions <- input$annoRecOptions
    customOptions <- input$annoCustomOptions
    singleHemi <- input$singleHemi
    
    # Create annotation file
    annoFile <- createAnnoFile(rawData, recOptions, customOptions, singleHemi)

    # App interactions
    enable("downloadAnnoFile")
    output$downloadAnnoFile <- downloadHandler(
      filename = function() {
        paste("AnnotationFile.csv")
      },
      content = function(file) {
        disable("downloadAnnoFile")
        write.csv(annoFile, file, row.names = FALSE)
      }
    )

    # Save globals
    global$dataFrames$annoFile <- annoFile
    
    endSection("Annotation File Creation")
  })

  ### 3.4 Raw Data for Merge ----
  observeEvent(input$rawDataInput, {
    startSection("Raw Data Merge Import")
    
    # Load globals
    annoCheck2 <- global$info$annoCheck2
    annoFile <- global$dataFrames$annoFile
    recAnnoCols <- global$info$recAnnoCols
    customAnnoCols <- global$info$customAnnoCols

    # Inputs
    rawDataPath <- input$rawDataInput$datapath

    # Read in raw data for merge
    rawData <- read.csv(rawDataPath)

    # Get cols for inputs
    annoCols <- c(recAnnoCols, customAnnoCols)
    numericCols <- colnames(rawData)[!colnames(rawData) %in% baseCols & 
                                       sapply(rawData, is.numeric) &
                                       !colnames(rawData) %in% annoCols]
    

    annoCheck1 <- T
    if (annoCheck1 & annoCheck2) {
      # Merge data and annots
      mergeItems <- mergeData(rawData, annoFile, numericCols)
      fullData <- mergeItems[[1]]
      numericCols <- mergeItems[[2]]
      
      # App interactions
      enable("downloadCheckpointData")

      output$downloadCheckpointData <- downloadHandler(
          filename = function() {
            paste("CheckpointData.csv")
          },
          content = function(file) {
            disable("downloadCheckpointData")
            fwrite(fullData, file, row.names = FALSE)
          }
      )
      
      showNavTabs(rootID = "Variable Selection", 
                  tabIDs = c("Layout Selection",
                           "Data and Transformation",
                           "Color Scheme",
                           "Variable Set and Save"))

      updateSelectizeInput(session, "valueVar", choices = numericCols, selected = numericCols[1])
      updateSelectizeInput(session, "AoI", choices = annoCols, selected = annoCols[1])
      
      # Save globals
      global$dataFrames$preVar <- fullData
      global$info$numericCols <- numericCols
      global$variables$valueVar <- numericCols[1]
    }

    # Save globals
    global$dataFrames$rawData <- rawData
    global$info$annoCheck1 <- annoCheck1
    global$info$numericCols <- numericCols
    
    endSection("Raw Data Merge Import")
  })
  
  ### 3.5 Annotation File Import ----
  observeEvent(input$annoFileInput, {
    startSection("Annotation File Merge Import")
    
    # Load globals
    annoCheck1 <- global$info$annoCheck1
    rawData <- global$dataFrames$rawData
    numericCols <- global$info$numericCols

    # Inputs
    annoDataPath <- input$annoFileInput$datapath

    # Read in annotation file for merge
    annoFile <- read.csv(annoDataPath)

    # Get cols for inputs
    recAnnoCols <- colnames(annoFile)[colnames(annoFile) %in% recCols]
    customAnnoCols <- colnames(annoFile)[!colnames(annoFile) %in% c(recCols,"include","fileName","hemi")]
    annoCols <- c(recAnnoCols, customAnnoCols)
    
    numericCols <- colnames(rawData)[!colnames(rawData) %in% baseCols & 
                                       sapply(rawData, is.numeric) &
                                       !colnames(rawData) %in% annoCols]
    
    annoCheck2 <- T
    if (annoCheck1 & annoCheck2) {
      # Merge data and annots
      mergeItems <- mergeData(rawData, annoFile, numericCols)
      fullData <- mergeItems[[1]]
      numericCols <- mergeItems[[2]]
      
      # App interactions
      enable("downloadCheckpointData")
      output$downloadCheckpointData <- downloadHandler(
          filename = function() {
            paste("CheckpointData.csv")
          },
          content = function(file) {
            disable("downloadCheckpointData")
            fwrite(fullData, file, row.names = FALSE)
          }
      )
  
      showNavTabs(rootID = "Variable Selection", 
                  tabIDs = c("Layout Selection",
                           "Data and Transformation",
                           "Color Scheme",
                           "Variable Set and Save"))
  
      updateSelectizeInput(session, "valueVar", choices = numericCols, selected = numericCols[1])
      updateSelectizeInput(session, "AoI", choices = annoCols, selected = annoCols[1])
      
      # Save globals
      global$dataFrames$preVar <- fullData
      global$info$numericCols <- numericCols
      global$variables$valueVar <- numericCols[1]
    }

    # Save globals
    global$variables$AoI <- annoCols[1]
    global$dataFrames$annoFile <- annoFile
    global$info$annoCheck2 <- annoCheck2
    global$info$recAnnoCols <- recAnnoCols
    global$info$customAnnoCols <- customAnnoCols
    
    endSection("Annotation File Merge Import")
  })

  ### 3.6 Checkpoint Load ----
  observeEvent(input$checkpointInput, {
    startSection("Checkpoint Load")
    
    # Inputs
    fullDataPath <- input$checkpointInput$datapath

    # Load checkpoint data
    fullData <- read.csv(fullDataPath)
    
    checkpointData <- loadCheckpointData(fullData)
    fullData1 <- checkpointData[[1]]
    numericCols <- checkpointData[[2]]
    recAnnoCols <- checkpointData[[3]]
    customAnnoCols <- checkpointData[[4]]
    annoCols <- c(recAnnoCols, customAnnoCols)
    

    # App interactions
    updateSelectizeInput(session, "AoI", choices = annoCols, selected = annoCols[1])
    updateSelectizeInput(session, "rowVars", choices = c("none",annoCols[1]), selected = "none")
    updateSelectizeInput(session, "valueVar", choices = c("none",numericCols), selected = numericCols[1])

    # Save globals
    global$dataFrames$preVar <- fullData1
    global$info$numericCols <- numericCols
    global$info$recAnnoCols <- recAnnoCols
    global$info$customAnnoCols <- customAnnoCols
    global$variables$AoI <- annoCols[1]
    global$variables$rowVars <- "none"
    global$variables$valueVar <- numericCols[1]
    global$variables$colorPalette <- "Viridis"
    output$colorOptions <- renderUI({
      tagList(
        selectizeInput(width = "100%","viridisOptions","Select Viridis Options",
                       choices = c("viridis","inferno","magma","plasma",
                                   "cividis","mako","rocket","turbo"), selected = "inferno",
                       multiple = F,options = list(create = F))
      )
    })
    
    global$variables$twoCol1 <- "white"
    global$variables$twoCol2 <- "red"
    global$variables$threeCol1 <- "white"
    global$variables$threeCol2 <- "orange"
    global$variables$threeCol3 <- "red"

    showNavTabs(rootID = "Variable Selection", 
                tabIDs = c("Layout Selection",
                         "Data and Transformation",
                         "Color Scheme",
                         "Variable Set and Save"))
    endSection("Checkpoint Load")
  })

  ### 3.7 Variable Load ----
  observeEvent(input$varInput, {
    startSection("Variable Load")

    # Inputs
    varInputPath <- input$varInput$datapath

    # Read in variable file
    variables <- readRDS(varInputPath)

    # Layout
    updateSelectizeInput(session, "AoI", choices = recCols, selected = variables$AoI)
    updateSelectizeInput(session, "regionLevel", selected = variables$regionLevel)
    updateSelectizeInput(session, "legendVar", selected = variables$legendVar)
    updateSelectizeInput(session, "rowVars", selected = variables$rowVars)
    updateSelectizeInput(session, "colVars", selected = variables$colVars)
    
    # Data Transformation
    updateSelectizeInput(session, "valueVar", selected = variables$valueVar)
    update_switch(session = session, "logged", value = variables$logged)
    update_switch(session = session, "multiPercent", value = variables$multiPercent)
    update_switch(session = session, "deviPercent", value = variables$deviPercent)
    update_switch(session = session, "trim", value = variables$trim)
    updateSelectizeInput(session, "regionsToRemove", selected = variables$regionsToRemove)

    # Color Scheme
    update_switch(session = session, "invert", value = variables$invert)
    updateNumericInput(session, "minVal", value = variables$minVal)
    updateNumericInput(session, "maxVal", value = variables$maxVal)

    # Color Palette Options
    updateSelectizeInput(session, "colorPalette", selected = variables$colorPalette)
    
    # Ensure that old version of var save get updated color options
    if (is.null(variables$viridisOptions)){ variables$viridisOptions <- "inferno"}
    if (is.null(variables$twoCol1)){ variables$twoCol1 <- "white"}
    if (is.null(variables$twoCol2)){ variables$twoCol2 <- "red"}
    if (is.null(variables$threeCol1)){ variables$threeCol1 <- "white"}
    if (is.null(variables$threeCol2)){ variables$threeCol2 <- "orange"}
    if (is.null(variables$threeCol3)){ variables$threeCol3 <- "red"}
    
    # Render correct color palette
    output$colorOptions <- renderUI({
      if (variables$colorPalette == "Viridis"){
        tagList(
          selectizeInput(width = "100%","viridisOptions","Select Viridis Options",
                         choices = c("viridis","inferno","magma","plasma",
                                     "cividis","mako","rocket","turbo"),selected = variables$viridisOptions,
                         multiple = F,options = list(create = F))
        )
          
      } else if (variables$colorPalette == "2 Color"){
        tagList(
          colourInput(width = "100%","twoCol1", "Select colour", variables$twoCol1),
          colourInput(width = "100%","twoCol2", "Select colour", variables$twoCol2)
        )
          
      } else if (variables$colorPalette == "3 Color"){
        tagList(
          colourInput(width = "100%","threeCol1", "Select colour", variables$threeCol1),
          colourInput(width = "100%","threeCol2", "Select colour", variables$threeCol2),
          colourInput(width = "100%","threeCol3", "Select colour", variables$threeCol3)
        )
      } else {
        tagList(
          selectizeInput(width = "100%","viridisOptions","Select Viridis Options",
                         choices = c("viridis","inferno","magma","plasma",
                                     "cividis","mako","rocket","turbo"),selected = "inferno",
                         multiple = F,options = list(create = F))
        )
      }
    })
    
    # Save globals
    saveVars <- variables
    global$variables <- variables

    endSection("Variable Load")
  })

  ## 4.0 Variable Selection Update ----
  # Should run whenever a variable input changes and when app starts
  observe({
    startSection("Variable Selection Update")
    
    # Do not run until variables actually exist
    if (!is.null(input$AoI)){
      # Create variable list using updated inputs
      variables <- list(
        #layout
        AoI = input$AoI, regionLevel = input$regionLevel, legendVar = input$legendVar, 
        rowVars = input$rowVars, colVars = input$colVars,
        #dataTrans
        valueVar = input$valueVar, logged = input$logged, multiPercent = input$multiPercent, deviPercent = input$deviPercent,
        trim = input$trim, regionsNoData = c("fiber tracts"),
        #scaleOptions
        invert = input$invert, minVal = input$minVal, maxVal = input$maxVal,
        colorPalette = input$colorPalette,
        #colorScheme
        viridisOptions = input$viridisOptions,
        twoCol1 = input$twoCol1, twoCol2 = input$twoCol2,
        threeCol1 = input$threeCol1, threeCol2 = input$threeCol2, threeCol3 = input$threeCol3,
        #stats
        varInterest = input$varInterest, strata = input$strata, statMethod = input$statMethod, sigLevel = input$sigLevel,
        customStrata = input$customStrata, statValueVar = input$statValueVar
      )
      
      # Update legendVar if regionLevel changed
      if (input$regionLevel != global$variables$regionLevel){
        print("Update regionLevel means update legendVar")
        if (input$regionLevel == "daughter"){
          updateSelectizeInput(session, "legendVar", choices = c("current level","parent","major"), selected = c("current level","major"))
        } else if (input$regionLevel == "parent"){
          updateSelectizeInput(session, "legendVar", choices = c("current level","major"), selected = c("current level","major"))
        } else if (input$regionLevel == "minor"){
          updateSelectizeInput(session, "legendVar", choices = c("current level","major"), selected = c("current level","major"))
        } else {
          updateSelectizeInput(session, "legendVar", choices = c("current level"), selected = c("current level"))
        }
        
      }
  
      # Update rowVars input if AOI changed
      if (any(input$AoI != global$variables$AoI)){
        print("Update AoI means update RowVars")
        updateSelectizeInput(session, "rowVars", choices = c("none",input$AoI), selected = "none")
      }
      
      # Make sure all color inputs have something
      if (is.null(variables$viridisOptions)){ variables$viridisOptions <- "inferno"}
      if (is.null(variables$twoCol1)){ variables$twoCol1 <- "white"}
      if (is.null(variables$twoCol2)){ variables$twoCol2 <- "red"}
      if (is.null(variables$threeCol1)){ variables$threeCol1 <- "white"}
      if (is.null(variables$threeCol2)){ variables$threeCol2 <- "orange"}
      if (is.null(variables$threeCol3)){ variables$threeCol3 <- "red"}
  
      # Update Color Inputs if overall type changed
      if (input$colorPalette != global$variables$colorPalette){
        print("Color palette update")
        output$colorOptions <- renderUI({
          if (input$colorPalette == "Viridis"){
            tagList(
              selectizeInput(width = "100%","viridisOptions","Select Viridis Options",
                             choices = c("viridis","inferno","magma","plasma",
                                         "cividis","mako","rocket","turbo"),
                             selected = variables$viridisOptions,
                             multiple = F,options = list(create = F))
            )
              
          } else if (input$colorPalette == "2 Color"){
            tagList(
              colourInput(width = "100%","twoCol1", "Select colour", 
                          variables$twoCol1),
              colourInput(width = "100%","twoCol2", "Select colour", 
                          variables$twoCol2)
            )
              
          } else {
            tagList(
              colourInput(width = "100%","threeCol1", "Select colour", 
                          variables$threeCol1),
              colourInput(width = "100%","threeCol2", "Select colour", 
                          variables$threeCol2),
              colourInput(width = "100%","threeCol3", "Select colour", 
                          variables$threeCol3)
            )
          }
        })
      }
    
      # Save to global
      global$variables <- variables
    } else {
      print("Variables don't exist yet, so ignore")
    }

    endSection("Variable Selection Update")
  })

  
  ## 5.0 Variable Set ----

  #Set Variables
  observeEvent(input$varSetButton, {
    startSection("Variable Set")
    
    # Load globals
    preVar <- global$dataFrames$preVar
    variables <- global$variables 

    # Inputs
    strataInput <- input$strata
    
    # Calculate value table
    AoIStrata <- variables$AoI
    valueTable <- calculateValueTable(preVar,variables,AoIStrata)

    # Get statistical analysis options
    statAnnoCols <- recCols[!recCols %in% c("mouse","batch")]
    statAnnoCols <- statAnnoCols[statAnnoCols %in% colnames(preVar)]

    # Update statistical analysis options
    StatsItems <- UpdateStats(preVar, statAnnoCols, strataInput, variables)
    officialStats <- StatsItems[[1]]
    strataOptions <- StatsItems[[2]]
    variables <- StatsItems[[3]]
    
    # App interactions
    activateItems(c("saveVar","downloadAnnoFile","normDataDownload"))
    
    output$saveVar <- downloadHandler(
      filename = function() {
        paste("VariableSave.rds")
      },
      content = function(file) {
        disable("saveVar")
        saveRDS(variables, file)
      }
    )
    
    showNavTabs(rootID = "Data Preview", 
                tabIDs = c("Normal Heatmap",
                           "Anatomical Heatmap"))

    #update text above min and max
    output$minMaxText <- renderText({
        paste0("Min Value: ", min(valueTable$preLimVal, na.rm = T), 
              " | Max Value: ", max(valueTable$preLimVal, na.rm = T))
    })
    
    output$normDataDownload <- downloadHandler(
        filename = function() {
            paste("normData.csv")
        },
        content = function(file) {
          disable("normDataDownload")
          write.csv(valueTable, file, row.names = FALSE)
        }
    )

    showNavTabs(rootID = "Statistical Analysis",
                tabIDs = c("Model Selection",
                            "Statistical Anatomical Heatmap",
                            "Per Region Boxplot"))
    
    updateSelectizeInput(session, "strata", choices = strataOptions, selected = variables$strata)
    updateSelectizeInput(session, "varInterest", choices = officialStats, selected = variables$varInterest)

    # The variables are set, so send the user straight to the normal heatmap
    showNavTabs(rootID = "root", tabIDs = c("Data Preview"))

    # Save globals
    global$dataFrames$fullData <- valueTable
    global$variables <- variables
    
    endSection("Variable Set")
    
  })


  ## 6.0 Data Views ----

  ### 6.1 Basic Heatmap ----
  observeEvent(input$normHeatmap, {
    startSection("Normal Heatmap")
    
    # Load globals
    fullData <- global$dataFrames$fullData
    variables <- global$variables
    
    # Make base plots
    theHeatmaps <- makeNormalHeatmaps(fullData, variables)
    data <- theHeatmaps[[1]]
    displayBasePlot <- theHeatmaps[[2]]
    fullPlot <- theHeatmaps[[3]]

    heatmapWidth <- length(colnames(data))/10
    
    # App interactions
    
    # Save and show the display plot
    ggsave(paste0(session$token,"_NormalHeatmap.png"), plot = displayBasePlot, width = heatmapWidth, height = 5,limitsize = FALSE)
    output$normPlot <- renderImage({
      list(src = paste0(session$token,"_NormalHeatmap.png"), alt = "plot wasn't made")
    }, deleteFile = T)
    
    # Enable and make the download for the full plot
    enable("normHeatmapDownload")
    output$normHeatmapDownload <- downloadHandler(
      filename = function() {
        paste0("heatmap", ".pdf")
      },
      content = function(file) {
        disable("normHeatmapDownload")
        ggsave(file, plot = fullPlot, width = heatmapWidth, height = 10, limitsize = FALSE)
        print("Done Downloading Heatmap")
      }
    )

    # Enable and make the download for the base plot
    enable("normJustPlotDownload")
    output$normJustPlotDownload <- downloadHandler(
      filename = function() {
        paste0("justPlot", ".pdf")
      },
      content = function(file) {
        disable("normJustPlotDownload")
        ggsave(file, plot = displayBasePlot, width = heatmapWidth, height = 10, limitsize = FALSE)
        print("Done Downloading Heatmap")
      }
    )
    
    endSection("Normal Heatmap")
  })


  ### 6.2 Anatomical Heatmap ----
  observeEvent(input$anatomicalHeatmap, {
    startSection("Anatomical Heatmap")
  
    # Load globals
    valueTable <- global$dataFrames$fullData
    variables <- global$variables

    # Inputs
    slices <- input$slices
    
    valueTable1 <- valueTable
    
    # Create Anatomical Plots
    theHeatmap <- makeAnatomicalHeatmaps(valueTable, slices, variables)
    heatmapHeight <- length(unique(valueTable$x)) + 2
    heatmapWidth <- length(slices) + 1
    
    ggsave(paste0(session$token,"_AnatomicalHeatmap.png"), theHeatmap, height = heatmapHeight, width = heatmapWidth,limitsize = FALSE)
    
    # App interactions
    output$anatomicalPlot <- renderImage({
      list(src = paste0(session$token,"_AnatomicalHeatmap.png"), alt = "plot wasn't made")
    }, deleteFile = T)
    
    enable("anatomicalHeatmapDownload")
    output$anatomicalHeatmapDownload <- downloadHandler(
      filename = function() {
        paste0("anatomicalHeatmap", ".svg")
      },
      content = function(file) {
        disable("anatomicalHeatmapDownload")
        ggsave(file, theHeatmap, height = length(unique(valueTable$x)) + 2, width = length(slices) + 1,limitsize = FALSE)
        print("Done Downloading Heatmap")
      }
    )
    
    endSection("Anatomical Heatmap")
  })
  
  
  ## 7.0 Stats ----
  
  ### 7.1 Calculations ----
  observeEvent(input$runStats,{
    startSection("Statistical Analysis")
    
    # Load globals
    preVar <- global$dataFrames$preVar
    variables <- global$variables

    # Inputs
    sigLevel <- input$sigLevel
    flipComp <- input$flipComp
    statMethod <- input$statMethod
    
    # Create variables for use in stats functions
    VoI <- variables$varInterest
    regionsToRemove <- variables$regionsToRemove
    strataCols <- variables$strata
    AoIStrata <- strataCols
    comparisonCols <- VoI
    
    # Adjust strata cols if empty
    if (is.null(strataCols)){
      strataCols <- "None"
    }
    
    # Calculate value table for stats
    fullData <- calculateValueTable(preVar,variables,AoIStrata,VoI)

    # Clean the dataframe including adjusting columns for the stats function
    MBH_data <- CleanData(fullData, comparisonCols, strataCols, flipComp)
    
    # Set up stats options based on method chosen
    if(statMethod == "2nd Gen"){
      null.interval <- sigLevel
      fdr_cutoff <- NULL
    } else {
      fdr_cutoff <- sigLevel
      null.interval <- NULL
    }
    
    # Run stats
    statsTable <- GetRegionStats(MBH_data, null.interval = null.interval, fdr_cutoff = fdr_cutoff)
    
    # Make sure stats table ran correctly
    if (any(statsTable != F)){

      # Check for errors in the table and stop if there are any
      if (length(statsTable)< 2){
        stop(print(statsTable))
      }

      # Get and add manual mean difference value to stats table
      statsTable <- MeanDifferenceValueCalc(statsTable, MBH_data)
      
      # App interactions
      updateSelectizeInput(session, "statBoxRegion", choices = unique(MBH_data$region))
      
      # Save globals
      # Both significance plots read these back, so they must be stored here
      global$dataFrames$statsData <- statsTable
      global$dataFrames$MBH_data <- MBH_data
      
      activateItems(c("statAnatomicalHeatmap","runBoxPlot","statsDownload"))
      
      # The stats finished, so send the user straight to the anatomical heatmap
      showNavTabs(rootID = "Statistical Analysis",
                  tabIDs = c("Statistical Anatomical Heatmap"))
      
      enable("statsDownload")
      output$statsDownload <- downloadHandler(
        filename = function() {
          paste("statsTable.csv")
        },
        content = function(file) {
          disable("statsDownload")
          write.csv(statsTable, file, row.names = FALSE)
        }
      )
      
      enable("statsDataDownload")
      output$statsDataDownload <- downloadHandler(
        filename = function() {
          paste("preStats.csv")
        },
        content = function(file) {
          disable("statsDataDownload")
          write.csv(MBH_data, file, row.names = FALSE)
        }
      )
      output$statsError <- output$text <- renderText({ "" })
    } else {
      output$statsError <- output$text <- renderText({ "Not enough overlap in regions and Variable of Interest and/or Strata." })
    }
    
    endSection("Statistical Analysis")
  })

  
  ### 7.2 Anatomical Heatmaps ----
  observeEvent(input$statAnatomicalHeatmap,{
    startSection("Statistical Anatomical Heatmap")
    
    # Load globals
    valueTable <- global$dataFrames$statsData
    variables <- global$variables

    # Inputs
    slices <- input$statSlices
    compTextSize <- input$compTextSize
    statValue <- input$statValueVar
    
    # Figure out what value to pull for the plot
    if (statValue == "Est. Difference"){
      statValue <- "estimate"
      isSig = F
    } else if (statValue == "Adj. P-Value"){
      statValue <- "p.adj"
      isSig = F
    } else if (statValue == "Mean Difference"){
      statValue <- "contrastValue"
      isSig = F
    } else {
      statValue <- "Significant"
      isSig <- T
    }
    
    # Adjust value table
    valueTable <- as.data.frame(valueTable)
    valueTable <- valueTable[,colnames(valueTable) %in% c("contrast","strata","region","hemi",statValue)]
    
    colnames(valueTable)[colnames(valueTable) == statValue] <- "postLimVal"
    colnames(valueTable)[colnames(valueTable) == "contrast"] <- "x"
    
    # The stats table carries the strata as one "strata" column, so key the
    # label off whether that column actually holds anything. Looking for the
    # original strata column names here never matched, which dropped the
    # strata from the facet labels and collapsed the rows together.
    if (!is.null(valueTable$strata) && any(nzchar(valueTable$strata))){
      valueTable$annots <- paste0(valueTable$x," | ",valueTable$strata)
    } else {
      valueTable$annots <- valueTable$x
    }
    
    valueTable$annots <- gsub("sex_","",valueTable$annots)
    valueTable$removeRegion <- F
    
    # split Strata col into cols later in code (find annoCols and swap to strataCols)
    annoCols <- "annots"
    theHeatmap <- makeAnatomicalHeatmaps(valueTable, slices, variables, annoCols, isSig = isSig)
    
    # One facet row per annotation, which is the contrast crossed with the
    # strata, so size the figure off those rather than the contrasts alone
    heatmapHeight <- length(unique(valueTable$annots)) + 2
    heatmapWidth <- length(slices) + 1
    
    ggsave(paste0(session$token,"_StatAnatomicalHeatmap.png"), theHeatmap, height = heatmapHeight, width = heatmapWidth,limitsize = FALSE)
    
    output$statAnatomicalPlot <- renderImage({
      list(src = paste0(session$token,"_StatAnatomicalHeatmap.png"), alt = "plot wasn't made")
    }, deleteFile = T)
    
    enable("statAnatomicalHeatmapDownload")
    output$statAnatomicalHeatmapDownload <- downloadHandler(
      filename = function() {
        paste0("statsAnatomicalHeatmap", ".pdf")
      },
      content = function(file) {
        disable("statAnatomicalHeatmapDownload")
        ggsave(file, theHeatmap, height = heatmapHeight, width = heatmapWidth,limitsize = FALSE)
        print("Done Downloading Heatmap")
      }
    )
    
    endSection("Statistical Anatomical Heatmap")
  })
  
  
  ### 7.3 Box Plots ----
  observeEvent(input$runBoxPlot,{
    startSection("Box Plot")
    
    MBH_data <- global$dataFrames$MBH_data
    res_df <- global$dataFrames$statsData
    
    region <- input$statBoxRegion
    
    #remove duplicates in MBH_data
    MBH_data <- MBH_data[!(duplicated(MBH_data$region) & duplicated(MBH_data$value) & duplicated(MBH_data$hemi)), ]
    
    region_id <- region
    theBoxPlot <- MakeBoxPlot(MBH_data,res_df,region)
    
    plotWidth <- length(unique(res_df$contrast)) + length(unique(res_df$strata)) + 2
    ggsave(paste0(session$token,"_StatBoxPlot.png"), theBoxPlot, width = plotWidth,limitsize = FALSE)
    
    output$boxPlot <- renderImage({
      list(src = paste0(session$token,"_StatBoxPlot.png"), alt = "plot wasn't made")
    }, deleteFile = T)
    
    enable("statBoxplotDownload")
    output$statBoxplotDownload <- downloadHandler(
      filename = function() {
        paste0("statBoxPlot", ".pdf")
      },
      content = function(file) {
        disable("statBoxplotDownload")
        ggsave(file, theBoxPlot, width = plotWidth,limitsize = FALSE)
        print("Done Downloading Box Plot")
      }
    )
    
    endSection("Box Plot")
  })

}


# . ----
# 4.0 Run App ----
shinyApp(ui = ui, server = server)



