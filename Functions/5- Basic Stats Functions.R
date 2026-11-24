# ___________________ ----

UpdateStats <- function(preVar, statAnnoCols, strataInput, variables){
  # If stats have not yet been set, set them here
  if (strataInput == "notYet"){
    officialStats <- c()
    if (length(statAnnoCols) > 0){
      for (i in 1:length(statAnnoCols)){
        if (length(unique(preVar[,statAnnoCols[i]])) > 1){
          officialStats <- c(officialStats,statAnnoCols[i])
        }
      }
      if (length(officialStats) > 1){
        variables$strata <- officialStats[1]
        variables$varInterest <- officialStats[2]
        strataOptions <- c("None",officialStats)
          
      } else if (length(officialStats) == 1){
        variables$strata <- "None"
        variables$varInterest <- officialStats[1]
        strataOptions <- c("None")
      }
    }

  # Update stats options
  } else {
    officialStats <- c()
    for (i in 1:length(statAnnoCols)){
      if (length(unique(preVar[,statAnnoCols[i]])) > 1){
        officialStats <- c(officialStats,statAnnoCols[i])
      }
    }
    strataOptions <- c("None",officialStats)
  }

  return(list(officialStats, strataOptions, variables))
}


# ___________________ ----

CleanData <- function(fullData, comparisonCols, strataCols, flipComp){
  MBH_data <- fullData
  MBH_data <- as.data.frame(MBH_data)

  MBH_data$hemi <- tolower(MBH_data$hemi)
  MBH_data <- MBH_data[!MBH_data$removeRegion,]

  # Notes for later:
  # MBH_data: the data frame for MBH (expected columns: value, hemi, comparison, region, strata)
  # strata: c("treatment","mpi","genotype","marker", customStrata) (paste with sep = "_")

  # Make sure sex column is in the MBH data if it is in the comparisons or strata
  if ("sex" %in% c(comparisonCols)){
    MBH_data <- add_column(MBH_data, sexCol = "sex", .before = "sex")
    sexLocation <- which(comparisonCols == "sex")
    comparisonCols <- append(comparisonCols, "sexCol", after=sexLocation-1)
  } else if ("sex" %in% c(strataCols)){
    MBH_data <- add_column(MBH_data, sexCol = "sex", .before = "sex")
    sexLocation <- which(strataCols == "sex")
    strataCols <- append(strataCols, "sexCol", after=sexLocation-1)
  }

  # Combine comparison and strata columns if there are more than 1 of either, so that stats can be run on the combinations
  if (length(comparisonCols) > 1) {
    MBH_data$comparison <- apply(MBH_data[,comparisonCols], 1, paste, collapse = "_")
  } else {
    MBH_data$comparison <- MBH_data[,comparisonCols]
  }
  MBH_data <- as.data.frame(MBH_data)

  # Same for strata, but also account for none option
  if (length(strataCols) > 1) {
    MBH_data$strata <- apply(MBH_data[,strataCols], 1, paste, collapse = "_")
  } else if (strataCols != "None"){
    MBH_data$strata <- MBH_data[,strataCols]
  }
  MBH_data <- as.data.frame(MBH_data)

  # Rename preLimVal for clarification
  MBH_data$value <- MBH_data$preLimVal

  # Only keep necessary cols and make a seperate dataset because its had issues pulling the old dataset before
  MBHKeepCols <- c("region","hemi","comparison","strata","value","sex","batch")
  if("sex" %in% comparisonCols | "sex" %in% strataCols){
    MBHKeepCols <- c("region","hemi","comparison","strata","value","batch")
  }
  MBH_data1 <- MBH_data[,colnames(MBH_data) %in% MBHKeepCols]

  # Make sure to save the dataset as backup
  MBH_data <- MBH_data1

  # Order the levels of the comparison factor based on if it needs to be flipped or not
  orderOfVoI <- levels(as.factor(MBH_data1$comparison))
  if (flipComp == T){
    orderOfVoI <- rev(orderOfVoI)
  }
  MBH_data1$comparison <- factor(MBH_data1$comparison,levels = orderOfVoI)

  # Make a unique ID for each region, hemi, and comparison combination for stats
  MBH_data1$uniqueID <- paste0(MBH_data1$region,"_",MBH_data1$hemi,"_",MBH_data1$comparison)

  return(MBH_data1)
}


# ___________________ ----
MeanDifferenceValueCalc <- function(statsTable, MBH_data){

  # Get comparisons
  comparisonsBase <- data.frame("contrast" = statsTable$contrast,
                                "tochange" = statsTable$contrast,
                                "region" = statsTable$region,
                                "hemi" = statsTable$hemi)
  comparisonsBase <- tidyr::separate(data = comparisonsBase, col = tochange, 
                                      into = c("var1","var2"), sep = " - ")

  comparisonDiff <- MBH_data %>% group_by(uniqueID) %>%
    summarise("meanValue" = mean(value),
              "region" = unique(region),
              "hemi" = unique(hemi),
              "var1" = unique(comparison),
              "var2" = unique(comparison))

  varValues <- left_join(comparisonsBase[,c("contrast","var1","var2","region","hemi")],
                          comparisonDiff[,c("meanValue","region","hemi","var1")], 
                          by = c("region","hemi","var1"))
  colnames(varValues) <- c("contrast","var1","var2","region","hemi","var1Value")

  varValues2 <- left_join(comparisonsBase[,c("contrast","var1","var2","region","hemi")],
                          comparisonDiff[,c("meanValue","region","hemi","var2")], 
                          by = c("region","hemi","var2"))
  colnames(varValues2) <- c("contrast","var1","var2","region","hemi","var2Value")

  varValues <- left_join(varValues,varValues2)

  varValues$contrastValue <- varValues$var1Value - varValues$var2Value

  varValues <- varValues[,c("contrast","region","hemi","contrastValue")]

  # Add manual mean difference value to stats table
  statsTable <- left_join(statsTable,varValues)

  # Remove empty contrasts
  statsTable <- statsTable[! is.na(statsTable$contrast),]

  return(statsTable)
}