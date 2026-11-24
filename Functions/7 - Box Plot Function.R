# This function takes the results file produced by the N2U analysis and plots a region of interest

# imports: dplyr, ggplot2, ggsci, ggbeeswarm, ggprism
# these are loaded in app.R

# Expected parameters
## MBH_data: the raw data for MBH
## res_df: the results from the analysis portion of the Shiny app
## region_id: the desired region to plot this is based on a UI input
## color_palette: a user supplied color palette, if one isn't provided the nejm palette is the default

## Random example data.frame can be created by running create_test_file.R and then model_regions.R

# Convenience function
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

MakeBoxPlot <- function(MBH_data,res_df,region, color_palette=NULL){
  if(! any(colnames(MBH_data) == "strata")){
    MBH_data$strata = ""
  }
  MBH_data$strata = as.character(MBH_data$strata)
  plot.MBH = MBH_data[MBH_data$region == region,]
  plot.MBH$loop_strata = paste0(plot.MBH$region,plot.MBH$hemi,plot.MBH$strata)
  sig_df = res_df[res_df$region == region,]
  sig_df$loop_strata = paste0(sig_df$region,sig_df$hemi,sig_df$strata)
  sig_df$strata = as.character(sig_df$strata)
  plot.MBH$strata = gsub("sex_","",plot.MBH$strata)
  plot.MBH$comparison = gsub("sex_","",plot.MBH$comparison)
  if(length(unique(MBH_data$strata))>1){
    plot.MBH = concat.split(data = plot.MBH, split.col = "strata", sep = "_")
  }
  plot.MBH$comparison = as.character(plot.MBH$comparison)
  
  plot.region = ggplot(data = plot.MBH, aes(x = hemi,
                                            y = value,
                                            color = comparison,
                                            group=interaction(hemi,comparison))) +
    geom_quasirandom(dodge.width=.75) +
    geom_boxplot(fill=NA,outliers = F,show.legend = F) +
    #We should probably add an option for users to set the name of this
    ylab("Value") +
    xlab("Hemisphere") +
    ggtitle(region) +
    theme_prism()
  
  if(length(which(grepl("strata_",fixed=T,colnames(plot.MBH)))) == 1){ 
    plot.MBH$strata_1 = as.character(plot.MBH$strata_1)
    plot.region = plot.region + facet_rep_grid(~strata_1, scales = "free")
  }
  
  if(length(which(grepl("strata_",fixed=T,colnames(plot.MBH)))) == 2){ 
    plot.MBH$strata_1 = as.character(plot.MBH$strata_1)
    plot.MBH$strata_2 = as.character(plot.MBH$strata_2)
    plot.region = plot.region + facet_rep_grid(strata_2~strata_1, scales = "free")
  }
  
  if(length(which(grepl("strata_",fixed=T,colnames(plot.MBH)))) == 3){ 
    plot.MBH$strata_1 = as.character(plot.MBH$strata_1)
    plot.MBH$strata_2 = as.character(plot.MBH$strata_2)
    plot.MBH$strata_3 = as.character(plot.MBH$strata_3)
    plot.region = plot.region + facet_rep_grid(strata_2+strata_3~strata_1, scales = "free")
  }
  
  if(length(which(grepl("strata_",fixed=T,colnames(plot.MBH)))) == 4){ 
    plot.MBH$strata_1 = as.character(plot.MBH$strata_1)
    plot.MBH$strata_2 = as.character(plot.MBH$strata_2)
    plot.MBH$strata_3 = as.character(plot.MBH$strata_3)
    plot.MBH$strata_4 = as.character(plot.MBH$strata_4)
    plot.region = plot.region + facet_rep_grid(strata_2+strata_3~strata_1 + strata_4, scales = "free")
  }
  
  if(is.null(color_palette)){
    plot.region = plot.region + scale_color_nejm()
  } else {
    plot.region = plot.region + scale_color_manual(values = color_palette)
  }
  
  if(!is.null(MBH_data$batch)){
    plot.region = plot.region + aes(alpha = batch)
  }
  
  if(!is.null(MBH_data$sex)){
    plot.region = plot.region + aes(shape = sex)
  }
  
  sig_df$label = case_when(sig_df$Significant == "NS" ~ "",
                           sig_df$p.adj == 0 ~ "*", # This portion of the case_when should only happen with second gen p-values
                           sig_df$p.adj < 0.001 ~ "***",#the next sections will not be triggered by sgpv because of the "NS"
                           sig_df$p.adj < 0.01 ~ "**",#These only apply to FDRs
                           sig_df$p.adj < 0.05 ~ "*",
                           TRUE ~ "")
  
  max.mbh = plot.MBH[order(plot.MBH$loop_strata,plot.MBH$value,decreasing=T),]
  max.mbh = max.mbh[!duplicated(max.mbh$loop_strata),]
  max.mbh$strata = as.character(max.mbh$strata)
  max.mbh = left_join(max.mbh,sig_df)
  max.mbh$label[is.na(max.mbh$label)] = ""
  
  plot.region = plot.region +
    geom_text(data=max.mbh,aes(x = hemi,y = 1.1*value,label=label),size=8,color="black") +
    coord_cartesian(clip="off")
  
  return(plot.region)
}
