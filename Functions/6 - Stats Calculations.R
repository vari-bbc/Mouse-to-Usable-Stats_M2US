# This function takes the MBH file and tests all the strata in a comparison of interest (eg treatment)
# The returned output is a data.frame that is to be merged onto the plotting data.frame once the user has selected the strata to plot

# imports: dplyr, MASS,  emmeans, sgpv, data.table, lme4
# these are loaded in app.R

#Parameters
## MBH_data: the data frame for MBH (expected columns: value, hemi, comparison, region, strata)
## One of the following should be set
#### null.interval: the null interval for SGPV (e.g. 0.05 for +/-5% of the SD)
#### fdr_cutoff: the significance level after adjustments (eg < 5% FDR)

## Random example data.frame can be created by running create_test_file.R

#Methods, models are run stratifying by region and the 'strata' column which is typically
#the concatenation of some combination of genotype, MPI, treatment, and/or marker.
#custom options are also allowed and may be passed to this function via the 'strata' column
#If sex or batch are present they are automatically adjusted for, but this function is specifically looking for pan-sex effects
#That is, any differences have the same magnitude in males and females, just perhaps a different starting average (no interactions)


# Convenience functions
# normalize the ranks for use with sgpv
rankit = function(t) {qnorm((rank(t)-0.5)/length(t))}

# The actual modelling function
## df = data.frame
## We should maybe consider making this a limma model some day
ranked_anova <- function(df,MBH_data){
  MBH_data$strata = as.character(MBH_data$strata)
  MBH_data$comparison = as.character(MBH_data$comparison)
  
  df$comparison = as.character(df$comparison)
  
  #Conditionals to determine final model
  if(! is.null(MBH_data$sex) & ! any(grepl("sex",fixed=T,MBH_data$strata)) & ! any(grepl("sex",fixed=T,MBH_data$comparison))){
    if(!is.null(MBH_data$batch)){
      model = try(lm(rankit(value) ~ comparison + batch+sex ,data=df),silent = T)
    } else { # else for batch check, so this only has sex
      model = try(lm(rankit(value) ~ comparison + sex ,data=df),silent = T)
    } #close on the batch checks when sex is present but not in starta
  } else { # if no sex is present
    if(!is.null(MBH_data$batch)){ # check for batch
      model = try(lm(rankit(value) ~ comparison + batch ,data=df),silent = T)
    } else { # else for batch check; no batch or sex adjustments
        model = try(lm(rankit(value) ~ comparison ,data=df),silent = T)
    } #close on the batch checks
  } # closed conditional
  res=NULL
  if(! inherits(model,"try-error")){
    res.p = try(data.frame(summary(emmeans(model,pairwise~comparison)$contrasts)))
    if(! inherits(res.p,"try-error")){
      res.CI = data.frame(confint(emmeans(model,pairwise~comparison)$contrasts))
      res = left_join(res.p,res.CI)
    } else {
      res.CI = data.frame(confint(emmeans(model,pairwise~comparison)$contrasts))
      res = left_join(res.p,res.CI)
    }

  }

  res$Model = "Non-parametric ranked ANOVA"

  return(res)

}


# We may deploy this in a future version, but it's difficult to automate due to
# Checking of assumptions. This really should only b implemented if we're certain beta is a good fit
# But there maybe a strangely behaving distirbution with these types of data, so ranked anova is just safer
# beta.model <- function(df){
#   if(any(df$value == 0 | df$value == 1)){df$value = lemon(df$value,nrow(df))}
#
#   if(length(unique(df$batch))>1){
#     fit = try(betareg(value ~ compare_var*hemi + sex + batch,data=df),silent=T)
#   } else {
#     fit = try(betareg(value ~ compare_var*hemi + sex,data=df),silent=T)
#   }
#
#   if(! inherits(fit,"try-error")){
#     res.p = data.frame(summary(emmeans(fit,pairwise~compare_var|hemi,nesting=NULL)$contrasts))
#     res.CI = data.frame(confint(emmeans(fit,pairwise~compare_var|hemi,nesting=NULL)$contrasts))
#     res = left_join(res.p,res.CI)
#     res$SGPV = sgpvalue(est.lo = res$asymp.LCL,
#                         est.hi = res$asymp.UCL,
#                         null.hi = log(1.05),
#                         null.lo = log(0.95))$p.delta
#   }
#
#   res$Model = "Beta regression"
#   return(res)
#
# }


GetRegionStats <- function(MBH_data, null.interval=NULL, fdr_cutoff = NULL){
  if(! any(colnames(MBH_data) == "strata")){
        MBH_data$strata = ""
  }
  
  if(length(which(table(MBH_data$strata,MBH_data$comparison) == 0)) >0){
    return("Error: Perfect separation between comparison and strata")
  } else {
    if(!is.null(null.interval) & ! is.null(fdr_cutoff)){
      return("Error: Please populate only one of null.interval or fdr_cutoff")
    } else {
      
      MBH_data$loop_strata = paste0(MBH_data$region,MBH_data$hemi,MBH_data$strata)
      MBH_data$strata = as.character(MBH_data$strata)
      MBH_data$comparison = as.character(MBH_data$comparison)
      MBH_data$value = as.numeric(MBH_data$value)
      output.res = rbindlist(lapply(unique(MBH_data$loop_strata), function(d) data.frame(loop_strata = d ,ranked_anova(MBH_data[MBH_data$loop_strata == d,],MBH_data))),fill=T)
      
      if (length(colnames(output.res)) > 2){
        if(!is.null(null.interval)){
          #If any of these failed to converge, we make them obviously large
          output.res$lower.CL[is.na(output.res$lower.CL)] <- -99
          output.res$upper.CL[is.na(output.res$upper.CL)] <- 99
          
          #remove this if wrong
          # output.res <- output.res[!is.na(output.res$p.value),]
          
          output.res$p.adj = sgpvalue(est.lo = output.res$lower.CL,
                              est.hi = output.res$upper.CL,
                              null.hi = qnorm(.5+null.interval/100),
                              null.lo = qnorm(.5-null.interval/100))$p.delta
  
          output.res$Significant = case_when(output.res$p.adj == 0 & output.res$estimate > 0 ~ "Significant Increase",
                                                  output.res$p.adj == 0 & output.res$estimate < 0 ~ "Significant Decrease",
                                                  output.res$p.adj >0 ~ "NS",
                                                  TRUE ~ "Could not be computed")
        }
  
        if(!is.null(fdr_cutoff)){
          output.res$p.adj = p.adjust(output.res$p.value,method = "BH")
          ## This variable allows us to print the brain heatmaps with specific colors for these groups
          ## Recommended colors for both sets are "#B2182B", "#2166AC", "white", "grey 80"
          output.res$Significant = case_when(output.res$p.adj < fdr_cutoff & output.res$estimate > 0 ~ "Significant Increase",
                                                 output.res$p.adj < fdr_cutoff & output.res$estimate < 0 ~ "Significant Decrease",
                                                 output.res$p.adj > fdr_cutoff ~ "NS",
                                                 TRUE ~ "Could not be computed")
        }
        output.res = left_join(subset(output.res,select = c(loop_strata, contrast, estimate, SE, p.value, lower.CL, upper.CL, p.adj, Significant)),unique(subset(MBH_data, select= c(strata,loop_strata, region, hemi))),by="loop_strata")
        output.res = subset(output.res,select=-c(loop_strata))
        return(output.res)
      } else {
        return(F)
      }
    }
  }
}



