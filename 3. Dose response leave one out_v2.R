
#################################################################################################
################################## Leave one-out ################################################
#################################################################################################


# 1.2. Check libraries, install missing packages, update old packages, and then load required packages
libs <- c("xlsx", "dosresmeta","dplyr","ggplot2","mixmeta","readxl", "writexl", "rms","splines", "metafor", "meta", " pathviewr")
missing <- !libs %in% installed.packages()
if (any(missing)) {
  install.packages(libs[missing])
}
lapply(libs, library, character.only = TRUE)


#################################################################################################
############################## Import data ######################################################
#################################################################################################


###ACM ###
ACM <- read_excel("C:/Users/kowe2286/OneDrive - The University of Sydney (Staff)/SR and MA of step counts/Meta-analysis/Data extraction.xlsx", sheet = "ACM")
ACM_one_study <- ACM %>% filter(One_study==1)


### CVD ### 
CVD <- read_excel("C:/Users/kowe2286/OneDrive - The University of Sydney (Staff)/SR and MA of step counts/Meta-analysis/Data extraction.xlsx", sheet = "CVD")
CVDmort_one_study <- CVD %>% filter(Mortality == 1 & One_study==1)
CVDinc_one_study <- CVD %>% filter(Mortality == 0 & One_study==1)


### Diabetes ### 
Diab <- read_excel("C:/Users/kowe2286/OneDrive - The University of Sydney (Staff)/SR and MA of step counts/Meta-analysis/Data extraction.xlsx", sheet = "Diabetes")
Diab_one_study <- Diab %>% filter(One_study==1)

### Cancer ###
Cancer <- read_excel("C:/Users/kowe2286/OneDrive - The University of Sydney (Staff)/SR and MA of step counts/Meta-analysis/Data extraction.xlsx", sheet = "Cancer")
Cancer_mort_one_study <- Cancer %>% filter(Mortality == 1)
Cancer_inc_one_study <- Cancer %>% filter(Mortality == 0)

### Cancer ###
Cognition <- read_excel("C:/Users/kowe2286/OneDrive - The University of Sydney (Staff)/SR and MA of step counts/Meta-analysis/Data extraction.xlsx", sheet = "Cognition")
Cognition_one_study <- Cognition %>% filter(One_study==1)

### Depressive symptoms ###
Depressive <- read_excel("C:/Users/kowe2286/OneDrive - The University of Sydney (Staff)/SR and MA of step counts/Meta-analysis/Data extraction.xlsx", sheet = "Mental health")

### Falls ###
Falls <- read_excel("C:/Users/kowe2286/OneDrive - The University of Sydney (Staff)/SR and MA of step counts/Meta-analysis/Data extraction.xlsx", sheet = "Falls")
Falls_one_study <- Falls %>% filter(One_study==1)




#################################################################################################
########################## Function - leave one out #############################################
#################################################################################################


#Function for linear and non-linear models with leave-one-out analysis
leaveoneout <- function(ds, xlname){
  #Standard error         
  ds$se <- (abs((log(ds$adjhr)-log(ds$lb)))/qnorm(0.975)+abs((log(ds$adjhr)-log(ds$ub)))/qnorm(0.975))/2
  ds$se <- ifelse(ds$adjhr==1,NA,ds$se)
  ds$dose1 <- ds$dose -2000 
  
  # Get all unique Study_ids
  study_ids <- unique(ds$Study_id)
  
  # Create a list to store results for each iteration
  all_results <- list()
  
  # Create dataframes to store prediction results for combined plots
  combined_lin_predictions <- data.frame()
  combined_nonlin_predictions <- data.frame()
  
  # Leave-one-out analysis
  for(exclude_id in study_ids) {
    cat(paste("\n\n======= EXCLUDING STUDY_ID:", exclude_id, "=======\n\n"))
    
    # Create subset excluding current Study_id
    ds_subset <- ds[ds$Study_id != exclude_id, ]
    
    # Run analysis on this subset
    iteration_results <- run_analysis(ds_subset, exclude_id)
    
    # Store results
    all_results[[paste0("excl_", exclude_id)]] <- iteration_results
    
    # Add this iteration's predictions to combined dataframes
    iter_lin_preds <- iteration_results$predictions %>%
      select(dose, lin.fit) %>%
      mutate(excluded_study = paste0("Excl. ", exclude_id))
    
    iter_nonlin_preds <- iteration_results$predictions %>%
      select(dose, spl.fit) %>%
      mutate(excluded_study = paste0("Excl. ", exclude_id))
    
    combined_lin_predictions <- bind_rows(combined_lin_predictions, iter_lin_preds)
    combined_nonlin_predictions <- bind_rows(combined_nonlin_predictions, iter_nonlin_preds)
  }
  
  # Generate combined plots
  
  # Set up a color palette with distinct colors
  n_studies <- length(study_ids)
  color_palette <- colorRampPalette(c("black", "black", "black", "black", "black", "black", "black"))(n_studies)
  
  # Assign colors to each excluded study
  study_colors <- setNames(color_palette, paste0("Excl. ", study_ids))
  
  # Combined non-linear plot
  combined_nonlin_plot <- ggplot() +
    # Add lines for all iterations
    geom_line(data = combined_nonlin_predictions, 
              aes(x = dose, y = spl.fit, color = excluded_study), 
              size = 1) +
    scale_color_manual(values = study_colors, name = "Excluded Study") +
    scale_y_continuous(trans = "log", breaks = c(1.0, 0.8, 0.6, 0.4, 0.2), limits = c(0.25, 1.5)) +
    scale_x_continuous(breaks = c(0, 2000, 4000, 6000, 8000, 10000, 12000), limits = c(0, 12200), expand = c(0, 0)) +
    labs(x = "Steps per day", 
         y = "Hazard ratio", 
         title = "  ") +
    theme_classic() +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_vline(xintercept = 2000, linetype="solid") + 
    theme(legend.position="none")
  
  # For linear model, same approach
  combined_lin_plot <- ggplot() +
    # Add lines for all iterations
    geom_line(data = combined_lin_predictions, 
              aes(x = dose, y = lin.fit, color = excluded_study), 
              size = 1) +
    scale_color_manual(values = study_colors, name = "Excluded Study") +
    scale_y_continuous(trans = "log", breaks = c(1.0, 0.8, 0.6, 0.4, 0.2), limits = c(0.25, 1.3)) +
    scale_x_continuous(breaks = c(0, 2000, 4000, 6000, 8000, 10000, 12000), limits = c(0, 12200), expand = c(0, 0)) +
    labs(x = "Steps per day", 
         y = "Hazard ratio", 
         title = "  ") +
    theme_classic() +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_vline(xintercept = 2000, linetype="solid") + 
    theme(legend.position="none")
  
  # Save the combined plots
  dataout = "C:/Users/kowe2286/OneDrive - The University of Sydney (Staff)/FS1_Migration/kowe2286/Meta-analysis steps and mortality/Tables/"
  
  ggsave(paste0(dataout, xlname, "_combined_nonlinear.png"), combined_nonlin_plot, width = 10, height = 7)
  ggsave(paste0(dataout, xlname, "_combined_linear.png"), combined_lin_plot, width = 10, height = 7)
  
  # Return all results
  return(list(
    individual_analyses = all_results,
    combined_nonlinear_plot = combined_nonlin_plot,
    combined_linear_plot = combined_lin_plot
  ))
}

# Helper function to run a single analysis iteration
run_analysis <- function(ds, excluded_id) {
  #Calculate the within-study correlations  
  addS <- lapply(split(ds, ds$Study_id), function(x) 
    covar.logrr(y=loghr, v=I(se^2), cases=case, n=n, type=type, data=x)) 
  sub <- subset(ds, !is.na(se)) 
  sub <- sub[order(sub$Study_id),]
  
  # Define knots
  knots <- quantile(sub$dose1, c(.10, .50, .90))
  
  #Fit linear 
  lin <- mixmeta(loghr ~ 0 + dose1, 
                 random= ~ 0 + 1 |Study_id,
                 method="ml",
                 data=sub,
                 control=list(addSlist=addS))
  print(summary(lin))
  
  #Fit non-linear
  nonlin <- mixmeta(loghr ~ 0 + rcs(dose1, knots), 
                    random= ~ 0 + 1 |Study_id,
                    method="ml",
                    data=sub,
                    control=list(addSlist=addS))
  print(summary(nonlin))
  
  #Compare fit
  print(waldtest(b = coef(nonlin), Sigma = vcov(nonlin), Terms = 2:2))
  
  #Predict values
  res_mm_ml <- data.frame(dose1 = -2000:10000) %>% 
    cbind(spl = exp(predict(nonlin, newdata = ., ci=TRUE))) %>%
    cbind(lin = exp(predict(lin, newdata = ., ci=TRUE)))
  res_mm_ml$dose <- res_mm_ml$dose1 + 2000
  
  count <- ds %>% 
    group_by(dose) %>% 
    summarise(n = n()) %>%
    filter(dose != 0)
  
  #Linear figure
  lin_plot <- ggplot(res_mm_ml, aes(dose, lin.fit)) +
    geom_line() +
    geom_ribbon(aes(ymin = lin.ci.lb, ymax = lin.ci.ub), alpha = .1) +
    scale_y_continuous(trans = "log", breaks=c(1.0, 0.8, 0.6, 0.4, 0.2),limits = c(0.25, 1.15)) +
    scale_x_continuous(breaks=c(0, 2000, 4000, 6000, 8000, 10000, 12000),limits = c(0, 12200), expand = c(0, 0)) +
    labs(x = "Steps per day", y = "Hazard ratio", 
         title = paste("Linear model (excluding Study_id:", excluded_id, ")")) +
    theme_classic() + 
    geom_hline(yintercept = 1, linetype="dashed") + 
    geom_vline(xintercept = 2000, linetype="solid") + 
    geom_rug(data=count, aes(x=dose, alpha = n), inherit.aes = FALSE, size = 1.5, show.legend = FALSE) 
  
  #Non linear figure
  nonlin_plot <- ggplot(res_mm_ml, aes(dose, spl.fit)) +
    geom_line() +
    geom_ribbon(aes(ymin = spl.ci.lb, ymax = spl.ci.ub), alpha = .1) +
    scale_y_continuous(trans = "log", breaks=c(1.0, 0.8, 0.6, 0.4, 0.2),limits = c(0.25, 1.15)) +
    scale_x_continuous(breaks=c(0, 2000, 4000, 6000, 8000, 10000, 12000, 12000),limits = c(0, 12200), expand = c(0, 0)) +
    labs(x = "Steps per day", y = "Hazard ratio",
         title = paste("Non-linear model (excluding Study_id:", excluded_id, ")")) +
    theme_classic() +
    geom_hline(yintercept = 1, linetype="dashed") + 
    geom_vline(xintercept = 2000, linetype="solid") + 
    geom_rug(data=count, aes(x=dose, alpha = n), inherit.aes = FALSE, size = 1.5, show.legend = FALSE) 
  
  return(list(
    lin_plot = lin_plot, 
    nonlin_plot = nonlin_plot,
    lin_model = lin,
    nonlin_model = nonlin,
    predictions = res_mm_ml
  ))
}


################################################################

leaveoneout(ACM_one_study, "ACM_leave_one_out")
leaveoneout(CVDmort_one_study, "ACM_leave_one_out")
leaveoneout(CVDinc_one_study, "ACM_leave_one_out")
leaveoneout(Diab_one_study, "Diab_leave_one_out")
leaveoneout(Cancer_mort_one_study, "Canmort_leave_one_out")
##leaveoneout(Cancer_inc_one_study, "Caninc_leave_one_out")
## leaveoneout(Cognition_one_study, "Cog_leave_one_out")
leaveoneout(Depressive, "Dep_leave_one_out")
leaveoneout(Falls_one_study, "Falls_leave_one_out")




 


