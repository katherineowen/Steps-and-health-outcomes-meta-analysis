# 1.2. Check libraries, install missing packages, update old packages, and then load required packages
libs <- c("xlsx", "dosresmeta","dplyr","ggplot2","mixmeta","readxl", "writexl", "rms","splines", "metafor", "meta",  "ggpubr")
missing <- !libs %in% installed.packages()
if (any(missing)) {
  install.packages(libs[missing])
}
lapply(libs, library, character.only = TRUE)


library(ggpubr)

#################################################################################################
############################## Import data ######################################################
#################################################################################################


###ACM ###
ACM <- read_excel("C:/Users/kowe2286/OneDrive - The University of Sydney (Staff)/SR and MA of step counts/Meta-analysis/Data extraction.xlsx", sheet = "ACM")

#Main analysis - one study per dataset
ACM_one_study <- ACM %>% filter(One_study==1)


#Subgroup analysis - by age
ACM_younger <- ACM_one_study %>% filter(Age_category==1)
ACM_older <- ACM_one_study %>% filter(Age_category==2)

#Subgroup analysis - by device
ACM_acc <- ACM_one_study %>% filter(Device==1)
ACM_ped <- ACM_one_study %>% filter(Device==2)

#Subgroup analysis - by device location
ACM_acc_hip_waist <- ACM_one_study %>% filter(Device==1 & Device_placement_update==1)
ACM_acc_wrist <- ACM_one_study %>% filter(Device==1 & Device_placement==2)

#Sensitivity - removing those that dont meet NOS 5 and 6 (Hansen 2020)
ACM_nos <- ACM_one_study %>% filter(Study_id!=7)


### CVD ### 
CVD <- read_excel("C:/Users/kowe2286/OneDrive - The University of Sydney (Staff)/SR and MA of step counts/Meta-analysis/Data extraction.xlsx", sheet = "CVD")

#Main analysis - one study per dataset
#CVD Mortality 
CVDmort_one_study <- CVD %>% filter(Mortality == 1 & One_study==1)

#CVD Incidence
CVDinc_one_study <- CVD %>% filter(Mortality == 0 & One_study==1)

#Subgroup analysis
CVDinc_younger <- CVDinc_one_study %>% filter(Age_category==1)
CVDinc_older <- CVDinc_one_study %>% filter(Age_category==2)

#Subgroup analysis - by device
CVDinc_acc <- CVDinc_one_study %>% filter(Device==1)
CVDinc_ped <- CVDinc_one_study %>% filter(Device==2)


### Diabetes ### 
Diab <- read_excel("C:/Users/kowe2286/OneDrive - The University of Sydney (Staff)/SR and MA of step counts/Meta-analysis/Data extraction.xlsx", sheet = "Diabetes")

#Main analysis - one study per dataset
Diab_one_study <- Diab %>% filter(One_study==1)

#Subgroup analysis - by age
Diab_younger <- Diab_one_study %>% filter(Age_category==1)
Diab_older <- Diab_one_study %>% filter(Age_category==2)


### Cancer ###
Cancer <- read_excel("C:/Users/kowe2286/OneDrive - The University of Sydney (Staff)/SR and MA of step counts/Meta-analysis/Data extraction.xlsx", sheet = "Cancer")

#Main analysis - one study per dataset
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


##############################################################################################


#Function for linear and non-linear models
dosres <- function(ds, xlname, title, k_value, n_value, p_value, i2_value, lin_plot_name = NULL, plot_name = NULL){

  

#Standard error         
ds$se <- (abs((log(ds$adjhr)-log(ds$lb)))/qnorm(0.975)+abs((log(ds$adjhr)-log(ds$ub)))/qnorm(0.975))/2
ds$se <- ifelse(ds$adjhr==1,NA,ds$se)

ds$dose1 <- ds$dose -2000 


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


#Write summary table to excel (by 1000 steps)
excel <- data.frame(dose1 <- seq(-2000, 10000, 1000)) %>% 
  cbind(spl = exp(predict(nonlin, newdata = ., ci=TRUE))) %>%
  cbind(lin = exp(predict(lin, newdata = ., ci=TRUE)))

excel$dose <- excel$dose1 + 2000


excel$HR_nonlin <- paste(round(excel$spl.fit, 2), " (", round(excel$spl.ci.lb, 2), ", ", round(excel$spl.ci.ub, 2), ")", sep="")
excel$HR_lin <- paste(round(excel$lin.fit, 2), " (", round(excel$lin.ci.lb, 2), ", ", round(excel$lin.ci.ub, 2), ")", sep="")


dataout = "C:/Users/kowe2286/OneDrive - The University of Sydney (Staff)/FS1_Migration/kowe2286/Meta-analysis steps and mortality/Tables/"
filename = glue::glue("{xlname}.xlsx")  
write_xlsx(excel, paste0(dataout, filename))



count<-ds %>% 
  group_by(dose) %>% 
  summarise(n = n()) %>%
  filter(dose != 0)


caption_text <- paste0("k=", k_value, "; n=", n_value, "; p-nonlinearity ", p_value, "; I2 = ", i2_value, "%")


#Linear figure
lin_plot <- ggplot(res_mm_ml, aes(dose, lin.fit)) +
  geom_line() +
  geom_ribbon(aes(ymin = lin.ci.lb, ymax = lin.ci.ub), alpha = .22, fill = "black") +
  scale_y_continuous(trans = "log", breaks=c(1.2, 1.0, 0.8, 0.6, 0.4, 0.2),limits = c(0.25, 1.4)) +
  scale_x_continuous(breaks=c(0, 2000, 4000, 6000, 8000, 10000, 12000),limits = c(0, 12400), expand = c(0, 0)) +
  labs(x = "Steps per day", y = "Hazard ratio", title = title, caption = caption_text) +
  theme_classic() + 
  geom_hline(yintercept = 1, linetype="dashed") + 
  geom_vline(xintercept = 2000, linetype="solid") + 
  geom_rug(data=count, aes(x=dose, alpha = n), inherit.aes = FALSE, size = 1.5, show.legend = FALSE) +
  theme(plot.caption = element_text(hjust = 0))

# Assign name to plot object if provided
if (!is.null(lin_plot_name)) {
  assign(lin_plot_name, lin_plot, envir = .GlobalEnv)
}



#Non linear figure
nonlin_plot <-ggplot(res_mm_ml, aes(dose, spl.fit)) +
  geom_line() +
  geom_ribbon(aes(ymin = spl.ci.lb, ymax = spl.ci.ub), alpha = .22, fill = "black") +
  scale_y_continuous(trans = "log", breaks=c(1.2, 1.0, 0.8, 0.6, 0.4, 0.2),limits = c(0.25, 1.4)) +
  scale_x_continuous(breaks=c(0, 2000, 4000, 6000, 8000, 10000, 12000, 12000),limits = c(0, 12400), expand = c(0, 0)) +
  labs(x = "Steps per day", y = "Hazard ratio", title = title, caption = caption_text) +
  theme_classic() +
  geom_hline(yintercept = 1, linetype="dashed") + 
  geom_vline(xintercept = 2000, linetype="solid") + 
  geom_rug(data=count, aes(x=dose, alpha = n), inherit.aes = FALSE, size = 1.5, show.legend = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


# Assign name to plot object if provided
if (!is.null(plot_name)) {
  assign(plot_name, nonlin_plot, envir = .GlobalEnv)
}


#optimal steps
res_mm_ml1 <- res_mm_ml %>%
  select(dose,spl.fit)

optimal_number <- find_curve_elbow(res_mm_ml1, export_type = "row_num", plot_curve = TRUE)

return(list(optimal_number, lin_plot, nonlin_plot))


}


##########################################################################
#################### All-cause mortality #################################
##########################################################################

dosres(ACM_one_study, "ACM_overall", 
       title = "a) All cause mortality",
       k_value = 14,
       n_value = "123,504",
       p_value = "<0.001",
       i2_value = 41.7,
       lin_plot_name = "acm_linear",
       plot_name = "acm_nonlinear")

plots_af <- ggarrange(acm_nonlinear, cvd_ind_nonlinear, 
                      cvd_mort_linear, can_inc_linear,
                      can_mort_linear, diab_linear, 
                      ncol=2, nrow=3)


ggsave("C:/Users/kowe2286/OneDrive - The University of Sydney (Staff)/FS1_Migration/kowe2286/Meta-analysis steps and mortality/plots_af.png", 
       plots_af, dpi = 300,
       width = 8, height = 12)

plots_gj <- ggarrange(cog_nonlinear, dep_linear,
                      fall_nonlinear, 
                      ncol=2, nrow=2)


ggsave("C:/Users/kowe2286/OneDrive - The University of Sydney (Staff)/FS1_Migration/kowe2286/Meta-analysis steps and mortality/plots_gj.png", 
       plots_gj, dpi = 300,
       width = 8, height = 8)



#by age
dosres(ACM_younger, "ACM_young", 
       title = "b) Younger adults",
       k_value = 8,
       n_value = "102,420",
       p_value = "<0.001",
       i2_value = 49.5,
       lin_plot_name = "acm_young_linear",
       plot_name = "acm_young_nonlinear")

dosres(ACM_older, "ACM_old", 
       title = "a) Older adults",
       k_value = 6,
       n_value = "21,084",
       p_value = "0.051",
       i2_value = 38.0,
       lin_plot_name = "acm_old_linear",
       plot_name = "acm_old_nonlinear")

plots_acm_age <- ggarrange(acm_old_linear, acm_young_nonlinear, ncol=2, nrow=1)
ggsave("C:/Users/kowe2286/OneDrive - The University of Sydney (Staff)/FS1_Migration/kowe2286/Meta-analysis steps and mortality/plots_acm_age.png", 
       plots_acm_age, dpi = 300,
       width = 8, height = 4)


#By device
dosres(ACM_acc, "ACM_acc", 
       title = "a) Accelerometer",
       k_value = 10,
       n_value = "117,124",
       p_value = "<0.001",
       i2_value = 51.2,
       lin_plot_name = "acm_acc_linear",
       plot_name = "acm_acc_nonlinear")

dosres(ACM_ped, "ACM_ped", 
       title = "b) Pedometer",
       k_value = 4,
       n_value = "6,380",
       p_value = "0.68",
       i2_value = 1.8,
       lin_plot_name = "acm_ped_linear",
       plot_name = "acm_ped_nonlinear")


plots_acm_device <- ggarrange(acm_acc_nonlinear, acm_ped_linear, ncol=2, nrow=1)
ggsave("C:/Users/kowe2286/OneDrive - The University of Sydney (Staff)/FS1_Migration/kowe2286/Meta-analysis steps and mortality/plots_acm_device.png", 
       plots_acm_device, dpi = 300,
       width = 8, height = 4)


dosres(ACM_nos, "ACM_nos", 
       title = " ",
       k_value = 13,
       n_value = "121,321",
       p_value = "<0.001",
       i2_value = 45.9,
       lin_plot_name = "acm_nos_linear",
       plot_name = "acm_nos_nonlinear")



##########################################################################
################# CVD incidence and mortality ############################
##########################################################################

dosres(CVDmort_one_study, "CVD_mortality_overall", 
       title = "c) CVD mortality",
       k_value = 3,
       n_value = "85,538",
       p_value = "0.29",
       i2_value = 85.1,
       lin_plot_name = "cvd_mort_linear",
       plot_name = "cvd_mort_nonlinear")


dosres(CVDinc_one_study, "CVD_incidence_overall", 
       title = "b) CVD incidence",
       k_value = 6,
       n_value = "103,291",
       p_value = "0.01",
       i2_value = 40.2,
       lin_plot_name = "cvd_ind_linear",
       plot_name = "cvd_ind_nonlinear")


#By age
dosres(CVDinc_older, "CVD_inc_older", 
       title = "a) Older adults",
       k_value = 3,
       n_value = "29,683",
       p_value = "0.08",
       i2_value = 23.2,
       lin_plot_name = "cvd_old_linear",
       plot_name = "cvd_old_nonlinear")

dosres(CVDinc_younger, "CVD_inc_younger", 
       title = "b) Younger adults",
       k_value = 3,
       n_value = "81,666",
       p_value = "<0.001",
       i2_value = 0.0,
       lin_plot_name = "cvd_young_linear",
       plot_name = "cvd_young_nonlinear")
       

plots_cvd_age <- ggarrange(cvd_old_linear, cvd_young_nonlinear, ncol=2, nrow=1)
ggsave("C:/Users/kowe2286/OneDrive - The University of Sydney (Staff)/FS1_Migration/kowe2286/Meta-analysis steps and mortality/plots_cvd_age.png", 
       plots_cvd_age, dpi = 300,
       width = 8, height = 4)


#By device
dosres(CVDinc_acc, "CVD_inc_acc",
       title = "a) Accelerometer",
       k_value = 4,
       n_value = "108,177",
       p_value = "<0.01",
       i2_value = 51.6,
       lin_plot_name = "cvd_acc_linear",
       plot_name = "cvd_acc_nonlinear")

dosres(CVDinc_ped, "CVD_inc_ped",
       title = "b) Pedometer",
       k_value = 2,
       n_value = "3,172",
       p_value = "0.69",
       i2_value = 0.0,
       lin_plot_name = "cvd_ped_linear",
       plot_name = "cvd_ped_nonlinear")

plots_cvd_device <- ggarrange(cvd_acc_nonlinear, cvd_ped_linear, ncol=2, nrow=1)
ggsave("C:/Users/kowe2286/OneDrive - The University of Sydney (Staff)/FS1_Migration/kowe2286/Meta-analysis steps and mortality/plots_cvd_device.png", 
       plots_cvd_device, dpi = 300,
       width = 8, height = 4)


##########################################################################
############################ Cancer ######################################
##########################################################################


dosres(Cancer_inc_one_study, "Cancer_incidnce_overall", 
       title = "d) Cancer incidence",
       k_value = 2,
       n_value = "100,505",
       p_value = "0.80",
       i2_value = 73.7,
       lin_plot_name = "can_inc_linear",
       plot_name = "can_inc_nonlinear")



dosres(Cancer_mort_one_study, "Cancer_mortality_overall", 
       title = "e) Cancer mortality",
       k_value = 3,
       n_value = "105,660",
       p_value = "0.23",
       i2_value = 69.0,
       lin_plot_name = "can_mort_linear",
       plot_name = "can_mort_nonlinear")



##########################################################################
############################ Diabetes ####################################
##########################################################################


dosres(Diab_one_study, "Diab_overall", 
       title = "f) Type 2 diabetes",
       k_value = 4,
       n_value = "49,609",
       p_value = "0.62",
       i2_value = 48.5,
       lin_plot_name = "diab_linear",
       plot_name = "diab_nonlinear")


#By age
dosres(Diab_older, "Diab_older", 
       title = "a) Older adults",
       k_value = 2,
       n_value = "7,893",
       p_value = "0.09",
       i2_value = 2.2,
       lin_plot_name = "diab_old_linear",
       plot_name = "diab_old_nonlinear")

dosres(Diab_younger, "Diab_younger", 
       title = "b) Younger adults",
       k_value = 2,
       n_value = "36,039",
       p_value = "0.09",
       i2_value = 66.8,
       lin_plot_name = "diab_young_linear",
       plot_name = "diab_young_nonlinear")


plots_diab_age <- ggarrange(diab_old_linear, diab_young_linear, ncol=2, nrow=1)
ggsave("C:/Users/kowe2286/OneDrive - The University of Sydney (Staff)/FS1_Migration/kowe2286/Meta-analysis steps and mortality/plots_diab_age.png", 
       plots_diab_age, dpi = 300,
       width = 8, height = 4)



##########################################################################
############################ Cognition ###################################
##########################################################################

dosres(Cognition_one_study, "Cognition_overall", 
       title = "g) Dementia",
       k_value = 2,
       n_value = "79,973",
       p_value = "<0.001",
       i2_value = 12.1,
       lin_plot_name = "cog_linear",
       plot_name = "cog_nonlinear")

##########################################################################
################## Depressive symptoms ###################################
##########################################################################

dosres(Depressive, "Depressive_symptoms", 
       title = "h) Depressive symptoms",
       k_value = 3,
       n_value = "77,565",
       p_value = "0.39",
       i2_value = 36.2,
       lin_plot_name = "dep_linear",
       plot_name = "dep_nonlinear")


##########################################################################
############################ Falls #######################################
##########################################################################

dosres(Falls_one_study, "Falls", 
       title = "i) Falls",
       k_value = 4,
       n_value = "38,882",
       p_value = "<0.01",
       i2_value = 49.3,
       lin_plot_name = "fall_linear",
       plot_name = "fall_nonlinear")






#Publication bias
ACM_one_study$se <- (abs((log(ACM_one_study$adjhr)-log(ACM_one_study$lb)))/qnorm(0.975)+abs((log(ACM_one_study$adjhr)-log(ACM_one_study$ub)))/qnorm(0.975))/2
ACM_one_study$se <- ifelse(ACM_one_study$adjhr==1,NA,ACM_one_study$se)

sub <- subset(ACM_one_study, !is.na(se)) 

regtest(x=sub$adjhr, sei=sub$se)
funnel(sub$adjhr, sub$se, refline=1, xlab="Effect size (HR)", ylab = "Standard error", contour.levels=NULL)







