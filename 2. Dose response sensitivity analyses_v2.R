
#################################################################################################
########################### Sensitivity analyses ################################################
#################################################################################################


## Sensitivity 1 - Include all studies - multilevel


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

CVD <- read_excel("C:/Users/kowe2286/OneDrive - The University of Sydney (Staff)/SR and MA of step counts/Meta-analysis/Data extraction.xlsx", sheet = "CVD")

CVD_mortality <- CVD %>% filter(Mortality == 1 )
CVD_incidence <- CVD %>% filter(Mortality == 0 )

Diab <- read_excel("C:/Users/kowe2286/OneDrive - The University of Sydney (Staff)/SR and MA of step counts/Meta-analysis/Data extraction.xlsx", sheet = "Diabetes")


#################################################################################################
################################## 1. Multilevel model ##########################################
#################################################################################################


sensitivity <- function(ds, k_value, p_value, i2_value){
  

#Standard error         
ds$se <- (abs((log(ds$adjhr)-log(ds$lb)))/qnorm(0.975)+abs((log(ds$adjhr)-log(ds$ub)))/qnorm(0.975))/2
ds$se <- ifelse(ds$adjhr==1,NA,ds$se)

ds$dose1 <- ds$dose -2000 


#Calculate the within-study correlations  
# addS <- lapply(split(ds, ds$id), function(x) 
#  covar.logrr(y=loghr, v=I(se^2), cases=case, n=n, type=type, data=x)) 
sub <- subset(ds, !is.na(se)) 
sub <- sub[order(sub$id),]


# Define knots
knots <- quantile(sub$dose1, c(.10, .50, .90))


#Fit linear mixed model
lin_mixed <- mixmeta(loghr ~ 0 + dose1, 
               random= ~ 0 + 1 |id/Study_id,
               S=I(se^2), 
               method="ml",
               data=sub)
print(summary(lin_mixed))


#Fit non-linear mixed model
nonlin_mixed <- mixmeta(loghr ~ 0 + rcs(dose1, knots), 
                  random= ~ 0 + 1 |id/Study_id,
                  S=I(se^2), 
                  method="ml",
                  data=sub)

print(summary(nonlin_mixed))


#Compare fit
print(waldtest(b = coef(nonlin_mixed), Sigma = vcov(nonlin_mixed), Terms = 2:2))

#Predict values
res_mm_ml <- data.frame(dose1 = -2000:10000) %>% 
  cbind(spl = exp(predict(nonlin_mixed, newdata = ., ci=TRUE))) %>%
  cbind(lin = exp(predict(lin_mixed, newdata = ., ci=TRUE)))

res_mm_ml$dose <- res_mm_ml$dose1 + 2000

count<-ds %>% 
  group_by(dose) %>% 
  summarise(n = n()) %>%
  filter(dose != 0)

caption_text <- paste0("k=", k_value, "; p-nonlinearity ", p_value, "; I2 = ", i2_value, "%")

#Linear figure
lin_mixed_plot <- ggplot(res_mm_ml, aes(dose, lin.fit)) +
  geom_line() +
  geom_ribbon(aes(ymin = lin.ci.lb, ymax = lin.ci.ub), alpha = .22, fill = "black") +
  scale_y_continuous(trans = "log", breaks=c(1.0, 0.8, 0.6, 0.4, 0.2),limits = c(0.25, 1.3)) +
  scale_x_continuous(breaks=c(0, 2000, 4000, 6000, 8000, 10000, 12000),limits = c(0, 12200), expand = c(0, 0)) +
  labs(x = "Steps per day", y = "Hazard ratio", caption = caption_text) +
  theme_classic() + 
  geom_hline(yintercept = 1, linetype="dashed") + 
  geom_vline(xintercept = 2000, linetype="solid") + 
  geom_rug(data=count, aes(x=dose, alpha = n), inherit.aes = FALSE, size = 1.5, show.legend = FALSE)+
  theme(plot.caption = element_text(hjust = 0)) 



#Non linear figure
nonlin_mixed_plot <-ggplot(res_mm_ml, aes(dose, spl.fit)) +
  geom_line() +
  geom_ribbon(aes(ymin = spl.ci.lb, ymax = spl.ci.ub), alpha = .22, fill = "black") +
  scale_y_continuous(trans = "log", breaks=c(1.0, 0.8, 0.6, 0.4, 0.2),limits = c(0.25, 1.3)) +
  scale_x_continuous(breaks=c(0, 2000, 4000, 6000, 8000, 10000, 12000, 12000),limits = c(0, 12200), expand = c(0, 0)) +
  labs(x = "Steps per day", y = "Hazard ratio", caption = caption_text) +
  theme_classic() +
  geom_hline(yintercept = 1, linetype="dashed") + 
  geom_vline(xintercept = 2000, linetype="solid") + 
  geom_rug(data=count, aes(x=dose, alpha = n), inherit.aes = FALSE, size = 1.5, show.legend = FALSE) +
  theme(plot.caption = element_text(hjust = 0))


return(list(lin_mixed_plot, nonlin_mixed_plot))


}


sensitivity(ACM, 
            k_value = 18,
            p_value = "<0.001",
            i2_value = 61.0)

sensitivity(CVD_mortality, 
            k_value = 5,
            p_value = "0.52",
            i2_value = 72.9)

sensitivity(CVD_incidence, 
            k_value = 7,
            p_value = "<0.01",
            i2_value = 59.0)

sensitivity(Diab, 
            k_value = 5,
            p_value = "0.73",
            i2_value = 52.1)








