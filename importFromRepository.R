library(dplyr)

user <- 'Tom'
if (user=='Tom'){
  dataSource <- "C:/git_repos/Anopheles-model/downloads/"  
  scriptPath  = 'C:/git_repos/Anopheles-model/'
}

# READ THE DATA 
repo <- read.csv(file=paste0(dataSource,"Bionomics Americas.csv"))

# DEFINE:
#           taxon- the data will be averaged across subgroups for each level of taxon
#                     by default the groups are species   
#           subgroup- the data within each level of taxon will be aggregated by subgroup
#                     by default the subgroups are sites
 
repo$taxon <- repo$species
repo$subgroup <- repo$site

# create data frames for the parameter values disaggregated by taxon and by subgroup within taxon
dPar <- as.data.frame(table(repo$taxon,repo$subgroup))
dPar <- dPar[dPar$Freq>0,1:2]
names(dPar) <- c("taxon", "subgroup")
taxon <- as.data.frame(table(repo$taxon))[,1:2]
names(taxon) <- c("taxon", "input records")


# FUNCTIONS
weightedMean <- function(x1,w1,x2=0,w2=0,x3=0,w3=0){
  x1[which(is.na(w1))] <- 0
  w1[which(is.na(w1))] <- 0
  x2[which(is.na(w2))] <- 0
  w2[which(is.na(w2))] <- 0
  x3[which(is.na(w3))] <- 0
  w3[which(is.na(w3))] <- 0
  return((x1*w1 + x2*w2 +x3*w3)/(w1+w2+w3))
}

# Arithmetic mean from multiple records in the same subgroup 

summary_mean <-function(df,var){
  df$var <- df[[var]]  
  summary <- df %>% 
                  group_by(taxon,subgroup) %>% 
                  summarise(varmean = mean(var, na.rm = TRUE))
return(summary)}

# Arithmetic mean from multiple records in the same taxon 
taxon_mean <-function(df,var){
  df$var <- df[[var]]  
  summary <- df %>% 
    group_by(taxon) %>% 
    summarise(varmean = mean(var, na.rm = TRUE))
  return(summary)}

# Calculate ratio from multiple records in the same subgroup with totals and n values 
summary_ratio <-function(df,total_var,n_var){
  df$total_var <- df[[total_var]]
  df$total_var[is.na(df$total_var)] <- 0
  df$n_var <- df[[n_var]]
  df$n_var[is.na(df$n_var)] <- 0
  summary <- as.data.frame(df %>%  
                             group_by(taxon,subgroup) %>% 
                             summarise(sum_total_var = sum(total_var), 
                                       sum_n_var = sum(n_var)))
  summary$ratio <- summary$sum_n_var/summary$sum_total_var
return(summary)}

# Compute the inputs used by the R package separately for each taxon and subgroup

# VECTOR BIOLOGY

# Parous rate             

dPar$parity <- summary_ratio(repo,'parity_total','parity_n')$ratio
dPar$gonotrophic_cycle <- summary_mean(repo,'gonotrophic_cycle_days')$varmean
dPar$DailySurvival <- summary_mean(repo,'daily_survival_rate_percent')$varmean/100 
Survival_per_cycle_data <- c()
Sac_rate_data <- c()

# HUMAN BITING AND RESTING
Indoor_Outdoor_biting_data <- which(repo$indoor_hbr_sampling !="" & repo$outdoor_hbr_sampling !="") 
Indoor_resting_total_biting <- which(repo$indoor_resting_sampling !="" 
                                     & repo$indoor_hbr_sampling !="" 
                                     & repo$outdoor_hbr_sampling !="")
Indoor_resting_biting_data <- which(repo$indoor_resting_sampling !="" 
                                    & repo$indoor_hbr_sampling !="")
Outdoor_human_animal_data  <- which(repo$outdoor_hbr_sampling !="" 
                                    & repo$outdoor_hbr_sampling !="")

# IB/(IB+OB) 
repo$proportion_indoor_biting <- repo$indoor_hbr /(repo$indoor_hbr + repo$outdoor_hbr)
dPar$proportion_indoor_biting <- summary_mean(repo,'proportion_indoor_biting')$varmean
taxon$proportion_indoor_biting <- taxon_mean(dPar,'proportion_indoor_biting')$varmean

# IR/(totalB)
repo$indoor_resting_by_all_biting <- repo$indoor_total/(repo$indoor_hbr + repo$outdoor_hbr)
dPar$indoor_resting_by_all_biting <- summary_mean(repo,'indoor_resting_by_all_biting')$varmean
taxon$indoor_resting_by_all_biting <- taxon_mean(dPar,'indoor_resting_by_all_biting')$varmean

# IR/IB 
repo$indoor_resting_by_biting <- repo$indoor_total/repo$indoor_hbr
dPar$indoor_resting_by_biting <- summary_mean(repo,'indoor_resting_by_biting')$varmean
taxon$indoor_resting_by_biting <- taxon_mean(dPar,'indoor_resting_by_biting')$varmean









HBI         
       
 
             
OB/OA        
Pr(IR)           
sac rate 
Survival per cycle 






# VECTOR INFECTION RATE

dPar$sr_dissection <- summary_ratio(repo,'sr_dissection_total','sr_dissection_n')$ratio
dPar$sr_dissection_prop <- summary_mean(repo,'sr_dissection_percent')$varmean/100
dPar$sr_csp <- summary_ratio(repo,'sr_csp_total','sr_csp_n')$ratio
dPar$sr_csp_prop <- summary_mean(repo,'sr_csp_percent')$varmean/100
dPar$sporozoite_rate_prop <- rowMeans(subset(dPar, select = c(sr_dissection_prop,sr_csp_prop)), na.rm = TRUE)
dPar$sporozoite_rate <- rowMeans(subset(dPar, select = c(sr_dissection,sr_csp)), na.rm = TRUE)
dPar$sporozoite_rate[is.na(dPar$sporozoite_rate) || (dPar$sporozoite_rate > 1)] <- dPar$sporozoite_rate_prop
  
taxon$sporozoite_rate <- taxon_mean(dPar,'sporozoite_rate')$varmean







sr_data <- which(!is.na(repo$sr_dissection_percent) 
                | !is.na(repo$sr_dissection_total_n)
                | !is.na(repo$sr_csp_percent)
                | !is.na(repo$sr_csp_n))
repo$sporozoite_rate <- with(repo,weightedMean(sr_dissection_percent,sr_dissection_total,
                                          sr_csp_percent,sr_csp_total))/100

preference3 <- summary_mean(repo,'sr_dissection_percent')$varmean/100
preference4 <- summary_mean(repo,'sr_csp_percent')$varmean/100


repo$oocyst_rate <- repo$oocyst_percent/100


# HOST PREFERENCE DATA


repo$HBI <- with(repo,weightedMean(indoor_host,indoor_host_total,
                                   outdoor_host,outdoor_host_total, 
                                   combined_host,combined_host_total))
                 
# Select which records to use 

control <- which(repo$control_type != "")


# AVERAGE OVER ALL RECORDS FOR EACH TAXON

taxon$DailySurvival <- taxon_mean(dPar,'DailySurvival')$varmean 
taxon$parity <- taxon_mean(dPar,'parity')$varmean
taxon$gonotrophic_cycle <- taxon_mean(dPar,'gonotrophic_cycle')$varmean 

  
        


# read species specific bionomic parameters except biting rhythms
url1 <- 'http://raw.githubusercontent.com/SwissTPH/Anopheles-model/master/raw_data/bionomicsData.xlsx'  
bionomicsData <- as.data.frame(read.xls(url1,sheet=1))  