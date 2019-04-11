library(dplyr)


user <- 'Tom'
if (user=='Tom'){
  dataSource <- "C:/git_repos/Anopheles-model/downloads/"  
  scriptPath  = 'C:/git_repos/Anopheles-model/'
}

# READ THE DATA 
africa <- read.csv(file=paste0(dataSource,"Bionomics Africa.csv"))
america <- read.csv(file=paste0(dataSource,"Bionomics Americas.csv"))
asia <- read.csv(file=paste0(dataSource,"Bionomics Asia-Pacific.csv"))
repo <- rbind(africa,america,asia)

# DEFINE:
#           taxon- the data will be averaged across subgroups for each level of taxon
#                     by default the groups are species   
#           subgroup- the data within each level of taxon will be aggregated by subgroup
#                     by default the subgroups are sites

# STANDARDISE THE ORTHOGRAPHY OF THE TAXA
# map onto the taxonomy file to standardise spelling of species and species complexes
# species complexes are mapped onto the nominate species
species_mapping <- as.data.frame(read.csv(file=paste0(dataSource,"species_mapping.csv"),sep=";"))  
taxonomy <- as.data.frame(read.csv(file=paste0(dataSource,"taxonomy.csv"),sep=";"))
taxonomy$taxon <- as.character(taxonomy$taxon)

standardise_taxon <- function(repo=repo){
  repo$Taxon_no <- species_mapping$Taxon_no[match(repo$species, species_mapping$species)]
  repo$taxon <- taxonomy$taxon[match(repo$Taxon_no, taxonomy$Taxon_no)]
return(repo)}

repo <- standardise_taxon(repo)

# subgrouping within species can be changed here
repo$subgroup <- repo$site

# create data frame for the parameter values disaggregated by taxon and by subgroup within taxon
dPar <- as.data.frame(table(repo$taxon,repo$subgroup))
dPar <- dPar[dPar$Freq>0,1:2]
names(dPar) <- c("taxon", "subgroup")


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
  df$counted <- ifelse(is.na(df$var),0,1) 
  summary <- df %>% 
                  group_by(taxon,subgroup) %>% 
                  summarise(varmean = mean(var, na.rm = TRUE),
                            obs= sum(counted))
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
  df$n_var <- df[[n_var]]
  df$n_var[is.na(df$n_var) | is.na(df$total_var)] <- NA
  df$total_var[is.na(df$n_var)] <- NA
  df$counted <- ifelse(is.na(df$total_var),0,1) 
  summary <- as.data.frame(df %>%  
                             group_by(taxon,subgroup) %>% 
                             summarise(sum_total_var = sum(total_var, na.rm = TRUE), 
                                       sum_n_var = sum(n_var, na.rm = TRUE),
                                       obs= sum(counted)))
  summary$ratio <- summary$sum_n_var/summary$sum_total_var
return(summary)}

# Calculate a weighted average of proportions given as r/n and of values expressed as percentages
mean_proportion <- function(df,total_var,n_var,percent_var){
  df$total_var <- df[[total_var]]
  df$n_var <- df[[n_var]]
  df$percent_var <- df[[percent_var]]
  df1 <- summary_ratio(df,total_var,n_var)
  df$proportion_var <- ifelse(is.na(df$total_var),NA,df$percent_var/100) 
  df2 <- summary_mean(df,'proportion_var')
  mean_proportion <- (df1$ratio * df1$obs +df2$varmean * df2$obs)/(df1$obs + df2$obs)
return(mean_proportion)}  

# Compute the inputs used by the R package separately for each taxon and subgroup

# VECTOR BIOLOGY

# Parous rate   

# clean up data with very small totals
repo$parity_total[repo$parity_total<6] <- NA
repo$parity_n[is.na(repo$parity_n)] <- with(repo, parity_percent[is.na(parity_n)]*parity_total[is.na(parity_n)]/100)
dPar$parity <- summary_ratio(repo,'parity_total','parity_n')$ratio
dPar$DailySurvival <- summary_mean(repo,'daily_survival_rate_percent')$varmean/100 

# HUMAN BITING AND RESTING

# transform both ratios and percentages into proportions
repo$proportion_indoor_biting <- with(repo,ifelse(indoor_outdoor_biting_units=="I:O",
                                    indoor_biting/(1+indoor_biting),indoor_biting/100))
dPar$proportion_indoor_biting <- summary_mean(repo,'proportion_indoor_biting')$varmean

# Use the ratio of resting to biting data to estimate endophily, if there are no direct estimates
# indoor resting (collection per man hour) by total HBR
repo$indoor_total[repo$indoor_resting_sampling != 'HRI'] <- NA
repo$indoor_resting_by_all_biting <- repo$indoor_total/(repo$indoor_hbr + repo$outdoor_hbr)
dPar$indoor_resting_by_all_biting <- summary_mean(repo,'indoor_resting_by_all_biting')$varmean

# ratio of indoor resting (collection per man hour) by indoor HBR
repo$indoor_resting_by_biting <- repo$indoor_total/repo$indoor_hbr
repo$indoor_resting_by_biting[repo$resting_unit !='per man hour'] <- NA
dPar$indoor_resting_by_biting <- summary_mean(repo,'indoor_resting_by_biting')$varmean

# proportion indoor resting estimated directly (usually from exit traps)
repo$proportion_indoor_resting <- repo$indoor_total/(repo$indoor_total + repo$outdoor_total)
dPar$proportion_indoor_resting <- summary_mean(repo,'proportion_indoor_resting')$varmean

# DURATION OF GONOTROPHIC CYCLE
# 
# Exclude estimates of gonotrophic cycle based only on fed:gravid ratio 
repo$gonotrophic_cycle_days<- with(repo,ifelse(pubmed_id == 16739419 | pubmed_id == 20482824,NA,gonotrophic_cycle_days))
# Calculate gravid_fed as a proportion   
repo$gravid_fed <- repo$indoor_gravid/(repo$indoor_fed+repo$indoor_gravid)

#identify papers from the MAP repository using sac rates: these will be merged in with the additional data
sac_data <- read.csv(file=paste0(dataSource,"sac_rate.csv"),sep=";")
sac_data <- standardise_taxon(sac_data)
dPar$gonotrophic_cycle <- summary_mean(repo,'gonotrophic_cycle_days')$varmean
dPar$gravid_fed <- summary_mean(repo,'gravid_fed')$varmean


# VECTOR INFECTION RATE

dPar$sr_dissection <- summary_ratio(repo,'sr_dissection_total','sr_dissection_n')$ratio
dPar$sr_dissection_prop <- summary_mean(repo,'sr_dissection_percent')$varmean/100
dPar$sr_csp <- summary_ratio(repo,'sr_csp_total','sr_csp_n')$ratio
dPar$sr_csp_prop <- summary_mean(repo,'sr_csp_percent')$varmean/100
dPar$sporozoite_rate_prop <- rowMeans(subset(dPar, select = c(sr_dissection_prop,sr_csp_prop)), na.rm = TRUE)
dPar$sporozoite_rate <- rowMeans(subset(dPar, select = c(sr_dissection,sr_csp)), na.rm = TRUE)
# use values derived from n and total if available, otherwise the summary values
dPar$sporozoite_rate[is.na(dPar$sporozoite_rate) || (dPar$sporozoite_rate > 1)] <- dPar$sporozoite_rate_prop[is.na(dPar$sporozoite_rate) || (dPar$sporozoite_rate > 1)]

dPar$oocyst_rate <- summary_ratio(repo,'oocyst_total','oocyst_n')$ratio
dPar$oocyst_prop <- summary_mean(repo,'oocyst_percent')$varmean/100
# use values derived from n and total if available, otherwise the summary values
dPar$oocyst_rate[is.na(dPar$oocyst_rate)] <- dPar$oocyst_prop[is.na(dPar$oocyst_rate)]

# HOST PREFERENCE DATA

# Anthropophilous index from human vs buffalo (or other animal) collections treated as a proportion

repo$AIprop <- with(repo,ifelse(host_unit=='AI',0.01*outdoor_host,NA))
dPar$AIprop <- summary_mean(repo,'AIprop')$varmean

# Human blood index

repo$HBI_total <- rowSums(subset(repo, select = c(indoor_host_total,outdoor_host_total,combined_host_total)), na.rm = TRUE)
repo$HBI_n <- rowSums(subset(repo, select = c(indoor_host_n,outdoor_host_n,combined_host_n)), na.rm = TRUE)
repo$HBI_perc <- rowMeans(subset(repo, select = c(indoor_host,outdoor_host,combined_host)), na.rm = TRUE)
# set data to missing if the recorded values are not HBI
repo$HBI_total[repo$host_unit[substr(repo$host_unit,1,3) != 'HBI']] <- NA 
repo$HBI_total[repo$HBI_total==0] <- NA 
repo$HBI_perc[repo$host_unit[substr(repo$host_unit,1,3) != 'HBI']] <- NA 
dPar$Chi <- summary_ratio(repo,'HBI_total','HBI_n')$ratio
# use values derived from n and total if available, otherwise the summary values
dPar$Chi[is.na(dPar$Chi)] <- (summary_mean(repo,'HBI_perc')$varmean/100)[is.na(dPar$Chi)] 

# MERGE IN ADDITIONAL DATA 
added <- read.csv(file=paste0(dataSource,"additional_bionomics.csv"),sep=";")
added <- standardise_taxon(added)
added$subgroup <- paste0(added$country,added$site)

library(reshape2)
suppl <- dcast(added, taxon+subgroup~variable, value.var='value', mean)
# concatenate the data frames while suppressing the message about coercing the factors to characters
#at this point any duplicated data should be removed
suppressWarnings(dPar_all <- bind_rows(dPar,suppl,sac_data))

# AVERAGE OVER ALL RECORDS FOR EACH TAXON

taxon <- as.data.frame(table(dPar_all$taxon))[,1:2]
names(taxon) <- c("taxon", "input records")
taxon$vector <- taxonomy$Vector[match(taxon$taxon, taxonomy$taxon)]
taxon$DailySurvival <- taxon_mean(dPar_all,'DailySurvival')$varmean 
taxon$parity <- taxon_mean(dPar_all,'parity')$varmean
taxon$sporozoite_rate <- taxon_mean(dPar_all,'sporozoite_rate')$varmean
taxon$oocyst_rate <- taxon_mean(dPar_all,'oocyst_rate')$varmean
taxon$Endophagy <- taxon_mean(dPar_all,'proportion_indoor_biting')$varmean
taxon$indoor_resting_by_biting <- taxon_mean(dPar_all,'indoor_resting_by_biting')$varmean
taxon$indoor_resting_by_all_biting <- taxon_mean(dPar_all,'indoor_resting_by_all_biting')$varmean
taxon$proportion_indoor_resting <- taxon_mean(dPar_all,'proportion_indoor_resting')$varmean
taxon$Chi <- taxon_mean(dPar_all,'Chi')$varmean
taxon$AIprop <- taxon_mean(dPar_all,'AIprop')$varmean
taxon$sac_rate <- taxon_mean(dPar_all,'sac_rate')$varmean
taxon$gravid_fed <- taxon_mean(dPar_all,'gravid_fed')$varmean
taxon$theta_f <- taxon_mean(dPar_all,'gonotrophic_cycle')$varmean #Total duration of the gonotrophic cycle 
taxon$tau <- 1/taxon_mean(dPar_all,'gravid_fed')$varmean  #Duration of the resting period
# If both the sac rate and the fed:gravid ratio are available then the duration of the gonotrophic cycle can be estimated 
# using the Charlwood (2016) formula at the level of the taxon.
taxon$theta_f[is.na(taxon$theta_f)] <- with(taxon,tau[is.na(theta_f)]+(1/sac_rate[is.na(theta_f)])-1) 


##################################################################################     
# Add in the biting rhythm

url2 <- 'http://raw.githubusercontent.com/SwissTPH/Anopheles-model/master/raw_data/bitingRhythm.xlsx'  

# read species specific biting rhythms

bitingRhythms <- read.csv(file=paste0(dataSource,"biting_rhythm.csv"),sep=";") 

source ('calculate_Pii.r')

# AVERAGING OF BITING RHYTHMS FOR ESTIMATING pii 

average_biting_rhythms <- function(species, Location=NULL,	StudyID=NULL,	country=NULL,	Site=NULL,	Season=NULL){
  rspecies <- which(biting_rhythm$species == species)
  rLocation <- rStudyID <- rCountry <- rSite <- rSeason <- rspecies
  if (!is.null(Location)) {rLocation <- which(biting_rhythm$Location %in% Location)} 
  if (!is.null(StudyID)) {rStudyID <- which(biting_rhythm$StudyID %in% StudyID)}
  if (!is.null(Country)) {rCountry <- which(biting_rhythm$country %in% country)}
  if (!is.null(Site)) {rSite <- which(biting_rhythm$Site %in% Site)}
  if (!is.null(Season)) {rSeason <- which(biting_rhythm$Season %in% Season)}
  rrhythms <- biting_rhythm[Reduce(intersect,list(rspecies,rLocation,rStudyID,rCountry,rSite,rSeason)),8:20]  
  average_rhythm <- NULL
  if (length(rrhythms[1]) == 1) average_rhythm <- rrhythms
  if (length(rrhythms[1]) > 1) average_rhythm <- colMeans(as.numeric(rrhythms), na.rm=TRUE)
  return(average_rhythm)  
}


# for the moment assume the same rhythm indoor and outdoor
# change this to use separate indoor outdoor data where available
taxon$Indoor <- taxon$Endophagy * average_biting_rhythms(species=vectorName)
taxon$Outdoor <- (1 - taxon$Endophagy) * average_biting_rhythms(species=vectorName)


############### FROM HERE ON NEEDS AMENDMENT

# POPULATE PROPORTION INDOOR RESTING 

# Select endophily data
bionomicsData$variable[bionomicsData$variable == 'IR/IB'] <- 'IR/(totalB)'
EndophilyData <- bionomicsData[ which(bionomicsData$variable=='IR/(totalB)'), ]
table(EndophilyData$species,EndophilyData$Site)

Endo <- with(dPar,taxon,)
Endo$logendo <- log(dPar$proportion_indoor_resting+0.001)
library(nlme)
model2 <- with(dPar,lme(logendo~0 + taxon,random=~1|subgroup,method="REML"))
summary(model2)
logEndophily<-model2$coefficients$fixed
Endophily<- exp(as.vector(logEndophily))
speciesVector <- substr(names(logEndophily), 8, nchar(names(logEndophily)))
j <- 1
i <- 1
while (i < length(speciesVector)){
  if (taxon[[j]][[1]][[1]]==speciesVector[i]) {
    taxon[[j]]$Endophily <- Endophily[i]
    i <- i+1
  }
  j <- j+1
}

# Check relationship between indoor resting from experimental huts and ratios indoor resting: biting
with(taxon, plot(log(proportion_indoor_resting),log(indoor_resting_by_biting)))
# use the endophily data as estimates of indoor resting where experimental hut data are not available
# assuming the above relationship of endophily and indoor resting is supported by the plot
taxon$proportion_indoor_resting <- taxon$Endophily/(1 + taxon[[i]]$Endophily)



# POPULATE taxon WITH ENDOPHAGY DATA
logitendo <- log(EndophagyData$value/(1.01-EndophagyData$value))
model1 <- lme(logitendo~0 + species,random=~1|Site,method="REML")
summary(model1)
logitEndophagy<-model1$coefficients$fixed
taxon$Endophagy<- 1/(1+exp(-logitEndophagy))





