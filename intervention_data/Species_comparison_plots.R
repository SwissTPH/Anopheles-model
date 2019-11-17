#=============================================================================

rm(list = ls())

# ---------------------------------------------------------
user <- 'Tom'
# Tom

#user <- 'Olivier' 
# Olivier




### Required libraries and functions
library(plyr)
library(dplyr)
library(ggplot2)
library(scales)


library(lubridate)

# ---------------------------------------------------------
# A: Set the working directory and paths for plots and scripts


if (user=='Tom'){
  setwd("C:/git_repos/om-workflow/vector_specific_parameterisation/albimanus/species_comparison")
  pathToSavePlots = 'C:/git_repos/om-workflow/vector_specific_parameterisation/albimanus/species_comparison/'
  scriptPath_omworkflow = 'C:/git_repos/om-workflow/'
  scriptPath  = 'C:/git_repos/om-workflow/vector_specific_parameterisation/species_comparisons/'
}

if (user=='Olivier'){
	setwd('C:/OB/git_repos/om-workflow/vector_specific_parameterisation/albimanus/species_comparison')
	pathToSavePlots = 'C:/OB/git_repos/om-workflow/vector_specific_parameterisation/albimanus/species_comparison/'
	scriptPath_omworkflow = 'C:/OB/git_repos/om-workflow/'
	scriptPath  = 'C:/OB/git_repos/om-workflow/vector_specific_parameterisation/species_comparisons/'	
}


speciesList<- c('albimanus','gambiaess')

# set colour levels
# define colourBlind safe colors
recomCol = c("#D55E00", "#009E73", "#0072A7","#CC79A7") # slightly more blue pink/purple #C879C8 to replace pink#CC79A7
my_recomCol = c("#D55E00", "#009E73", "#0072A7","#C879C8")
my_col_Int=c(my_recomCol)#,"#009E73")

labeli2 <- function(variable, value){
  value <- droplevels(value)
  names_li <- list("albimanus"="An. albimanus", "gambiaess"="An. gambiae s.l.")
  return(names_li[value])
}

hour <- as.POSIXct(strptime(c("2016-12-17 18:30","2016-12-17 19:30","2016-12-17 20:30","2016-12-17 21:30","2016-12-17 22:30",
                              "2016-12-17 23:30","2016-12-18 00:30","2016-12-18 01:30","2016-12-18 02:30",
                              "2016-12-18 03:30","2016-12-18 04:30","2016-12-18 05:30","2016-12-18 06:30")
                            ,format = "%Y-%m-%d %H:%M"))
# PLOTS OF ACTIVITY CYCLES 
Activity_cycle<-data.frame() 
for (species in speciesList){
  temp <- read.csv(file=paste0(scriptPath_omworkflow,"vector_specific_parameterisation/",species,"/ActivityCycle_time.csv"), skip = 1, header = FALSE, sep=" ")            
  temp$Percent_in_bed <- temp$V4
  temp$hour <- hour
  temp$Country<-rep(species,length(temp$hour))
  keeps <- c("Country","hour", "Percent_in_bed")
  Activity_cycle <- rbind(Activity_cycle,temp[keeps])
}                          

activity <- ggplot(data = Activity_cycle, 
              aes(y = Percent_in_bed, x = hour, group = Country, colour = Country))+  
  theme_bw(base_size = 24)+
  geom_smooth(se = FALSE, size=2) +
  theme(legend.position=c(0.8,0.2))+
  ylab("Percentage of population in sleeping spaces")+
  xlab("Time")+
  scale_color_manual(name="Country", values=my_col_Int, labels = c("Haiti", "Tanzania"))+
  scale_x_datetime(breaks = date_breaks("2 hour"), labels = date_format("%H:%M"))
  theme(axis.text.x = element_text(size = rel(2.0),color="black"))+
  theme(axis.text.y = element_text(size = rel(2.0),color="black"))+
  theme(axis.title.x = element_text(size = rel(2.0),color="black"))+
  theme(axis.title.y = element_text(size = rel(2.0),color="black"))

ggsave(activity, file=paste(pathToSavePlots, "Activity_cycle.png", sep=""), width = 200, height = 200, units = "mm")


png(filename = "Activity_cycleOB.png", width = 1.5*620, height = 1.5*620, pointsize = 12, bg = "white", res = 144)
op <- par(mfrow = c(1, 1), mar= c(4, 4, 1, 1) + 0.1, oma=c(0,0,0,0))
plot(c(0), type="n", ylim=c(0,100), xlim=as.numeric(c(as.POSIXlt('2016-12-17 18:30:00'),as.POSIXlt('2016-12-18 06:30:00'))), xlab="Time", ylab="Population in sleeping spaces (%)", yaxt="n", xaxt="n")
axis(1, at=as.numeric(c(as.POSIXlt('2016-12-17 18:30:00'), as.POSIXlt('2016-12-17 22:30:00'), as.POSIXlt('2016-12-18 02:30:00'), as.POSIXlt('2016-12-18 6:30:00'))), labels=c('18:30', '22:30', '02:30', '06:30'))
axis(2)
Datasub <- subset(Activity_cycle, Country=='albimanus')
lines(as.numeric(as.POSIXlt(Datasub$hour)), Datasub$Percent_in_bed, col=my_col_Int[1], lty=1, lwd=2)
Datasub <- subset(Activity_cycle, Country=='gambiaess')
lines(as.numeric(as.POSIXlt(Datasub$hour)), Datasub$Percent_in_bed, col=my_col_Int[2], lty=1, lwd=2)
legend(x=as.numeric(c(as.POSIXlt('2016-12-17 18:30:00'))), y=100, c("Haiti", "Tanzania"), col=c(my_col_Int[1:2]), lty=c(1), lwd=2, bg="white")
par(op)
dev.off()




# PLOT OF BITING CYCLES 
Biting_cycle<-data.frame() 
for (species in speciesList){
  temp <- read.csv(file=paste0(scriptPath_omworkflow,"vector_specific_parameterisation/",species,"/BitingCycle_time.csv"), skip = 1, header = FALSE, sep=" ")            
  totals <- tapply(temp[complete.cases(temp),5], temp[complete.cases(temp),6], sum)
  totals <- c(rep(totals[1],26))
  temp$Human_landing <- temp$V5/totals
  temp$Location <- temp$V6
  temp$hour <- c(hour,hour)
  temp$species<-temp$V2
  keeps <- c("species","hour", "Human_landing", "Location")
  Biting_cycle <- rbind(Biting_cycle,temp[keeps])
} 

biting <- ggplot(data = Biting_cycle, 
              aes(y = Human_landing, x = hour, group = Location, colour = Location))+  
  theme_bw(base_size = 18)+
  theme(legend.position=c(0.3,0.8))+
  scale_y_continuous(limits=c(0,max(Biting_cycle$Human_landing, na.rm = FALSE)))+
  geom_point(size=2)+
  geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3), se = FALSE, size=2)+
  ylab("Relative biting rate")+
  xlab("Time")+
  scale_color_manual(name="Location", values=my_col_Int, labels = c("Indoor", "Outdoor"))+
  scale_x_datetime(breaks = date_breaks("4 hour"), labels = date_format("%H:%M"))+
  theme(axis.text.x = element_text(size = rel(1.5),color="black"))+
  theme(axis.text.y = element_text(size = rel(1.5),color="black"))+
  theme(axis.title.x = element_text(size = rel(1.5),color="black"))+
  theme(axis.title.y = element_text(size = rel(1.5),color="black"))+
  facet_grid(.~ species, labeller=labeli2) 

ggsave(biting, file=paste(pathToSavePlots, "Biting_cycle.png", sep=""), width = 350, height = 200, units = "mm")




# PLOT OF BITING CYCLES OB 
Biting_cycleOB<-data.frame() 
for (species in speciesList){
  temp <- read.csv(file=paste0(scriptPath_omworkflow,"vector_specific_parameterisation/",species,"/BitingCycle_time.csv"), skip = 1, header = FALSE, sep=" ")            
  totals <- tapply(temp[complete.cases(temp),5], temp[complete.cases(temp),6], sum)
  totals <- c(rep(sum(totals),26)) #now the sum of both indoors and outdoors
  temp$Human_landing <- temp$V5/totals
  temp$Location <- temp$V6
  temp$hour <- c(hour,hour)
  temp$species<-temp$V2
  keeps <- c("species","hour", "Human_landing", "Location")
  Biting_cycleOB <- rbind(Biting_cycleOB,temp[keeps])
} 


png(filename = "Biting_cycleOB.png", width = 3*620, height = 1.5*620, pointsize = 12, bg = "white", res = 144)
op <- par(mfrow = c(1, 2), mar= c(0, 1, 1, 1) + 0.1, oma=c(4,4,2,0))
plot(c(0), type="n", ylim=c(0,0.5), xlim=as.numeric(c(as.POSIXlt('2016-12-17 18:30:00'),as.POSIXlt('2016-12-18 06:30:00'))), ylab="", xlab="", yaxt="n", xaxt="n")
axis(1, at=as.numeric(c(as.POSIXlt('2016-12-17 18:30:00'), as.POSIXlt('2016-12-17 22:30:00'), as.POSIXlt('2016-12-18 02:30:00'), as.POSIXlt('2016-12-18 6:30:00'))), labels=c('18:30', '22:30', '02:30', '06:30'))
axis(2)
Datasub <- subset(Biting_cycleOB, species=='albimanus' & Location=='Indoor')
lines(as.numeric(as.POSIXlt(Datasub$hour)), Datasub$Human_landing, col=my_col_Int[1], lty=1, lwd=2)
Datasub <- subset(Biting_cycleOB, species=='albimanus' & Location=='Outdoor')
lines(as.numeric(as.POSIXlt(Datasub$hour)), Datasub$Human_landing, col=my_col_Int[2], lty=1, lwd=2)
legend(x=as.numeric(c(as.POSIXlt('2016-12-17 22:30:00'))), y=1.8, c("indoor", "outdoor"), col=my_col_Int[1:2], lty=1, lwd=2, bg="white")
mtext(side=3, expression(paste(italic("An. albimanus"), sep="")), line=1, outer=FALSE)
mtext(side=2, "Proportion of overall bites", line=3)
plot(c(0), type="n", ylim=c(0,0.5), xlim=as.numeric(c(as.POSIXlt('2016-12-17 18:30:00'),as.POSIXlt('2016-12-18 06:30:00'))), ylab="", xlab="", yaxt="n", xaxt="n")
axis(1, at=as.numeric(c(as.POSIXlt('2016-12-17 18:30:00'), as.POSIXlt('2016-12-17 22:30:00'), as.POSIXlt('2016-12-18 02:30:00'), as.POSIXlt('2016-12-18 6:30:00'))), labels=c('18:30', '22:30', '02:30', '06:30'))
axis(2, labels=FALSE)
Datasub <- subset(Biting_cycleOB, species=='Gambiaess' & Location=='Indoor')
lines(as.numeric(as.POSIXlt(Datasub$hour)), Datasub$Human_landing, col=my_col_Int[1], lty=1, lwd=2)
Datasub <- subset(Biting_cycleOB, species=='Gambiaess' & Location=='Outdoor')
lines(as.numeric(as.POSIXlt(Datasub$hour)), Datasub$Human_landing, col=my_col_Int[2], lty=1, lwd=2)
mtext(side=3, expression(paste(italic("An. gambiae"), " s.s.", sep="")), line=1, outer=FALSE)
mtext(side=1, "Time", line=3, outer=TRUE)
par(op)
dev.off()
	



# PLOT OF EFFECT OF NETS BY TIME  
NetEffect<-data.frame() 
for (species in speciesList){
  temp <- read.csv(file=paste0(scriptPath_omworkflow,"vector_specific_parameterisation/",species,"/NetEffect_time.csv"), skip = 1, header = FALSE, sep=" ")            
  totals <- mean(temp[complete.cases(temp),5])
  temp$Exposure <- temp$V5/totals
  temp$Location <- temp$V6
  temp$hour <- c(hour,hour)
  temp$species<-temp$V2
  keeps <- c("species","hour", "Exposure", "Location")
  NetEffect <- rbind(NetEffect,temp[keeps])
} 

NetEffect_time <- ggplot(NetEffect, aes(y = Exposure, x = hour, color = Location)) +
#  geom_smooth(method = "lm", formula = y ~ splines::bs(x, 4), se = FALSE, size=2) +
  geom_smooth(se = FALSE, size=2) +
  xlab("Time")+
  ylab("Relative Exposure (bites per hour)")+
  theme_bw(base_size = 18)+
  theme(legend.position=c(0.3,0.8))+
  scale_color_manual(name="User type", values=my_col_Int, labels = c("net user", "non-user"))+
  scale_x_datetime(breaks = date_breaks("4 hour"), labels = date_format("%H:%M"))+
  theme(axis.text.x = element_text(size = rel(1.5),color="black"))+
  theme(axis.text.y = element_text(size = rel(1.5),color="black"))+
  theme(axis.title.x = element_text(size = rel(1.5),color="black"))+
  theme(axis.title.y = element_text(size = rel(1.5),color="black"))+
  facet_grid(.~ species, labeller=labeli2) 
ggsave(NetEffect_time, file=paste(pathToSavePlots, "NetEffect_time.png", sep=""), width = 350, height = 200, units = "mm")




# PLOT OF EFFECT OF NETS BY TIME OB 
NetEffectOB<-data.frame() 
for (species in speciesList){
  temp <- read.csv(file=paste0(scriptPath_omworkflow,"vector_specific_parameterisation/",species,"/NetEffect_time.csv"), skip = 1, header = FALSE, sep=" ")            
  #totals <- mean(temp[complete.cases(temp),5])
  totals <- tapply(temp[complete.cases(temp),5], temp[complete.cases(temp),6], sum)
  totals <- c(rep(totals[2],26)) #now the sum of Non-user (only)
  temp$Exposure <- temp$V5/totals
  temp$Location <- temp$V6
  temp$hour <- c(hour,hour)
  temp$species<-temp$V2
  keeps <- c("species","hour", "Exposure", "Location")
  NetEffectOB <- rbind(NetEffectOB,temp[keeps])
} 


png(filename = "NetEffect_timeOB.png", width = 3*620, height = 1.5*620, pointsize = 12, bg = "white", res = 144)
op <- par(mfrow = c(1, 2), mar= c(0, 1, 1, 1) + 0.1, oma=c(4,4,2,0))
plot(c(0), type="n", ylim=c(0,0.7), xlim=as.numeric(c(as.POSIXlt('2016-12-17 18:30:00'),as.POSIXlt('2016-12-18 06:30:00'))), ylab="", xlab="", yaxt="n", xaxt="n")
axis(1, at=as.numeric(c(as.POSIXlt('2016-12-17 18:30:00'), as.POSIXlt('2016-12-17 22:30:00'), as.POSIXlt('2016-12-18 02:30:00'), as.POSIXlt('2016-12-18 6:30:00'))), labels=c('18:30', '22:30', '02:30', '06:30'))
axis(2)
Datasub <- subset(NetEffectOB, species=='albimanus' & Location=='Net user')
lines(as.numeric(as.POSIXlt(Datasub$hour)), Datasub$Exposure, col=my_col_Int[1], lty=1, lwd=2)
Datasub <- subset(NetEffectOB, species=='albimanus' & Location=='Non-user')
lines(as.numeric(as.POSIXlt(Datasub$hour)), Datasub$Exposure, col=my_col_Int[2], lty=1, lwd=2)
legend(x=as.numeric(c(as.POSIXlt('2016-12-17 22:30:00'))), y=0.7, c("net user", "non-user"), col=my_col_Int[1:2], lty=1, lwd=2, bg="white")
mtext(side=3, expression(paste(italic("An. albimanus"), sep="")), line=1, outer=FALSE)
mtext(side=2, "Proportion of overall bites of a non-user", line=3)
plot(c(0), type="n", ylim=c(0,0.7), xlim=as.numeric(c(as.POSIXlt('2016-12-17 18:30:00'),as.POSIXlt('2016-12-18 06:30:00'))), ylab="", xlab="", yaxt="n", xaxt="n")
axis(1, at=as.numeric(c(as.POSIXlt('2016-12-17 18:30:00'), as.POSIXlt('2016-12-17 22:30:00'), as.POSIXlt('2016-12-18 02:30:00'), as.POSIXlt('2016-12-18 6:30:00'))), labels=c('18:30', '22:30', '02:30', '06:30'))
axis(2, labels=FALSE)
Datasub <- subset(NetEffectOB, species=='Gambiaess' & Location=='Net user')
lines(as.numeric(as.POSIXlt(Datasub$hour)), Datasub$Exposure, col=my_col_Int[1], lty=1, lwd=2)
Datasub <- subset(NetEffectOB, species=='Gambiaess' & Location=='Non-user')
lines(as.numeric(as.POSIXlt(Datasub$hour)), Datasub$Exposure, col=my_col_Int[2], lty=1, lwd=2)
mtext(side=3, expression(paste(italic("An. gambiae"), " s.s.", sep="")), line=1, outer=FALSE)
mtext(side=1, "Time", line=3, outer=TRUE)
par(op)
dev.off()
	

	

png(filename = "AlignmentCombinedOB.png", width = 2*620, height = 2.5*620, pointsize = 12, bg = "white", res = 144)
op <- par(mfrow = c(3, 2), mar= c(0, 1, 1, 1) + 0.1, oma=c(4,4,2,0))

plot(c(0), type="n", ylim=c(0,0.5), xlim=as.numeric(c(as.POSIXlt('2016-12-17 18:30:00'),as.POSIXlt('2016-12-18 06:30:00'))), ylab="", xlab="", yaxt="n", xaxt="n")
axis(1, at=as.numeric(c(as.POSIXlt('2016-12-17 18:30:00'), as.POSIXlt('2016-12-17 22:30:00'), as.POSIXlt('2016-12-18 02:30:00'), as.POSIXlt('2016-12-18 6:30:00'))), labels=FALSE)
axis(2)
Datasub <- subset(Biting_cycleOB, species=='albimanus' & Location=='Indoor')
lines(as.numeric(as.POSIXlt(Datasub$hour)), Datasub$Human_landing, col=my_col_Int[1], lty=1, lwd=2)
Datasub <- subset(Biting_cycleOB, species=='albimanus' & Location=='Outdoor')
lines(as.numeric(as.POSIXlt(Datasub$hour)), Datasub$Human_landing, col=my_col_Int[2], lty=1, lwd=2)
legend(x=as.numeric(c(as.POSIXlt('2016-12-17 22:30:00'))), y=0.5, c("indoor", "outdoor"), col=my_col_Int[1:2], lty=1, lwd=2, bg="white")
mtext(side=3, expression(paste(italic("An. albimanus")," in Haiti", sep="")), line=1, outer=FALSE)
mtext(side=2, "Proportion of overall bites", line=3)
text(x=as.numeric(c(as.POSIXlt('2016-12-18 6:30:00'))), y=0.5, "a", cex=2)
plot(c(0), type="n", ylim=c(0,0.5), xlim=as.numeric(c(as.POSIXlt('2016-12-17 18:30:00'),as.POSIXlt('2016-12-18 06:30:00'))), ylab="", xlab="", yaxt="n", xaxt="n")
axis(1, at=as.numeric(c(as.POSIXlt('2016-12-17 18:30:00'), as.POSIXlt('2016-12-17 22:30:00'), as.POSIXlt('2016-12-18 02:30:00'), as.POSIXlt('2016-12-18 6:30:00'))), labels=FALSE)
axis(2, labels=FALSE)
Datasub <- subset(Biting_cycleOB, species=='Gambiaess' & Location=='Indoor')
lines(as.numeric(as.POSIXlt(Datasub$hour)), Datasub$Human_landing, col=my_col_Int[1], lty=1, lwd=2)
Datasub <- subset(Biting_cycleOB, species=='Gambiaess' & Location=='Outdoor')
lines(as.numeric(as.POSIXlt(Datasub$hour)), Datasub$Human_landing, col=my_col_Int[2], lty=1, lwd=2)
mtext(side=3, expression(paste(italic("An. gambiae"), " s.s. in Tanzania", sep="")), line=1, outer=FALSE)
mtext(side=1, "Time", line=3, outer=TRUE)
text(x=as.numeric(c(as.POSIXlt('2016-12-18 6:30:00'))), y=0.5, "b", cex=2)

plot(c(0), type="n", ylim=c(0,100), xlim=as.numeric(c(as.POSIXlt('2016-12-17 18:30:00'),as.POSIXlt('2016-12-18 06:30:00'))), xlab="Time", ylab="", yaxt="n", xaxt="n")
mtext(side=2, "Population in sleeping spaces (%)", line=3)
axis(1, at=as.numeric(c(as.POSIXlt('2016-12-17 18:30:00'), as.POSIXlt('2016-12-17 22:30:00'), as.POSIXlt('2016-12-18 02:30:00'), as.POSIXlt('2016-12-18 6:30:00'))), labels=FALSE)
axis(2)
Datasub <- subset(Activity_cycle, Country=='albimanus')
lines(as.numeric(as.POSIXlt(Datasub$hour)), Datasub$Percent_in_bed, col=my_col_Int[1], lty=1, lwd=2)
#Datasub <- subset(Activity_cycle, Country=='gambiaess')
#lines(as.numeric(as.POSIXlt(Datasub$hour)), Datasub$Percent_in_bed, col=my_col_Int[2], lty=1, lwd=2)
#legend(x=as.numeric(c(as.POSIXlt('2016-12-17 18:30:00'))), y=100, c("Haiti", "Tanzania"), col=c(my_col_Int[1:2]), lty=c(1), lwd=2, bg="white")
text(x=as.numeric(c(as.POSIXlt('2016-12-18 6:30:00'))), y=100, "c", cex=2)
plot(c(0), type="n", ylim=c(0,100), xlim=as.numeric(c(as.POSIXlt('2016-12-17 18:30:00'),as.POSIXlt('2016-12-18 06:30:00'))), xlab="Time", ylab="", yaxt="n", xaxt="n")
axis(1, at=as.numeric(c(as.POSIXlt('2016-12-17 18:30:00'), as.POSIXlt('2016-12-17 22:30:00'), as.POSIXlt('2016-12-18 02:30:00'), as.POSIXlt('2016-12-18 6:30:00'))), labels=FALSE)
axis(2, labels=FALSE)
#Datasub <- subset(Activity_cycle, Country=='albimanus')
#lines(as.numeric(as.POSIXlt(Datasub$hour)), Datasub$Percent_in_bed, col=my_col_Int[1], lty=1, lwd=2)
Datasub <- subset(Activity_cycle, Country=='gambiaess')
lines(as.numeric(as.POSIXlt(Datasub$hour)), Datasub$Percent_in_bed, col=my_col_Int[1], lty=1, lwd=2)
#legend(x=as.numeric(c(as.POSIXlt('2016-12-17 18:30:00'))), y=100, c("Haiti", "Tanzania"), col=c(my_col_Int[1:2]), lty=c(1), lwd=2, bg="white")
text(x=as.numeric(c(as.POSIXlt('2016-12-18 6:30:00'))), y=100, "d", cex=2)

plot(c(0), type="n", ylim=c(0,0.7), xlim=as.numeric(c(as.POSIXlt('2016-12-17 18:30:00'),as.POSIXlt('2016-12-18 06:30:00'))), ylab="", xlab="", yaxt="n", xaxt="n")
axis(1, at=as.numeric(c(as.POSIXlt('2016-12-17 18:30:00'), as.POSIXlt('2016-12-17 22:30:00'), as.POSIXlt('2016-12-18 02:30:00'), as.POSIXlt('2016-12-18 6:30:00'))), labels=c('18:30', '22:30', '02:30', '06:30'))
axis(2)
Datasub <- subset(NetEffectOB, species=='albimanus' & Location=='Net user')
lines(as.numeric(as.POSIXlt(Datasub$hour)), Datasub$Exposure, col=my_col_Int[1], lty=1, lwd=2)
Datasub <- subset(NetEffectOB, species=='albimanus' & Location=='Non-user')
lines(as.numeric(as.POSIXlt(Datasub$hour)), Datasub$Exposure, col=my_col_Int[2], lty=1, lwd=2)
legend(x=as.numeric(c(as.POSIXlt('2016-12-17 22:30:00'))), y=0.7, c("net user", "non-user"), col=my_col_Int[1:2], lty=1, lwd=2, bg="white")
#mtext(side=3, expression(paste(italic("An. albimanus"), sep="")), line=1, outer=FALSE)
mtext(side=2, "Proportion of overall bites of a non-user", line=3)
text(x=as.numeric(c(as.POSIXlt('2016-12-18 6:30:00'))), y=0.7, "e", cex=2)
plot(c(0), type="n", ylim=c(0,0.7), xlim=as.numeric(c(as.POSIXlt('2016-12-17 18:30:00'),as.POSIXlt('2016-12-18 06:30:00'))), ylab="", xlab="", yaxt="n", xaxt="n")
axis(1, at=as.numeric(c(as.POSIXlt('2016-12-17 18:30:00'), as.POSIXlt('2016-12-17 22:30:00'), as.POSIXlt('2016-12-18 02:30:00'), as.POSIXlt('2016-12-18 6:30:00'))), labels=c('18:30', '22:30', '02:30', '06:30'))
axis(2, labels=FALSE)
Datasub <- subset(NetEffectOB, species=='Gambiaess' & Location=='Net user')
lines(as.numeric(as.POSIXlt(Datasub$hour)), Datasub$Exposure, col=my_col_Int[1], lty=1, lwd=2)
Datasub <- subset(NetEffectOB, species=='Gambiaess' & Location=='Non-user')
lines(as.numeric(as.POSIXlt(Datasub$hour)), Datasub$Exposure, col=my_col_Int[2], lty=1, lwd=2)
#mtext(side=3, expression(paste(italic("An. gambiae"), " s.s.", sep="")), line=1, outer=FALSE)
mtext(side=1, "Time", line=3, outer=TRUE)
text(x=as.numeric(c(as.POSIXlt('2016-12-18 6:30:00'))), y=0.7, "f", cex=2)
par(op)
dev.off()
	
	
	
	


vectorInterventions <- c('No Intervention', 'LLINs','Permethrin IRS', 'DDT IRS', 'Bendiocarb IRS', 'House Screening', 'LLINs.albimanus')
speciesList<- c('albimanus','gambiaess')
interventions <- c(vectorInterventions,'Larviciding','Case Management','RTS,S mass vaccination','Test and Treat')
EntoData<-data.frame() 
for (species in speciesList){
  temp <- read.csv(file=paste0(scriptPath_omworkflow,"vector_specific_parameterisation/",species,"/reduction_Rc_coverage.csv"), header=TRUE, sep=" ")            
  EntoData <- rbind(EntoData,temp)
}
EntoData$Coverage <- EntoData$Coverage*100

interventions <- c('Actellic 300 CS','Icon 10 CS','Actellic 50 EC','Permethrin IRS', 'DDT IRS', 'Bendiocarb IRS')

EntoData$Intervention <- revalue(EntoData$Intervention, c('Bendiocarb IRS'='Bendiocarb', 'DDT IRS'='DDT'))


# plot comparison of IRS
interventions <- c('Actellic 300 CS','Icon 10 CS','DDT', 'Bendiocarb')
irsplot <- ggplot(data = EntoData %>% filter(Intervention %in% interventions), aes(y = Reduction_Rc, x = Coverage, colour =Intervention))+
  geom_line(size = 2)+ scale_color_manual(name="Insecticide", values=my_col_Int)+
  theme_bw(base_size = 18)+
  theme(legend.position=c(0.3,0.8))+
  theme(axis.text.x = element_text(size = rel(1.5),color="black"))+
  theme(axis.text.y = element_text(size = rel(1.5),color="black"))+
  theme(axis.title.x = element_text(size = rel(1.5),color="black"))+
  theme(axis.title.y = element_text(size = rel(1.5),color="black"))+
  facet_grid(.~ required_species, labeller=labeli2)  +
  xlab("Coverage (%)") + ylab("Reduction in vectorial capacity (%)") 

ggsave(irsplot, file=paste(pathToSavePlots, "ComparisonEffectsIRS.png", sep=""), width = 350, height = 200, units = "mm")



png(filename = "ComparisonEffectsIRSOB.png", width = 3*620, height = 1.5*620, pointsize = 12, bg = "white", res = 144)
op <- par(mfrow = c(1, 2), mar= c(0, 1, 1, 1) + 0.1, oma=c(4,4,2,0))
plot(c(0), type="n", ylim=c(0,100), xlim=c(0,100), ylab="", xlab="", yaxt="n", xaxt="n")
axis(1)
axis(2)
EntoDatasub <- subset(EntoData, required_species=='albimanus' & Intervention=='Bendiocarb IRS albimanus but outside until bedtime')
lines(EntoDatasub$Coverage, EntoDatasub$Reduction_Rc, col=my_col_Int[1], lty=1, lwd=2)
EntoDatasub <- subset(EntoData, required_species=='albimanus' & Intervention=='Bendiocarb IRS albimanus')
lines(EntoDatasub$Coverage, EntoDatasub$Reduction_Rc, col=my_col_Int[1], lty=2, lwd=2)
EntoDatasub <- subset(EntoData, required_species=='Gambiaess' & Intervention=='Bendiocarb IRS gambiae but outside until bedtime')
lines(EntoDatasub$Coverage, EntoDatasub$Reduction_Rc, col=my_col_Int[2], lty=1, lwd=2)
EntoDatasub <- subset(EntoData, required_species=='Gambiaess' & Intervention=='Bendiocarb IRS gambiae')
lines(EntoDatasub$Coverage, EntoDatasub$Reduction_Rc, col=my_col_Int[2], lty=2, lwd=2)
legend(x=0, y=100, c(expression(paste(italic("An. albimanus"), ", people outside until bedtime", sep="")), expression(paste(italic("An. albimanus"), ", people inside" , sep="")),expression(paste(italic("An. gambiae"), " s.s., people outside until bedtime", sep="")),expression(paste(italic("An. gambiae"), " s.s., people inside", sep=""))), col=rep(c(my_col_Int[1], my_col_Int[2]), each=2), lty=c(rep(c(1,2),2)), lwd=2, bg="white")
mtext(side=3, "Bendiocarb", line=1, outer=FALSE)
mtext(side=2, "Reduction in vectorial capacity (%)", line=3)
plot(c(0), type="n", ylim=c(0,100), xlim=c(0,100), ylab="", xlab="", yaxt="n", xaxt="n")
axis(1)
axis(2, labels=FALSE)
EntoDatasub <- subset(EntoData, required_species=='albimanus' & Intervention=='DDT IRS albimanus but outside until bedtime')
lines(EntoDatasub$Coverage, EntoDatasub$Reduction_Rc, col=my_col_Int[1], lty=1, lwd=2)
EntoDatasub <- subset(EntoData, required_species=='albimanus' & Intervention=='DDT IRS albimanus')
lines(EntoDatasub$Coverage, EntoDatasub$Reduction_Rc, col=my_col_Int[1], lty=2, lwd=2)
EntoDatasub <- subset(EntoData, required_species=='Gambiaess' & Intervention=='DDT IRS gambiae but outside until bedtime')
lines(EntoDatasub$Coverage, EntoDatasub$Reduction_Rc, col=my_col_Int[2], lty=1, lwd=2)
EntoDatasub <- subset(EntoData, required_species=='Gambiaess' & Intervention=='DDT IRS gambiae')
lines(EntoDatasub$Coverage, EntoDatasub$Reduction_Rc, col=my_col_Int[2], lty=2, lwd=2)
mtext(side=3, "DDT", line=1, outer=FALSE)
mtext(side=1, "Coverage (%)", line=3, outer=TRUE)
par(op)
dev.off()
	




# plot comparison of screening and LLINs
interventions <- c('LLINs','House Screening')

screeningfig <- ggplot(data = EntoData %>% filter(Intervention %in% interventions), aes(y = Reduction_Rc, x = Coverage, colour =required_species))+
  geom_line(size = 2)+ scale_color_manual(name="Species", values=my_col_Int, labels = c("An. albimanus", "An. gambiae s.l."))+
  theme_bw(base_size = 18) +
  theme(legend.position=c(0.2,0.8)) +
  theme(axis.text.x = element_text(size = rel(1.5),color="black"))+
  theme(axis.text.y = element_text(size = rel(1.5),color="black"))+
  theme(axis.title.x = element_text(size = rel(1.5),color="black"))+
  theme(axis.title.y = element_text(size = rel(1.5),color="black"))+
  scale_x_continuous(limits = c(0, 100))+
  scale_y_continuous(limits = c(0, 100))+	
  #geom_abline(intercept = 0, slope = 1) +
  facet_grid(.~ Intervention)  + 
  xlab("Coverage (%)") + ylab("Reduction in vectorial capacity (%)") 

ggsave(screeningfig, file=paste(pathToSavePlots, "ComparisonEffectsLLINScreening.png", sep=""), width = 350, height = 200, units = "mm")


png(filename = "ComparisonEffectsLLINScreeningOB.png", width = 3*620, height = 1.5*620, pointsize = 12, bg = "white", res = 144)
op <- par(mfrow = c(1, 2), mar= c(0, 1, 1, 1) + 0.1, oma=c(4,4,2,0))
plot(c(0), type="n", ylim=c(0,100), xlim=c(0,100), ylab="", xlab="", yaxt="n", xaxt="n")
axis(1)
axis(2)
EntoDatasub <- subset(EntoData, required_species=='albimanus' & Intervention=='House Screening but outside until bedtime')
lines(EntoDatasub$Coverage, EntoDatasub$Reduction_Rc, col=my_col_Int[1], lty=1, lwd=2)
EntoDatasub <- subset(EntoData, required_species=='albimanus' & Intervention=='House Screening')
lines(EntoDatasub$Coverage, EntoDatasub$Reduction_Rc, col=my_col_Int[1], lty=2, lwd=2)
EntoDatasub <- subset(EntoData, required_species=='Gambiaess' & Intervention=='House Screening but outside until bedtime')
lines(EntoDatasub$Coverage, EntoDatasub$Reduction_Rc, col=my_col_Int[2], lty=1, lwd=2)
EntoDatasub <- subset(EntoData, required_species=='Gambiaess' & Intervention=='House Screening')
lines(EntoDatasub$Coverage, EntoDatasub$Reduction_Rc, col=my_col_Int[2], lty=2, lwd=2)
legend(x=0, y=100, c(expression(paste(italic("An. albimanus"), ", people outside until bedtime", sep="")), expression(paste(italic("An. albimanus"), ", people inside" , sep="")),expression(paste(italic("An. gambiae"), " s.s., people outside until bedtime", sep="")),expression(paste(italic("An. gambiae"), " s.s., people inside", sep=""))), col=rep(c(my_col_Int[1], my_col_Int[2]), each=2), lty=c(rep(c(1,2),2)), lwd=2, bg="white")
mtext(side=3, "House screening", line=1, outer=FALSE)
mtext(side=2, "Reduction in vectorial capacity (%)", line=3)
plot(c(0), type="n", ylim=c(0,100), xlim=c(0,100), ylab="", xlab="", yaxt="n", xaxt="n")
axis(1)
axis(2, labels=FALSE)
EntoDatasub <- subset(EntoData, required_species=='albimanus' & Intervention=='LLINs lambdacyhalothrin albimanus')
lines(EntoDatasub$Coverage, EntoDatasub$Reduction_Rc, col=my_col_Int[1], lty=1, lwd=2)
EntoDatasub <- subset(EntoData, required_species=='albimanus' & Intervention=='LLINs lambdacyhalothrin albimanus inside house before biting starts')
lines(EntoDatasub$Coverage, EntoDatasub$Reduction_Rc, col=my_col_Int[1], lty=2, lwd=2)
EntoDatasub <- subset(EntoData, required_species=='Gambiaess' & Intervention=='LLINs lambdacyhalothrin gambiae')
lines(EntoDatasub$Coverage, EntoDatasub$Reduction_Rc, col=my_col_Int[2], lty=1, lwd=2)
EntoDatasub <- subset(EntoData, required_species=='Gambiaess' & Intervention=='LLINs lambdacyhalothrin gambiae inside house before biting starts')
lines(EntoDatasub$Coverage, EntoDatasub$Reduction_Rc, col=my_col_Int[2], lty=2, lwd=2)
mtext(side=3, "ITNs", line=1, outer=FALSE)
mtext(side=1, "Coverage (%)", line=3, outer=TRUE)
par(op)
dev.off()

# plot human-side interventions
interventions <- c('Case Management','RTS,S mass vaccination','Test and Treat')
human <- ggplot(data = EntoData %>% filter(required_species == 'albimanus') %>% filter(Intervention %in% interventions), 
  aes(y = Reduction_Rc, x = Coverage, colour = Intervention))+
  theme_bw(base_size = 24)+
  stat_smooth(size = 2, method = "loess",se=FALSE) +
  scale_color_manual(name="Intervention", values=my_col_Int)+
  geom_abline(intercept = 0, slope = 1) + 
  theme(legend.position=c(0.3,0.8))+
  xlab("Coverage (%)") + ylab("Reduction in Reproduction Number (%)") 

ggsave(human, file=paste(pathToSavePlots, "HumanInterventions.png", sep=""), width = 200, height = 200, units = "mm")


png(filename = "HumanInterventionsOB.png", width = 1.5*620, height = 1.5*620, pointsize = 12, bg = "white", res = 144)
op <- par(mfrow = c(1, 1), mar= c(4, 4, 1, 1) + 0.1, oma=c(0,0,0,0))
plot(c(0), type="n", ylim=c(0,100), xlim=c(0,100), xlab="Coverage (%)", ylab="Reduction in reproduction number (%)")
EntoDatasub <- subset(EntoData, required_species=='albimanus' & Intervention=='Case Management')
lines(EntoDatasub$Coverage, EntoDatasub$Reduction_Rc, col=my_col_Int[1], lty=1, lwd=2)
EntoDatasub <- subset(EntoData, required_species=='albimanus' & Intervention=='RTS,S mass vaccination')
lines(EntoDatasub$Coverage, EntoDatasub$Reduction_Rc, col=my_col_Int[2], lty=1, lwd=2)
EntoDatasub <- subset(EntoData, required_species=='albimanus' & Intervention=='Test and Treat')
lines(EntoDatasub$Coverage, EntoDatasub$Reduction_Rc, col=my_col_Int[3], lty=1, lwd=2)
legend(x=0, y=100, c("Case management", "RTS,S mass vaccination", "Test and treat"), col=c(my_col_Int[1:3]), lty=c(1), lwd=2, bg="white")
par(op)
dev.off()
