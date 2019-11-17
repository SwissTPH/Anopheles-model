ds <- "C:/git_repos/Anopheles-model/raw_data/"

# Read the parameters from the table in Briet et al (2012), supplementary information
# For LLINs PN2 gambiae Akron extract the 5th row 

Pparams <- read.csv(file=paste0(ds,"Permanet_parameters.csv"),sep=";")
params <- Pparams[5,]

# From eq1 in Briet et al, supplementary information
# p is the insecticide concentration and h is the holed area. 
preprandialKillingEffect<- function(params,p,h){
  effect <- with(params,preprandialKillingEffect_baseFactor + 
          preprandialKillingEffect_holeFactor*exp(-h*preprandialKillingEffect_holeScalingFactor) +
          preprandialKillingEffect_insecticideFactor*(1-exp(-p* preprandialKillingEffect_insecticideScalingFactor))+ 
          preprandialKillingEffect_interactionFactor*exp(-h* preprandialKillingEffect_holeScalingFactor)*
          (1-exp(-p*preprandialKillingEffect_insecticideScalingFactor)))
  return(effect)
}

# From eq5 in Briet et al, supplement
enteringValue<-function(params,p){
  effect <- with(params,exp(log(entering_insecticideFactor)*(1-exp(-p*entering_insecticideScalingFactor))))
  return(effect)                                                             
}
                                 
# From eq6 in Briet et al, supplement
attackingValue<-function(params,p,h){
  effect <- with(params,attacking_baseFactor + 
                   attacking_holeFactor*exp(-h*attacking_holeScalingFactor) +
                   attacking_insecticideFactor*(1-exp(-p* attacking_insecticideScalingFactor))+ 
                   attacking_interactionFactor*exp(-h* attacking_holeScalingFactor)*
                   (1-exp(-p*attacking_insecticideScalingFactor)))
  return(effect)
}

# I don't find the equation for post-prandial killing, but the statement that the holes don't affect post-prandial killing
# and the available parameters suggest that the following is intended (even though one would expect that 
# holes should reduce post-prandial killing)
postprandialKillingEffect<-function(params,p){
  effect <- with(params,postprandialKillingEffect_baseFactor + 
                   postprandialKillingEffect_insecticideFactor*(1-exp(-p* postprandialKillingEffect_insecticideScalingFactor))+ 
                   (1-exp(-p*postprandialKillingEffect_insecticideScalingFactor)))
  return(effect)
}

# Parameter estimates for the Briet et al (2018) model
scalingfactor<-1
beta0.PEnt <- logit(0.99) 
beta1.PEnt <- -0.986
beta0.PAtt <- -4.069
beta1.PAtt <- 0.541
beta2.PAtt <- 0.807
beta3.PAtt <- -0.101		
beta0.PBmu <- -4.840
beta1.PBmu <- -0.128
beta2.PBmu <- 2.325
beta3.PBmu <- -0.113		
beta0.PCmu <- logit(0.01)
beta1.PCmu <- -0.114
beta2.PCmu <- 1.303
beta3.PCmu <- 0
logHolesMax<-log(192000+1)



#### set number of interpolation points and create empty vectors for results.
Pii <- 0.89
nips <- num.interpolation.points <- 20
alphai.reduction <- PBmu.pred <- PCmu.pred <- PEnt.pred <- PAtt.pred <- net.survival <- rep(NA,nips)
PBi <- PCi <- entering <- attacking <- preprandialKilling <- postprandialKilling <- rep(NA,nips)

for (ip in 1:nips) {
  
# For LLINs PN2 gambiae Akron with decay in holed surface area and insecticide and via functions that translate physical status into effect. Attrition is also modelled, and affects the coverage via the survival parameter (in multiple_feval_impact.r)

  insecticide.content <- Effect_exponentialDecay(initialEffect = 55, halflife = 1.5, duration = 3, nips, ip) #Central scenario in https://static-content.springer.com/esm/art%3A10.1186%2F1475-2875-12-401/MediaObjects/12936_2013_3728_MOESM2_ESM.xml
  net.survival[ip] <- Effect_smoothcompactDecay(initialEffect = 1, L = 20.7725, exponent = 18, duration = 3, nips, ip) #Central scenario in  https://static-content.springer.com/esm/art%3A10.1186%2F1475-2875-12-401/MediaObjects/12936_2013_3728_MOESM2_ESM.xml
  #<holeRate mean="1.80" sigma="0.80"/><ripRate mean="1.80" sigma="0.80"/><ripFactor value="0.3"/>
  holed.surface.area <- Effect_parabolicGrowth(a=0, b= 58.5949, duration = 3, nips, ip)
  #This approximates the central scenario hole index. It probably needs updating to holed area in cm2. 

  # For Olivier's 2012 parameterisation  
  
  entering[ip]<- enteringValue(params,insecticide.content)
  attacking[ip] <- attackingValue(params,insecticide.content,holed.surface.area)
  preprandialKilling[ip] <- preprandialKillingEffect(params,insecticide.content,holed.surface.area)  
  postprandialKilling[ip] <- postprandialKillingEffect(params,insecticide.content)  

  # For Olivier's 2018 parameterisation  

  holed.surface.area.cap <- min(holed.surface.area, 192000)
  logHoles <- log(holed.surface.area.cap + 1)
  
  PEnt.u<-ilogit(beta0.PEnt + scalingfactor*beta1.PEnt*log(0+1))		
  PAtt.u<-ilogit(beta0.PAtt + beta1.PAtt*logHolesMax +scalingfactor*beta2.PAtt*log(0+1) + beta3.PAtt*scalingfactor*logHolesMax*log(0+1))
  PBmu.u<-ilogit(beta0.PBmu + beta1.PBmu*logHolesMax +scalingfactor*beta2.PBmu*log(0+1) + beta3.PBmu*scalingfactor*logHolesMax*log(0+1))
  PCmu.u<-ilogit(beta0.PCmu + beta1.PCmu*logHolesMax +scalingfactor*beta2.PCmu*log(0+1) + beta3.PCmu*scalingfactor*logHolesMax*log(0+1))
  
  PEnt<-ilogit(beta0.PEnt + scalingfactor*beta1.PEnt*log(insecticide.content+1))		
  PAtt<-ilogit(beta0.PAtt + beta1.PAtt*logHoles +scalingfactor*beta2.PAtt*log(insecticide.content+1) + beta3.PAtt*scalingfactor*logHoles*log(insecticide.content+1))
  PBmu<-ilogit(beta0.PBmu + beta1.PBmu*logHoles +scalingfactor*beta2.PBmu*log(insecticide.content+1) + beta3.PBmu*scalingfactor*logHoles*log(insecticide.content+1))
  PCmu<-ilogit(min(600,beta0.PCmu + beta1.PCmu*logHoles +scalingfactor*beta2.PCmu*log(insecticide.content+1) + beta3.PCmu*scalingfactor*logHoles*log(insecticide.content+1)))
  
  PEnt.pred[ip] <- PEnt/PEnt.u
  PAtt.pred[ip] <- PAtt/PAtt.u
  alphai.reduction[ip] <- 1-PEnt.pred[ip]*PAtt.pred[ip]
  PBmu.pred[ip] <- 1-(1-PBmu)/(1-PBmu.u)	
  PCmu.pred[ip] <- 1-(1-PCmu)/(1-PCmu.u)
  PBi[ip]<- 0.95*(1-Pii*PBmu.pred[ip])
  PCi[ip]<- 0.95*(1-Pii*PCmu.pred[ip])
}  


#params[['LLINs PN2 gambiae Akron']][['alphai']][1] <- alphai[2]*(1 - Pii*alphai.reduction) 
#params[['LLINs PN2 gambiae Akron']][['PBi']][1] <- PBi*(1-Pii*PBmu.pred) 
#params[['LLINs PN2 gambiae Akron']][['PCi']][1] <- PCi*(1-Pii*PCmu.pred)
#params[['LLINs PN2 gambiae Akron']][['survival']] <- net.survival

library(ggplot2)

product2018 <- PEnt.pred * PAtt.pred * PBi * PCi
product2012 <- entering * attacking * (1 - preprandialKilling) * (1-postprandialKilling)
df <- data.frame(alphai.reduction,PBi,PCi,PEnt.pred,PAtt.pred,product2012,product2018,
                 net.survival,entering, attacking,preprandialKilling,postprandialKilling,interpolation=seq(1:20))
ggplot(data=df, aes(x=interpolation))+
  theme_bw(base_size = 14)+
  scale_y_continuous(limits=c(0,1))+
  geom_line(size=1,y=PAtt.pred, colour='blue')+
  geom_line(size=1,y=PEnt.pred, colour='orange')+
  geom_line(size=1,y=PCi, colour='red')+
  geom_line(size=1,y=PBi, colour='green')+
  geom_line(size=1,y=attacking, colour='blue',linetype = "dashed")+
  geom_line(size=1,y=entering, colour='orange',linetype = "dashed")+
  geom_line(size=1,y=1-postprandialKilling, colour='red',linetype = "dashed")+
  geom_line(size=1,y=1-preprandialKilling, colour='green',linetype = "dashed")+
  geom_line(size=1,y=product2018, colour='black') +
geom_line(size=1,y=product2012, colour='black',linetype = "dashed")

