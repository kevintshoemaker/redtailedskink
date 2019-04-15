#####################
# Western red-tailed skink project
#####################

# collaboration with Derek Hall

# total of 168 sites surveyed. [154 actually in the final dataset]
# total of 29 sites with skinks
# no covariates available for detection prob  [check on this!]
# covars: soils, elevation, slope, veg, moisture, geology
# Encounter history for 6 days of trapping at 120 sites
# 12 sites were sampled in multiple years
# Sometimes traps were checked every 2 or 3 days, the dash means traps not checked that day (but the trap was open!)
    

# Research Question: what environmental factors are associated with skink occupancy? (account for imperfect detection)



## TODO: run zero-inflated poisson model with offset to do model selection!
## try quadratic term on elevation -- middle elevations were preferred!  (see derek's paper in herp review)
## try interaction of aspect and moisture!

# TODO: try adding vegetation back into the final model
# TODO: try making elevation a threshold relationship instead of quadratic... [actually the quadratic works pretty well!]



####################
# Questions for Derek
####################

# Any covariates potentially related to detection probability? If so, which ones?
# What about interaction terms?

# If a site is occupied in one year, is it sure to be occupied the next year?  [YES]

# what does the "skink" covariate signify? [NOTE: there are some historical sites that we couldn't find them in this study --# skink category: yes means there was red tailed skink.  also has historic (not captured in this study). SKILT is great basin skink.]

# don't worry about vegetation?

# veg association- run Chi Squared test
# geology - run chisquare test

# all skinks occurred in volcanic...

# ecoregion.   Important. Tied to elevation.

# combine slope and aspect into a single variable? South-facingness? 



####################
# Notes
####################

# moisture: wet vs dry (mesic is wettish)
# if a site is occupied, it's occupied at other surveys. 

# try running ecoregion as originally specified
#    site 103 is transition
#    site 108 is Great Basin

# Kevin,
# 
# Thinking about skink ecology based on the literature and what I observed during the study, 
# western red-tailed skinks occur within certain elevation ranges, they like springs and mesic areas, 
# and they need cover (rocks, cracks, litter). Thus; I believe elevation, aspect, moisture and substrate 
# and to a lesser degree vegetation (because it is tied in with elevation and aspect and moisture) 
# are the most important variables. Does it make sense to look at all possible combinations of these 5 variables 
# (32 potential models with all interactions)? I suggest looking at all main effects, then the 
# 5 main effects I mentioned. I need some help trying to understand what other interactions 
# mean and which ones we should include. Maybe call me tomorrow or Thursday to discuss. I should be around both days. Thanks!
  
  
#####################
# Prepare workspace
#####################

rm(list=ls())

library(unmarked)
library(R2jags)
library(Hmisc)
library(MASS)
library(boot)
library(vcd)
library(AER)
library(car)
library(DHARMa)
library(pscl)    # for zero inflated model?

####################
# Load custom functions
####################


chisq <- function(fm) {
  umf <- getData(fm)
  y <- getY(umf)
  y[y>1] <- 1
  sr <- fm@sitesRemoved
  if(length(sr)>0)
    y <- y[-sr,,drop=FALSE]
  fv <- fitted(fm, na.rm=TRUE)
  y[is.na(fv)] <- NA
  sum((y-fv)^2/(fv*(1-fv)), na.rm=TRUE)
}

var = "ecoregion"
var = "elev"
predict.jags <- function(var,quad=F){
  namestable2 <- namestable
  lots = 1000
  quadvars = NULL
  allsamples <- mod$BUGSoutput$sims.list
  totsamples <- nrow(allsamples[[1]])
  
  origname <- namestable2$name2[namestable2$name1==var]
  notthisvar <- namestable2$name2[namestable2$name1!=var]
  
  varclasses <- sapply(1:ncol(sitecovs),function(t) class(sitecovs[,t]))
  names(varclasses) <- names(sitecovs)
  
  betanames <- character(nrow(namestable2))
  
  i=1
  for(i in 1:nrow(namestable2)){
    temp <- names(allsamples[grep(namestable2$name1[i],names(allsamples))])
    if(length(temp)>1){
      quadvars=c(quadvars,namestable2$name2[i])
      betanames <- c(betanames,"")
      namestable2 <- rbind(namestable2,paste(namestable2[i,],"2",sep=""))
      varclasses <- c(varclasses,"numeric")
      names(varclasses)[length(varclasses)] <- sprintf("%s2",namestable2$name2[i])
      betanames[c(i,length(betanames))] <- temp
    }else{
      betanames[i] <- temp
    }
  }
  
  t=4
  funct <- function(t){
    if(varclasses[t]=="factor"){
      tab <- table(sitecovs[,t])
      factor(names(tab[which.max(tab)]),levels=levels(sitecovs[,t]))
    }else{
      mean(sitecovs[,t],na.rm=T)
    }
  }
  
  allmeans <- lapply(1:ncol(sitecovs),function(t) funct(t))
  names(allmeans) <- names(sitecovs)
  
  isfact <- is.factor(sitecovs[,origname])
  if(isfact){
    newdata <- data.frame(
      var1 = factor(levels(sitecovs[,origname]),levels=levels(sitecovs[,origname]))
    )
    v=1
    for(v in 1:length(notthisvar)){
      if(notthisvar[v]%in%quadvars){
        newdata[,notthisvar[v]] = ifelse(varclasses[notthisvar[v]]!="factor",0, allmeans[notthisvar[v]][[1]])
        newdata[,sprintf("%s2",notthisvar[v])] = 0 #allmeans[notthisvar[v]][[1]]^2
      }else{
        newdata[,notthisvar[v]] = ifelse(varclasses[notthisvar[v]]!="factor",0, allmeans[notthisvar[v]][[1]])
      }
    }
  }else{
    if(origname%in%quadvars){
      varrang <- range(scale(sitecovs[,origname]))
      newdata <- data.frame(
        var1 = seq(varrang[1],varrang[2],length=100),
        var2 = seq(varrang[1],varrang[2],length=100)^2
      )
      for(v in 1:length(notthisvar)){
        if(notthisvar[v]%in%quadvars){
          newdata[,notthisvar[v]] = allmeans[notthisvar[v]]
          newdata[,sprintf("%s2",notthisvar[v])] = 0
        }else{
          newdata[,notthisvar[v]] = allmeans[notthisvar[v]]
        }
      }
    }else{
      varrang <- range(scale(sitecovs[,origname]))
      newdata <- data.frame(
        var1 = seq(varrang[1],varrang[2],length=100)
      )
      for(v in 1:length(notthisvar)){
        if(notthisvar[v]%in%quadvars){
          newdata[,notthisvar[v]] = allmeans[notthisvar[v]]
          newdata[,sprintf("%s2",notthisvar[v])] = 0
        }else{
          newdata[,notthisvar[v]] = allmeans[notthisvar[v]]
        }
      }
    }

  }
  names(newdata)[1] <- origname
  if(origname%in%quadvars) names(newdata)[2] <- sprintf("%s2",origname)
  
  newdata <- newdata[,match(namestable2$name2,names(newdata))]
  
  sims <- matrix(NA,ncol=nrow(newdata),nrow=lots)
  
  # logit(p.occupied[site]) <- baseocc.logit + beta.elev.o*elev[site] + beta.slope.o[slope[site]] +
  #   beta.substrate.o[substrate[site]] + beta.landform.o[landform[site]] + 
  #   beta.vegassoc.o[vegassoc[site]] + beta.ecoregion.o[ecoregion[site]] +
  #   beta.moisture.o[moisture[site]]
  i=1
  for(i in 1:nrow(sims)){
    thissamp <- sample(1:totsamples,1)
    j=1
    for(j in 1:ncol(sims)){
      linear.occ <- qlogis(allsamples[["baseocc"]][thissamp])
      k=1
      for(k in 1:length(betanames)){
        if(varclasses[namestable2$name2[k]]=="factor"){
          linear.occ <- linear.occ + allsamples[[betanames[k]]][thissamp,as.numeric(newdata[j,namestable2$name2[k]])]
        }else{
          #sc <- scale(sitecovs[,namestable2$name2[k]])
          #conv <- (newdata[j,namestable2$name2[k]]-attributes(sc)[["scaled:center"]])/attributes(sc)[["scaled:scale"]]
          linear.occ <- linear.occ + allsamples[[betanames[k]]][thissamp]*newdata[j,namestable2$name2[k]]
        }
      }
      sims[i,j] <- plogis(linear.occ)
      
    }
  }  
  
  quantiles <- sapply(1:ncol(sims), function(t) quantile(sims[,t],c(0,0.025,0.25,0.5,0.75,0.975,1))  )
  colnames(quantiles) <- newdata[,origname]
  
  filename <- sprintf("bivar_plot_%s.svg",var)
  
  #png(filename=filename,width=500,height=500)
  svg(filename=filename,width=4,height=4)
  
  if(isfact){
    xvals <- barplot(quantiles["50%",],ylim=c(0,1),ylab="P(occurrence)")
    errbar(xvals,quantiles["50%",],quantiles["97.5%",],quantiles["2.5%",],add=T)
  } else{
    sc <- scale(sitecovs[,origname])   # attributes(sc)[["scaled:center"]])/attributes(sc)[["scaled:scale"]]
    frac = diff(range(newdata[,origname]))/30
    plot(newdata[,origname],apply(quantiles,2,median),ylim=c(0,1),ylab="P(occurrence)",xlab=origname,type="l",lwd=3,
         xlim=range(newdata[,origname])+c(frac,-frac),xaxt="n")
    polygon(c(newdata[,origname],rev(newdata[,origname])),c(quantiles["2.5%",],rev(quantiles["97.5%",])),col=gray(0.7),border=NA)
    lines(newdata[,origname],apply(quantiles,2,median),lwd=3)
    temp <- seq(-2,2,length=5)
    axis(1,at=temp,labels=round((temp*attributes(sc)[["scaled:scale"]])+attributes(sc)[["scaled:center"]]))
    rug(sc)
  }
  
  
  dev.off()
  ### summarize spread and visualize effect size...
}



#####################
# Load data
#####################

sites <- read.csv("sites.csv",header=T,stringsAsFactors = F)
head(sites)

sites2 <- read.csv("sites.csv",header=T)
head(sites2)
summary(sites2)     # summarize factors

detections <- read.csv("detections.csv",na.strings = c("-","",NA),header=T,stringsAsFactors = F)
head(detections,10)

#####################
# Process data
#####################

caphist <- as.matrix(detections[,grepl("Visit",names(detections))])     # detection histories
caphist <- data.matrix(as.data.frame(gsub("[*]","0",caphist),stringsAsFactors = F))
rownames(caphist) <- detections$LocationCode      # sites corresponding to detection histories

effort <- detections$TrapNights
names(effort) <- detections$LocationCode

captures_total <- apply(caphist,1,sum,na.rm=T)
captures_bin <- ifelse(captures_total>0,1,0)

sum(captures_bin)      # 27 known-occupied sites

sum(captures_total)    # 43 total captures
 
sitelist <- unique(detections$LocationCode)       # unique surveyed sites


############
# process site-level data

names(sites)

sitecovs <- sites[,c("LocationCode","Elevation","X.Slope","Geology.Type","Aspect",
                     "NewSoilSubstrate","NEWLandform","NEW.Vegetation","Ecoregion",
                     "NEWMoisture")]

names(sitecovs) <- c("LocationCode","Elevation","Slope","Geology","Aspect",
                     "Substrate","Landform","Vegetation","Ecoregion",
                     "Moisture")

rownames(sitecovs) <- sitecovs$LocationCode

# sitecovs$LocationCode %in% sitelist    # check: make sure all relevant sites are in the site list
# sitelist %in% sitecovs$LocationCode


hist(sitecovs$Elevation)    # rescale this!
sitecovs$Elevation_sc <- scale(sitecovs$Elevation)

sitecovs$Slope   # doesn't seem to be useful- too many "varied" observations. Maybe call these NA for now? or make categorical? 
sitecovs$Slope[sitecovs$Slope%in%as.character(0:10)] <- "flat"
sitecovs$Slope[sitecovs$Slope%in%as.character(11:40)] <- "moderate"
sitecovs$Slope[sitecovs$Slope%in%as.character(41:99)] <- "steep"
sitecovs$Slope <- gsub("0-10","flat",sitecovs$Slope)
table(sitecovs$Slope)
sitecovs$Slope <- factor(sitecovs$Slope,levels=c("Varied","flat","moderate","steep"))

sitecovs$Aspect
sitecovs$Aspect[sitecovs$Aspect%in%as.character(100:270)] <- "south"
sitecovs$Aspect[sitecovs$Aspect%in%as.character(0:90)] <- "north"
sitecovs$Aspect[sitecovs$Aspect%in%as.character(271:360)] <- "north"
sitecovs$Aspect <- gsub("Varied SE to SW","south",sitecovs$Aspect)
sitecovs$Aspect <- gsub("50-120","Varied",sitecovs$Aspect)
table(sitecovs$Aspect)
sitecovs$Aspect <- factor(sitecovs$Aspect,levels=c("Varied","south","north"))

sitecovs$Substrate
# unique(sitecovs$Substrate)
table(sitecovs$Substrate)
# temp <- sapply(1:nrow(sitecovs),function(t) strsplit(sitecovs$Substrate[t],"/") )
# temp2 <- t(sapply(1:nrow(sitecovs),function(t) temp[[t]][1:2] ))
# temp3 <- sapply(1:nrow(sitecovs),function(t) paste(temp2[t,],collapse="/") )
# sitecovs$Substrate[grepl("Boulder",temp3)] <- "Rock"
# sitecovs$Substrate[grepl("Rock",temp3)] <- "Rock"
# sitecovs$Substrate[grepl("Gravel",temp3)] <- "Gravel"
# sitecovs$Substrate[grepl("Sand",temp3)] <- "Gravel"
# sitecovs$Substrate[grepl("Loam",temp3)] <- "Loam"
# sitecovs$Substrate[grepl("Litter",temp3)] <- "Litter"
# table(sitecovs$SoilSubstrate)
sitecovs$Substrate <- factor(sitecovs$Substrate,levels=c("Litter","Soil","Rock","Boulder"))


sitecovs$Landform
#unique(sitecovs$Landform)
table(sitecovs$Landform)
# temp <- sapply(1:nrow(sitecovs),function(t) strsplit(sitecovs$Landform[t],"/") )
# temp2 <- t(sapply(1:nrow(sitecovs),function(t) temp[[t]][1:2] ))
# temp3 <- sapply(1:nrow(sitecovs),function(t) paste(temp2[t,],collapse="/") )
# sitecovs$Landform[grepl("Hill",temp3)] <- "Hill"
# sitecovs$Landform[grepl("Mesa",temp3)] <- "Hill"
# sitecovs$Landform[grepl("Hilltop",temp3)] <- "Hill"
# sitecovs$Landform[grepl("Knoll",temp3)] <- "Hill"
# sitecovs$Landform[grepl("slope",temp3)] <- "Hill"
# sitecovs$Landform[grepl("slope",temp3)] <- "Hill"
# sitecovs$Landform[grepl("Mesa",temp3)] <- "Hill"
# sitecovs$Landform[grepl("Mountain",temp3)] <- "Hill"
# sitecovs$Landform[grepl("outcrop",temp3)] <- "Hill"
# sitecovs$Landform[grepl("Bajada",temp3)] <- "Valley"
# sitecovs$Landform[grepl("Wash",temp3)] <- "Valley"
# sitecovs$Landform[grepl("Canyon",temp3)] <- "Valley"
# sitecovs$Landform[grepl("valley",temp3)] <- "Valley"
# sitecovs$Landform[grepl("Valley",temp3)] <- "Valley"
# sitecovs$Landform[grepl("Draw",temp3)] <- "Valley"
# table(sitecovs$Landform)
sitecovs$Landform <- factor(sitecovs$Landform,levels=c("Slope","Drainage","Hill/Knoll","Ridge","Cliff"))



sitecovs$Vegetation
#unique(sitecovs$Vegetation)
table(sitecovs$Vegetation)
# temp <- sapply(1:nrow(sitecovs),function(t) strsplit(sitecovs$Vegetation[t],"\\/|\\-") )
# temp2 <- t(sapply(1:nrow(sitecovs),function(t) temp[[t]][1:2] ))
# temp3 <- sapply(1:nrow(sitecovs),function(t) paste(temp2[t,],collapse="/") )
# sitecovs$Vegetation[grepl("Woodland",temp3)] <- "Woodland"
# sitecovs$Vegetation[grepl("Shrubland",temp3)] <- "Shrubland"
# sitecovs$Vegetation[grepl("Miscellaneous",temp3)] <- "Other"
# sitecovs$Vegetation[grepl("Disturbed",temp3)] <- "Other"
# table(sitecovs$Vegetation)
sitecovs$Vegetation <- factor(sitecovs$Vegetation,levels=c("Woodland","Sagebrush","Blackbrush","Nevada Ephedra","White bursage"))
levels(sitecovs$Vegetation) <- c("Woodland","Sagebrush","Blackbrush","Ephedra","Bursage")

sitecovs$Ecoregion
# sitecovs$Ecoregion[sitecovs$Ecoregion=="Great Basin/Transition"] <- "Transition/Mojave"
# sitecovs$Ecoregion[sitecovs$Ecoregion=="Transition"] <- "Transition/Mojave"
# sitecovs$Ecoregion[sitecovs$Ecoregion=="Mojave"] <- "Transition/Mojave"
table(sitecovs$Ecoregion)
sitecovs$Ecoregion <- as.factor(sitecovs$Ecoregion)

sitecovs$Geology
# sitecovs$Ecoregion[sitecovs$Ecoregion=="Great Basin/Transition"] <- "Transition/Mojave"
# sitecovs$Ecoregion[sitecovs$Ecoregion=="Transition"] <- "Transition/Mojave"
# sitecovs$Ecoregion[sitecovs$Ecoregion=="Mojave"] <- "Transition/Mojave"
table(sitecovs$Geology)
sitecovs$Geology <- factor(sitecovs$Geology,levels=c("Volcanic","Sedimentary","Metamorphic"))

sitecovs$Moisture
#unique(sitecovs$Moisture)
table(sitecovs$Moisture)
# temp <- sapply(1:nrow(sitecovs),function(t) strsplit(sitecovs$Moisture[t],"/") )
# temp2 <- t(sapply(1:nrow(sitecovs),function(t) temp[[t]][1:2] ))
# temp3 <- sapply(1:nrow(sitecovs),function(t) paste(temp2[t,],collapse="/") )
# sitecovs$Moisture[grepl("Dry",temp3)] <- "Dry"
# sitecovs$Moisture[grepl("Mesic",temp3)] <- "Mesic"
# sitecovs$Moisture[grepl("Wet",temp3)] <- "Mesic"
# table(sitecovs$Moisture)
sitecovs$Moisture <- as.factor(sitecovs$Moisture)

# rename to conform with the BUGS model

# elev[site],slope[site],substrate[site], landform[site], vegassoc[site], ecoregion[site], moisture[site]

namestable <- data.frame(
  name1 = c(
    "elev",
    #"slope",
    "aspect",
    "substrate", 
    #"landform", 
    "vegassoc", 
    "ecoregion", 
    "moisture"
    ),
  name2 = c(
    "Elevation",
    #"Slope",
    "Aspect",
    "Substrate",
    #"Landform",
    "Vegetation",
    "Ecoregion",
    "Moisture"
    ),
  stringsAsFactors = F
)

##############
# Exploratory visualizations
##############

summary(sitecovs)

table(sitecovs$Ecoregion,sitecovs$Moisture)    # moisture is highly correlated with ecoregion  - GB is wetter

table(sitecovs$Ecoregion,sitecovs$Vegetation)   # GB has woodland and sagebrush, Mojave has blackbrush and bursage. transition is ephedra and blackbrush

table(sitecovs$Landform,sitecovs$Vegetation)     # landform not highly correlated, but may not be useful?

table(sitecovs$Vegetation,sitecovs$Substrate)    # not super correlated.

table(sitecovs$Slope,sitecovs$Aspect)

plot(sitecovs$Elevation~sitecovs$Slope)         # no correlation between slope and elevation

plot(sitecovs$Elevation~sitecovs$Ecoregion)    # ecoregion is highly correlated with elevation




# numeric covars: Elevation
# factor covars: Slope (3), Aspect (3), Substrate (3), Landform (2), Vegetation (3), Ecoregion (3), Moisture (2)


###################
# Explopratory analyses: (known) occupancy vs covariates..
####################

Occupied <- apply(caphist,1,function(t) ifelse(sum(t,na.rm=T)>0,1,0)  )


ndx <- match(names(Occupied),rownames(sitecovs))

# Slope
table(Occupied,sitecovs$Slope[ndx])    # more occupancy on flat and varied?..
chisq.test(Occupied,sitecovs$Slope[ndx])
fisher.test(Occupied,sitecovs$Slope[ndx])    # no relationship with slope

table(Occupied,sitecovs$Aspect[ndx])      # strong relation with aspect - south aspects strongly associated with occupancy
chisq.test(Occupied,sitecovs$Aspect[ndx])
fisher.test(Occupied,sitecovs$Aspect[ndx]) 

table(Occupied, sitecovs$Substrate[ndx])     # no relation with substrate
chisq.test(Occupied, sitecovs$Substrate[ndx])
fisher.test(Occupied, sitecovs$Substrate[ndx])

table(Occupied, sitecovs$Landform[ndx])     # no relation with landform
chisq.test(Occupied, sitecovs$Landform[ndx])
fisher.test(Occupied, sitecovs$Landform[ndx])

table(Occupied, sitecovs$Geology[ndx])     # no relation with Geology
chisq.test(Occupied, sitecovs$Geology[ndx])
fisher.test(Occupied, sitecovs$Geology[ndx])

table(Occupied, sitecovs$Vegetation[ndx])     # weak relation with veg assoc - more associated with sagebrush
chisq.test(Occupied, sitecovs$Vegetation[ndx])
fisher.test(Occupied, sitecovs$Vegetation[ndx])   

table(Occupied, sitecovs$Ecoregion[ndx])     # strong relationship with ecoregion - more occupancy in great basin
chisq.test(Occupied, sitecovs$Ecoregion[ndx])
fisher.test(Occupied, sitecovs$Ecoregion[ndx])

table(Occupied, sitecovs$Moisture[ndx])     # strong relationship with moisture
chisq.test(Occupied, sitecovs$Moisture[ndx])
fisher.test(Occupied, sitecovs$Moisture[ndx])


## fairly strong relationships with moisture, aspect, ecoregion, vegetation

# test the linker.  
#cbind(as.numeric(names(Occupied)),as.numeric(site_session_linker$sitenames1))

plot(sitecovs$Elevation[ndx]~as.factor(Occupied))           # no relationship with Elevation
summary(glm(Occupied ~ sitecovs$Elevation[ndx], family="binomial"))


plot(sitecovs$Elevation_sc[ndx],Occupied)
mod1 <- glm(Occupied ~ sitecovs$Elevation_sc[ndx] + I(sitecovs$Elevation_sc[ndx]^2), family="binomial")
summary(mod1)
cf <- coefficients(mod1)
curve(plogis(cf[1]+cf[2]*x+cf[3]*x^2),add=T)

###################
# BUGS code
###################

Bugsfile <- "Skinkcode_bugs.txt"
cat(
  "
  model{

##################  
    ### likelihood

  for(site in 1:nsites){
      log.effort[site] <- log(effort[site])
      log(expcount[site]) <- log.effort[site] + basecount.log         # for now, assume equal capture rate across all occupied sites
      expcount2[site] <- occupied[site] * expcount[site]              # zero inflation
      p[site] <- r/(r+expcount2[site])
      obscount[site] ~ dnegbin(p[site],r)   #dpois(expcount2[site])                         # data node
      obscount2[site] ~ dnegbin(p[site],r)  # dpois(expcount2[site])                        # for goodness of fit testing
      error.obs[site] <- (obscount[site]-expcount2[site])^2 
      error.sim[site] <- (obscount2[site]-expcount2[site])^2
  }
  
  toterr.obs <- sum(error.obs[])    # SSEobs
  toterr.sim <- sum(error.sim[])    # SSEsim

    ### occupancy
  for(site in 1:nsites){  
    logit(p.occupied[site]) <- baseocc.logit + 
                                  beta.elev.o*elev[site] + 
                                  beta.elev2.o*pow(elev[site],2) +
                                  beta.moisture.o[moisture[site]] +
                                   #beta.slope.o[slope[site]] +
                                  beta.aspect.o[aspect[site]] +
                                  beta.substrate.o[substrate[site]] + 
                                   #beta.landform.o[landform[site]] + 
                                  beta.vegassoc.o[vegassoc[site]] + 
                                  beta.ecoregion.o[ecoregion[site]] 
    occupied[site] ~ dbern(p.occupied[site])
  }


    ### priors
  baseocc ~ dunif(0,1)
  r ~ dunif(0.1,20)
  basecount ~ dunif(0,2)     # baseline expected count for an occupied site
  baseocc.logit <- log(baseocc/(1-baseocc))
  basecount.log <- log(basecount)

  beta.elev.o ~ dnorm(0,0.5)
  beta.elev2.o ~ dnorm(0,0.5)


  beta.aspect.o[1] <- 0
  for(i in 2:naspects){
    beta.aspect.o[i] ~ dnorm(0,0.5)
  }

  beta.substrate.o[1] <- 0
  for(i in 2:nsubstrates){
    beta.substrate.o[i] ~ dnorm(0,0.5)
  }

  # beta.landform.o[1] <- 0
  # for(i in 2:nlandforms){
  #   beta.landform.o[i] ~ dnorm(0,0.5)
  # }

  beta.vegassoc.o[1] <- 0
  for(i in 2:nvegassocs){
    beta.vegassoc.o[i] ~ dnorm(0,0.5)
  }

  beta.ecoregion.o[1] <- 0
  for(i in 2:necoregions){
    beta.ecoregion.o[i] ~ dnorm(0,0.5)
  }

  beta.moisture.o[1] <- 0
  for(i in 2:nmoistures){
    beta.moisture.o[i] ~ dnorm(0,0.5)
  }

}
  
  "
  ,file=Bugsfile
)



#############
# Prepare data for JAGS
#############

data.for.jags <- list(
  nsites = nrow(sitecovs),
  nsubstrates = length(levels(sitecovs$Substrate)),
  #nslopes = length(levels(sitecovs$Slope)),
  naspects = length(levels(sitecovs$Aspect)),
  #nlandforms = length(levels(sitecovs$Landform)),
  nvegassocs = length(levels(sitecovs$Vegetation)),
  necoregions = length(levels(sitecovs$Ecoregion)),
  nmoistures = length(levels(sitecovs$Moisture)),
  elev = as.numeric(scale(sitecovs$Elevation[ndx])),
  #slope = as.numeric(sitecovs$Slope[ndx]),
  aspect = as.numeric(sitecovs$Aspect[ndx]),
  substrate = as.numeric(sitecovs$Substrate[ndx]),
  #landform = as.numeric(sitecovs$Landform[site_session_linker$index]),
  vegassoc = as.numeric(sitecovs$Vegetation[ndx]),
  ecoregion = as.numeric(sitecovs$Ecoregion[ndx]),
  moisture = as.numeric(sitecovs$Moisture[ndx]),
  obscount = captures_total,
  effort = effort/180
)


inits.for.jags <- function(){
  list(
    occupied = rep(1,times=nrow(caphist)),   # sitecovs
    baseocc = 0.5,
    basecount = 0.005,
    r = 10
  )
}

inits.for.jags()


# monitor: baseocc, basedetect, beta.moisture.o, beta.ecoregion.o, beta.vegassoc.o, beta.landform.o, beta.substrate.o,
#  beta.substrate.o, beta.slope.o, beta.elev.o, occupied

params.to.monitor <- c(
  "baseocc", 
  "basecount", 
  "beta.moisture.o", 
  "beta.ecoregion.o", 
  "beta.vegassoc.o", 
  #"beta.landform.o", 
  "beta.substrate.o",
  #"beta.slope.o",
  "beta.aspect.o",
  "beta.elev.o",
  "beta.elev2.o",
  "occupied",
  "toterr.obs",
  "toterr.sim",
  "obscount2",
  "r"
)

niter = 40000
mod <- jags(data=data.for.jags, parameters.to.save = params.to.monitor, inits = inits.for.jags, 
            model.file = Bugsfile, n.iter=niter, DIC=FALSE, n.chains = 3,
            n.burnin = 20000, n.thin = 10)


##############
# Assess convergence
##############

mod2 <- as.mcmc(mod)

#str(mod2)

plot(mod2[,"basecount"],main=c("Baseline expected count"))
plot(mod2[,"baseocc"], main="Baseline occupancy")
#plot(mod2[,"beta.elev.o"])

gelman.diag(mod2[,"basecount"])
gelman.diag(mod2[,"baseocc"])    # good convergence

gelman.diag(mod2[,"beta.ecoregion.o[2]"])
gelman.diag(mod2[,"beta.elev2.o"])
gelman.diag(mod2[,"beta.elev.o"])
gelman.diag(mod2[,"beta.moisture.o[2]"])
gelman.diag(mod2[,"beta.aspect.o[2]"])
gelman.diag(mod2[,"beta.vegassoc.o[2]"])
gelman.diag(mod2[,"r"])

#############
# Goodness of fit tests (posterior predictive check)
#############

# Bayesian p value and plot
# goodness of fit plot- number of zeros, ones, twos etc... 

SSEobs <- mod$BUGSoutput$sims.list$toterr.obs
SSEsim <- mod$BUGSoutput$sims.list$toterr.sim
pval <- length(which(SSEsim > SSEobs))/length(SSEobs)

png("GOF1.png",width=500,height = 400)
plot(SSEsim~SSEobs,ylim=c(0,250),xlim=c(0,100))
abline(1,1,col="blue",lwd=2)
text(20,200,sprintf("Bayesian p-value: %s",pval))
dev.off()


nsites <- nrow(sites)
nsims <- nrow(simobs)
simobs <- mod$BUGSoutput$sims.list$obscount2[,1:nsites]

gof2 <- t(sapply(1:nsims,function(v) sapply(0:6,function(t) length(which(simobs[v,]==t)) )  ))

colnames(gof2) <- c(0:6)

head(gof2)

gof2 <- tidyr::gather(as.data.frame(gof2),key="count")
head(gof2)

gof2_obs <- sapply(0:6,function(t) length(which(captures_total==t)) )

png("GOF2.png",width=500,height = 400)
boxplot(gof2$value~gof2$count,xlab="Total skink observations",ylab="Number of surveyed sites")
text((0:6)+1,gof2_obs,"X",cex=1.5,col="blue")
dev.off()

###### final goodness of fit figure

svg("GOF_final.svg",w=5,h=5)
layout(matrix(c(1,2),nrow=2))
par(mai=c(0.9,0.9,0.1,0))
plot(SSEsim~SSEobs,ylim=c(0,250),xlim=c(0,100),xlab="SSE, obs", ylab="SSE, sim")
abline(1,1,col="blue",lwd=2)
text(20,200,sprintf("Bayesian p: %#.3g",pval))
par(mai=c(0.9,0.9,0,0))
boxplot(gof2$value~gof2$count,xlab="Tot captures",ylab="# surveyed sites")
text((0:6)+1,gof2_obs,"X",cex=1.25,col="blue")
dev.off()

##############
# Visualize posterior densities as histograms
##############

hist(mod$BUGSoutput$sims.list$basecount)
mean(mod$BUGSoutput$sims.list$basecount)  # count - 0.004 (per trap night...)

hist(mod$BUGSoutput$sims.list$r,breaks=20)

hist(mod$BUGSoutput$sims.list$baseocc)
mean(mod$BUGSoutput$sims.list$baseocc)     
median(mod$BUGSoutput$sims.list$baseocc)    

#hist(mod$BUGSoutput$sims.list$beta.elev.o)            # not much effect of elevation?

hist(mod$BUGSoutput$sims.list$beta.ecoregion.o[,2])     # effect of ecoregion: slightly lower occupancy for mojave
hist(mod$BUGSoutput$sims.list$beta.ecoregion.o[,3])     # very slightly lower occupancy for transition?

#hist(mod$BUGSoutput$sims.list$beta.landform.o[,2])     # not much effect of landform. valley slightly more likely to be occupied?

hist(mod$BUGSoutput$sims.list$beta.moisture.o[,2])      # strong effect of moisture

# hist(mod$BUGSoutput$sims.list$beta.slope.o[,2])      # lower occupancy for high slopes
# hist(mod$BUGSoutput$sims.list$beta.slope.o[,3])     # higher occupancy for low slopes

hist(mod$BUGSoutput$sims.list$beta.aspect.o[,2])      # strong effect of south aspect vs. varied
hist(mod$BUGSoutput$sims.list$beta.aspect.o[,3])     # slightly negative effect of north aspect

hist(mod$BUGSoutput$sims.list$beta.substrate.o[,2])    # no effect of rock
hist(mod$BUGSoutput$sims.list$beta.substrate.o[,3])   # no effect of gravel

# hist(mod$BUGSoutput$sims.list$beta.vegassoc.o[,2])  # positive effect of sagebrush
# hist(mod$BUGSoutput$sims.list$beta.vegassoc.o[,3])  # slight positive effect of blackbrush?
# hist(mod$BUGSoutput$sims.list$beta.vegassoc.o[,4])  # no effect of ephedra
# hist(mod$BUGSoutput$sims.list$beta.vegassoc.o[,5])  # no effect of bursage

hist(mod$BUGSoutput$sims.list$beta.elev.o[])   # no effect of elev
hist(mod$BUGSoutput$sims.list$beta.elev2.o[])  # highly significant curvature on elevation relationship


##########
# make histogram type figure

png("coefficients.png",width=900,height = 700)
layout(matrix(c(1,2,3,4,5,6,7,8),nrow=4,byrow=T))

hist(mod$BUGSoutput$sims.list$beta.elev.o[],main="Elevation",xlab="",freq=F)   # no effect of elev
lines(density(mod$BUGSoutput$sims.list$beta.elev.o[]),col="blue",lty=2,lwd=2)
temp <- HPDinterval(runjags::combine.mcmc(mod2)[,"beta.elev.o"])
arrows(temp["var1","lower"],0.05,temp["var1","upper"],0.05,col="red",lwd=3,code=3,angle=90,length=0.1)

hist(mod$BUGSoutput$sims.list$beta.elev2.o[],main="Elevation^2",xlab="",freq=F)  # highly significant curvature on elevation relationship
lines(density(mod$BUGSoutput$sims.list$beta.elev2.o[]),col="blue",lty=2,lwd=2)
temp <- HPDinterval(runjags::combine.mcmc(mod2)[,"beta.elev2.o"])
arrows(temp["var1","lower"],0.05,temp["var1","upper"],0.05,col="red",lwd=3,code=3,angle=90,length=0.1)

hist(mod$BUGSoutput$sims.list$beta.moisture.o[,2],main="Wet/Mesic",xlab="",freq=F)      # strong effect of moisture
lines(density(mod$BUGSoutput$sims.list$beta.moisture.o[,2]),col="blue",lty=2,lwd=2)
temp <- HPDinterval(runjags::combine.mcmc(mod2)[,"beta.moisture.o[2]"])
arrows(temp["var1","lower"],0.05,temp["var1","upper"],0.05,col="red",lwd=3,code=3,angle=90,length=0.1)

hist(mod$BUGSoutput$sims.list$beta.aspect.o[,2],main="South-facing slope",xlab="",freq=F)
lines(density(mod$BUGSoutput$sims.list$beta.aspect.o[,2]),col="blue",lty=2,lwd=2)
temp <- HPDinterval(runjags::combine.mcmc(mod2)[,"beta.aspect.o[2]"])
arrows(temp["var1","lower"],0.05,temp["var1","upper"],0.05,col="red",lwd=3,code=3,angle=90,length=0.1)

hist(mod$BUGSoutput$sims.list$beta.substrate.o[,2],main="Rocky substrate",xlab="Posterior dist for coefficient estimate",freq=F)
lines(density(mod$BUGSoutput$sims.list$beta.substrate.o[,2]),col="blue",lty=2,lwd=2)
temp <- HPDinterval(runjags::combine.mcmc(mod2)[,"beta.substrate.o[2]"])
arrows(temp["var1","lower"],0.05,temp["var1","upper"],0.05,col="red",lwd=3,code=3,angle=90,length=0.1)

hist(mod$BUGSoutput$sims.list$beta.ecoregion.o[,2],main="Mojave ecoregion",xlab="Posterior dist for coefficient estimate",freq=F)
lines(density(mod$BUGSoutput$sims.list$beta.ecoregion.o[,2]),col="blue",lty=2,lwd=2)
temp <- HPDinterval(runjags::combine.mcmc(mod2)[,"beta.ecoregion.o[2]"])
arrows(temp["var1","lower"],0.05,temp["var1","upper"],0.05,col="red",lwd=3,code=3,angle=90,length=0.1)

hist(mod$BUGSoutput$sims.list$beta.vegassoc.o[,2],main="Sagebrush vegetation",xlab="Posterior dist for coefficient estimate",freq=F)
lines(density(mod$BUGSoutput$sims.list$beta.vegassoc.o[,2]),col="blue",lty=2,lwd=2)
temp <- HPDinterval(runjags::combine.mcmc(mod2)[,"beta.vegassoc.o[2]"])
arrows(temp["var1","lower"],0.05,temp["var1","upper"],0.05,col="red",lwd=3,code=3,angle=90,length=0.1)

hist(mod$BUGSoutput$sims.list$beta.vegassoc.o[,3],main="Blackbrush vegetation",xlab="Posterior dist for coefficient estimate",freq=F)
lines(density(mod$BUGSoutput$sims.list$beta.vegassoc.o[,3]),col="blue",lty=2,lwd=2)
temp <- HPDinterval(runjags::combine.mcmc(mod2)[,"beta.vegassoc.o[3]"])
arrows(temp["var1","lower"],0.05,temp["var1","upper"],0.05,col="red",lwd=3,code=3,angle=90,length=0.1)

dev.off()

##########
# make dot type figure

svg("coefficients2.svg",width=5,height = 5)

par(mai=c(0.9,2,0,0))
plot(1,1,pty="n",main="",ylim=c(0,13),xlim=c(-4,4),xlab="Standardized regression coefficient",ylab="",pch="",bty="n",yaxt="n")
abline(v=0,lwd=3)
abline(h=seq(0.5,12.5,1),lty=2,lwd=0.5)
axlabs <- c()

#hist(mod$BUGSoutput$sims.list$beta.elev.o[],main="Elevation",xlab="",freq=F)   # no effect of elev
#lines(density(mod$BUGSoutput$sims.list$beta.elev.o[]),col="blue",lty=2,lwd=2)
temp <- HPDinterval(runjags::combine.mcmc(mod2)[,"beta.elev.o"],prob=0.9)
arrows(temp["var1","lower"],13,temp["var1","upper"],13,col="red",lwd=3,code=3,angle=90,length=0.1)
axlabs=c(axlabs,"Elevation")

#hist(mod$BUGSoutput$sims.list$beta.elev2.o[],main="Elevation^2",xlab="",freq=F)  # highly significant curvature on elevation relationship
#lines(density(mod$BUGSoutput$sims.list$beta.elev2.o[]),col="blue",lty=2,lwd=2)
temp <- HPDinterval(runjags::combine.mcmc(mod2)[,"beta.elev2.o"],prob=0.9)
arrows(temp["var1","lower"],12,temp["var1","upper"],12,col="red",lwd=3,code=3,angle=90,length=0.1)
axlabs=c(axlabs,"Elevation^2")

#hist(mod$BUGSoutput$sims.list$beta.moisture.o[,2],main="Wet/Mesic",xlab="",freq=F)      # strong effect of moisture
#lines(density(mod$BUGSoutput$sims.list$beta.moisture.o[,2]),col="blue",lty=2,lwd=2)
temp <- HPDinterval(runjags::combine.mcmc(mod2)[,"beta.moisture.o[2]"],prob=0.9)
arrows(temp["var1","lower"],11,temp["var1","upper"],11,col="red",lwd=3,code=3,angle=90,length=0.1)
axlabs=c(axlabs,"Wet/Mesic")

#hist(mod$BUGSoutput$sims.list$beta.aspect.o[,2],main="South-facing slope",xlab="",freq=F)
#lines(density(mod$BUGSoutput$sims.list$beta.aspect.o[,2]),col="blue",lty=2,lwd=2)
temp <- HPDinterval(runjags::combine.mcmc(mod2)[,"beta.aspect.o[2]"],prob=0.9)
arrows(temp["var1","lower"],10,temp["var1","upper"],10,col="red",lwd=3,code=3,angle=90,length=0.1)
axlabs=c(axlabs,"South-facing slope")

#hist(mod$BUGSoutput$sims.list$beta.aspect.o[,3],main="South-facing slope",xlab="",freq=F)
#lines(density(mod$BUGSoutput$sims.list$beta.aspect.o[,3]),col="blue",lty=2,lwd=2)
temp <- HPDinterval(runjags::combine.mcmc(mod2)[,"beta.aspect.o[3]"],prob=0.9)
arrows(temp["var1","lower"],9,temp["var1","upper"],9,col="red",lwd=3,code=3,angle=90,length=0.1)
axlabs=c(axlabs,"North-facing slope")

#hist(mod$BUGSoutput$sims.list$beta.substrate.o[,2],main="Soil substrate",xlab="Posterior dist for coefficient estimate",freq=F)
#lines(density(mod$BUGSoutput$sims.list$beta.substrate.o[,2]),col="blue",lty=2,lwd=2)
temp <- HPDinterval(runjags::combine.mcmc(mod2)[,"beta.substrate.o[2]"],prob=0.9)
arrows(temp["var1","lower"],8,temp["var1","upper"],8,col="red",lwd=3,code=3,angle=90,length=0.1)
axlabs=c(axlabs,"Soil substrate")

#hist(mod$BUGSoutput$sims.list$beta.substrate.o[,2],main="Soil substrate",xlab="Posterior dist for coefficient estimate",freq=F)
#lines(density(mod$BUGSoutput$sims.list$beta.substrate.o[,2]),col="blue",lty=2,lwd=2)
temp <- HPDinterval(runjags::combine.mcmc(mod2)[,"beta.substrate.o[3]"],prob=0.9)
arrows(temp["var1","lower"],7,temp["var1","upper"],7,col="red",lwd=3,code=3,angle=90,length=0.1)
axlabs=c(axlabs,"Rock substrate")

#hist(mod$BUGSoutput$sims.list$beta.substrate.o[,2],main="Soil substrate",xlab="Posterior dist for coefficient estimate",freq=F)
#lines(density(mod$BUGSoutput$sims.list$beta.substrate.o[,2]),col="blue",lty=2,lwd=2)
temp <- HPDinterval(runjags::combine.mcmc(mod2)[,"beta.substrate.o[4]"],prob=0.9)
arrows(temp["var1","lower"],6,temp["var1","upper"],6,col="red",lwd=3,code=3,angle=90,length=0.1)
axlabs=c(axlabs,"Boulder substrate")

#hist(mod$BUGSoutput$sims.list$beta.ecoregion.o[,2],main="Mojave ecoregion",xlab="Posterior dist for coefficient estimate",freq=F)
#lines(density(mod$BUGSoutput$sims.list$beta.ecoregion.o[,2]),col="blue",lty=2,lwd=2)
temp <- HPDinterval(runjags::combine.mcmc(mod2)[,"beta.ecoregion.o[2]"],prob=0.9)
arrows(temp["var1","lower"],5,temp["var1","upper"],5,col="red",lwd=3,code=3,angle=90,length=0.1)
axlabs=c(axlabs,"Mojave ecoregion")

#hist(mod$BUGSoutput$sims.list$beta.ecoregion.o[,2],main="Mojave ecoregion",xlab="Posterior dist for coefficient estimate",freq=F)
#lines(density(mod$BUGSoutput$sims.list$beta.ecoregion.o[,2]),col="blue",lty=2,lwd=2)
temp <- HPDinterval(runjags::combine.mcmc(mod2)[,"beta.ecoregion.o[3]"],prob=0.9)
arrows(temp["var1","lower"],4,temp["var1","upper"],4,col="red",lwd=3,code=3,angle=90,length=0.1)
axlabs=c(axlabs,"Transition ecoregion")

#hist(mod$BUGSoutput$sims.list$beta.vegassoc.o[,2],main="Sagebrush vegetation",xlab="Posterior dist for coefficient estimate",freq=F)
#lines(density(mod$BUGSoutput$sims.list$beta.vegassoc.o[,2]),col="blue",lty=2,lwd=2)
temp <- HPDinterval(runjags::combine.mcmc(mod2)[,"beta.vegassoc.o[2]"],prob=0.9)
arrows(temp["var1","lower"],3,temp["var1","upper"],3,col="red",lwd=3,code=3,angle=90,length=0.1)
axlabs=c(axlabs,"Sagebrush vegetation")

#hist(mod$BUGSoutput$sims.list$beta.vegassoc.o[,3],main="Blackbrush vegetation",xlab="Posterior dist for coefficient estimate",freq=F)
#lines(density(mod$BUGSoutput$sims.list$beta.vegassoc.o[,3]),col="blue",lty=2,lwd=2)
temp <- HPDinterval(runjags::combine.mcmc(mod2)[,"beta.vegassoc.o[3]"],prob=0.9)
arrows(temp["var1","lower"],2,temp["var1","upper"],2,col="red",lwd=3,code=3,angle=90,length=0.1)
axlabs=c(axlabs,"Blackbrush vegetation")

#hist(mod$BUGSoutput$sims.list$beta.vegassoc.o[,3],main="Blackbrush vegetation",xlab="Posterior dist for coefficient estimate",freq=F)
#lines(density(mod$BUGSoutput$sims.list$beta.vegassoc.o[,3]),col="blue",lty=2,lwd=2)
temp <- HPDinterval(runjags::combine.mcmc(mod2)[,"beta.vegassoc.o[4]"],prob=0.9)
arrows(temp["var1","lower"],1,temp["var1","upper"],1,col="red",lwd=3,code=3,angle=90,length=0.1)
axlabs=c(axlabs,"Ephedra vegetation")

#hist(mod$BUGSoutput$sims.list$beta.vegassoc.o[,3],main="Blackbrush vegetation",xlab="Posterior dist for coefficient estimate",freq=F)
#lines(density(mod$BUGSoutput$sims.list$beta.vegassoc.o[,3]),col="blue",lty=2,lwd=2)
temp <- HPDinterval(runjags::combine.mcmc(mod2)[,"beta.vegassoc.o[5]"],prob=0.9)
arrows(temp["var1","lower"],0,temp["var1","upper"],0,col="red",lwd=3,code=3,angle=90,length=0.1)
axlabs=c(axlabs,"Bursage vegetation")

axis(2,at=0:13,labels=rev(axlabs),las=2,col=NA)

dev.off()

#########
# Visualize covariate effects on occupancy
#########


#########
# Plot site-level occupancy

hist(mod$BUGSoutput$sims.list$occupied[,1])   # this site is probably not occupied...
occups <- sapply(1:nrow(sitecovs),function(t) (table(mod$BUGSoutput$sims.list$occupied[,t])/length(mod$BUGSoutput$sims.list$occupied[,t]))["1"] )
occups[which(is.na(occups))] <- 0
names(occups) <- sitecovs$LocationCode[ndx]

ns <- round(length(occups)/2)

svg("allsites.svg",width=12,height = 5)
layout(matrix(c(1,2),nrow=2))
par(mai=c(1,0.8,0,0))
barplot(sort(occups)[1:ns],names.arg = names(sort(occups))[1:ns],las=2,ylim=c(0,1),xlab="location code",ylab="Predicted occupancy")
barplot(sort(occups)[(ns+1):length(occups)],names.arg = names(sort(occups))[(ns+1):length(occups)],las=2,ylim=c(0,1),xlab="location code",ylab="Predicted occupancy")
dev.off()


########
# total occupancy

quantile(apply(mod$BUGSoutput$sims.list$occupied,1,sum)/nsites,c(0.05,0.5,0.95))


#######
# Circular barplot?

# Libraries
library(tidyverse)

# Create dataset
data=data.frame(
  id=as.numeric(1:nsites),
  individual=names(sort(occups)),
  value=sort(occups)
)

# ----- This section prepare a dataframe for labels ---- #
# Get the name and the y position of each label
label_data=data

# calculate the ANGLE of the labels
number_of_bar=nrow(label_data)
angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)

# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
label_data$hjust<-ifelse( angle < -90, 1, 0)

# flip angle BY to make them readable
label_data$angle<-ifelse(angle < -90, angle+180, angle)
# ----- ------------------------------------------- ---- #

#png("allsites2.png",width=1500,height = 1200)
svg("allsites2.svg",8,8)
# Start the plot
p = ggplot(data, aes(x=as.factor(id), y=value)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a blue color
  geom_bar(stat="identity", fill=alpha("skyblue", 0.7)) +
  
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-1,1.2) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm")      # Adjust the margin to make in sort labels are not truncated!
  ) +
  
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = 0) +
  
  # Add the labels, using the label_data dataframe that we have created before
  geom_text(data=label_data, aes(x=id, y=value+.1, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) 

p
dev.off()

#########
# Plot effects of covars on occupancy

namestable$name1

predict.jags("elev")
predict.jags("ecoregion")
predict.jags("aspect")
predict.jags("substrate")
#predict.jags("landform")
predict.jags("vegassoc")
predict.jags("moisture")
#predict.jags("slope")
#predict.jags("aspect")




##################################
######################
# ALTERNATIVE METHOD: UNMARKED (more standard occupancy analysis)

#library(help="unmarked")

#vignette("unmarked",package="unmarked")

ndx <- match(as.numeric(rownames(caphist)),sitecovs$LocationCode)

ndx2 <- !duplicated(rownames(caphist))   # remove duplicates?


umf1 <- unmarkedFrameOccu(y = caphist[ndx2,],
              siteCovs = sitecovs[ndx[ndx2],])

umf2 <- unmarkedFrameOccu(y = caphist,
                          siteCovs = sitecovs[ndx,])

###########
# Null model

mod1 <- occu(~1 ~1, umf2)
summary(mod1)
backTransform(mod1, 'state')       # occupancy of 0.219
backTransform(mod1, 'det')       # detection of 0.203

###########
# Elevation only

mod2 <- occu(~1 ~scale(Elevation), umf2)
summary(mod2)

mod3 <- occu(~1 ~Slope, umf2)
summary(mod3)

mod4 <- occu(~1 ~Aspect, umf2)
summary(mod4)

mod5 <- occu(~1 ~Moisture, umf2)
summary(mod5)

mod6 <- occu(~1 ~Slope+SoilSubstrate+Ecoregion, umf2)
summary(mod6)

mod7 <- occu(~1 ~Slope+Aspect+Ecoregion, umf2)
summary(mod6)
backTransform(linearComb(mod6, coefficients = c(1,0,0,0,0,0), type = 'state'))
backTransform(mod1, 'det') 

mod8 <- occu(~1 ~Slope*Aspect+Ecoregion, umf2)
summary(mod8)
backTransform(linearComb(mod6, coefficients = c(1,0,0,0,0,0), type = 'state'))
backTransform(mod1, 'det')


##########
# Model selection

mods <- fitList('psi(.)p(.)' = mod1,
                'psi(Elevation)p(.)' = mod2,
                'psi(Slope)p(.)' = mod3,
                'psi(Aspect)p(.)' = mod4,
                'psi(Moisture)p(.)' = mod5,
                'psi(Slope+Substrate+Ecoregion)p(.)' = mod6,
                'psi(Slope+Ecoregion+Aspect)p(.)' = mod7,
                'psi(Slope*Aspect+Ecoregion)p(.)' = mod8
                )
modSel(mods)


#########
# raw estimates

caphist2 <-caphist[Occupied==1,] 

firsts <- sapply(1:nrow(caphist2),function(t) which(caphist2[t,]==1)[1] )
nbouts <- ncol(caphist2)

toadd <- numeric(0)
for(i in 1:nrow(caphist2)){
  if((firsts[i]+1)<nbouts)  toadd <- c(toadd,caphist2[i,(firsts[i]+1):nbouts])
}

toadd <- toadd[!is.na(toadd)]

detection <- sum(toadd)/length(toadd)   #detection should be around 0.21  (like "unmarked" says)



################
# Alternative: Zero Inflated Poisson
################


# NOTe: use DHARMa to generate residuals and then perform likelihood-ratio goodness-of-fit tests


#############
# Diagnose which count distribution to use (remember this is not looking at residuals, just giving us a first pass to see which distribution is best) 

# test for poisson
fit <- vcd::goodfit(captures_total, type="poisson")           # goodness of fit test for poisson
summary(fit)                                                  # clearly NOT poisson
vcd::rootogram(fit)                                           # terrible fit
vcd::Ord_plot(captures_total)                                 # this diagnoses that negative binomial distribution is appropriate

# test for negative binomial
fit <- vcd::goodfit(captures_total, type="nbinom",method="ML")  # goodness of fit test for neg binom
summary(fit)                                                    # not a good fit- but better!  This is a strict test
vcd::rootogram(fit)                                             # not great fit
vcd::distplot(captures_total, type="nbinom")                    # Neg binom is okay

# visualize poisson fit
lambda <- MASS::fitdistr(captures_total,"poisson")$estimate
hist(captures_total,freq=F,ylim=c(0,1),breaks=seq(0,4,1))   # very skewed   #breaks=seq(0,1401,10)
curve(dpois(x,lambda),add=T,col="green",lwd=2,from=0,to=4)


# visualize nbinom fit
size <- MASS::fitdistr(captures_total,"negative binomial")$estimate["size"]
mu <- MASS::fitdistr(captures_total,"negative binomial")$estimate["mu"]
hist(captures_total,freq=F,ylim=c(0,1),breaks=seq(0,4,1))   # very skewed   #breaks=seq(0,1401,10)
curve(dnbinom(x,size,mu=mu),add=T,col="green",lwd=2,from=0,to=4)   # much better fit


############### 
# Preliminary linear model- poisson (BAD MODEL)

df <- sitecovs
df$totcaps <- 0
df$totcaps[ndx] <- captures_total

df$effort <- 0
df$effort[ndx] <- effort/180 

df$Elevation_sc <- scale(df$Elevation)

names(df)

mod.pois <- glm(totcaps ~ Elevation+Moisture+Aspect,data=df,maxit=100,family="poisson")  
summary(mod.pois)

boot::glm.diag.plots(mod.pois)             # not great
anova(mod.pois, test="Chisq")              # test vs null model. lots of information gain with moisture and aspect
deviance(mod.pois)/mod.pois$df.residual    # overdispersion okay?
dispersiontest(mod.pois)                   # fails the test... overdispersed!
car::influencePlot(mod.pois)


############### 
# Preliminary linear model- negative binomial


mod.nb <- MASS::glm.nb(totcaps ~ Elevation+Moisture+Aspect,data=df)     #,maxit=1000
summary(mod.nb)

boot::glm.diag.plots(mod.nb)             # some high influence points, otherwise okay
anova(mod.nb, test="Chisq")              # test vs null model. Seems to be some information in the "Slope" variable
deviance(mod.nb)/mod.nb$df.residual      # some overdispersion, but not too bad
car::influencePlot(mod.nb)                    # some highly influential points

resDHARMa = DHARMa::simulateResiduals(mod.nb)   # better way to examine residuals for glms
plot(resDHARMa)                          # much better! Residuals look okay
DHARMa::testUniformity(resDHARMa)        # passes the test!


par(mfrow = c(1,2))                      # examine residuals across the range of parameter space. Examine for possible non-linearities
DHARMa::plotResiduals(df$Elevation,  resDHARMa$scaledResiduals)
DHARMa::plotResiduals(df$Moisture,  resDHARMa$scaledResiduals)  # relationship driven by one point



##################################
## if the mod.nb doesn't work... ??


## option 1: use maximum likelihood model fitting directly?

# library(bbmle)
# library(emdbook)
# m1 <- mle2(Total.Elk~dnbinom(mu=exp(logmu),size=1/invk),
#            param=list(logmu~Wetlands+Slope+Grassland),
#            start=list(logmu=0,invk=1),
#            data=df)
# prof <- profile(m1)
# confint(prof)   # no significance?
# confint(m1,method="uniroot")   # alternative method for confidence intervals...


############
## option 2: use zero inflated poisson model?

mod.zip <- zeroinfl(totcaps ~ 1 | offset(effort) + Elevation_sc + I(Elevation_sc^2) + Moisture, data=df, dist = "poisson")     # zero inflated negative binomial
# offset(effort) + Aspect + Moisture  + Ecoregion


### model selection: not really needed- if we can get it in the Bayesian model, it can be included (just present the confidence interval.)
###  no interactions are possible

summary(mod.zip)

mod.zip$loglik

 # predprob(mod.zip)      #extract probabilities






VIF <- (-2*a$loglik)/mod.zinb$df.residual      # variance inflation factor 
VIF  # better...


AIC(mod.zinb, mod.nb)   # zero inflation model has much better AIC

vuong(mod.nb,mod.zinb)   # zero inflated model is superior





VIF
?zeroinfl



##### Play with bbmle 

load(system.file("vignetteData","orob1.rda",package="bbmle"))
summary(orob1)

m1B <- mle2(m~dbetabinom(prob,size=n,theta),
            param=list(prob~dilution),
            start=list(prob=0.5,theta=1),
            data=orob1)







################
# OLD CODE
################

# site_session_linker <- data.frame(
#   sitenames1 = rownames(caphist),
#   index = match(rownames(caphist),sitecovs$LocationCode),
#   stringsAsFactors = F
# )

#################
# notes


# With the assumption that occupied sites are always occupied, there were some strange results!




