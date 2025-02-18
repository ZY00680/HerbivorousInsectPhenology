######################
# Appendix S3 
# of
# Towards more predictive and generalizable herbivorous insect phenology models
# by 
# Zimo Yang, Elise Woodruff, David Held, and Nate B Hardy
# Artifical in Ecological Application
#######################
# Contents

# 1. General Data pre-process
# 2. LT taxonomy models
# 3. DD taxonomy models
# 4. DD approximation percent error
# 5. LT variation among development stages
# 6. Variation in LT variation among development stage; Not explicitly dicussed in main text due to poor performance
# 7. R2 calculation for taxonomy models
# 8. Phylogeny models and R2 calculation
# 9. Chi-square tests to describe database coverage
######################



######################
library(emmeans)
library(performance)
library(ggplot2)
library(sjPlot)
library(gridExtra)
library(dplyr)
library(brms)
library(geiger)
library(nlme)
library(bayestestR)

#######################
# Pre-process
######################

d <- read.csv('/home/zimo/Desktop/Zimo-PhenModels/phenoDBherbiesDataPack/phenoDBherbies.4.5.csv')

# Natural log transformation, refactor reference level
d <- d[-1,] # First row is column notes.
d$th_family <- as.factor(d$th_family)
d$Growth_Form <- as.factor(d$Growth_Form)
d$tMode <- as.factor(d$tMode)
d$Herb <- as.factor(d$Herb)
d$Wood <- as.factor(d$Wood)
d$th_family <- relevel(d$th_family, ref = "Fabaceae")
d$eggAdDD <- log(as.numeric(d$eggAdDD))
d$Host_Breadth <- log(as.numeric(d$Host_Breadth))
d$Body_Size <- log(as.numeric(d$Body_Size))
d$base_temp_eggAd <- as.numeric(d$base_temp_eggAd)
d$base_temp_Egg <- as.numeric(d$base_temp_Egg)
d$base_temp_Larvae <- as.numeric(d$base_temp_Larvae)
d$base_temp_Pupae <- as.numeric(d$base_temp_Pupae)
d$eggDD <- log(as.numeric(d$eggDD))
d$laDD <- log(as.numeric(d$laDD))
d$puDD <- log(as.numeric(d$puDD))
d$Latitude <- as.numeric(d$Latitude)
d$Longitude <- as.numeric(d$Longitude)
d$Experiment <- as.factor(d$Experiment)
d$lat.range <- as.numeric(d$lat.range)
d$lon.range <- as.numeric(d$lon.range)
d$X <- as.factor(d$X) # Experiment marker

# Add life history identifier
d$LH[d$Order %in% c("Hemiptera", "Phasmatodea", "Psocoptera", "Orthoptera")] <- "Hemi"
d$LH[d$Order %in% c("Coleoptera","Diptera","Hymenoptera","Lepidoptera","Thysanoptera")] <- "Holo"

# Filter out records with inaccurate location
d <- d %>% filter(lat.range <= 10) %>% # Latitude range threshold = 10 degree.
	filter(Estimated_Location != "E")# %>% # Exclude location guessed from researcher's affiliation.

# Scaling of continuous effects
d$absLatitude = scale(abs(d$Latitude))
d$Host_Breadth = scale(d$Host_Breadth)
d$Body_Size = scale(d$Body_Size)

# A second dataframe without artificial diet to allow hierarchical design of tested host plant taxa.
d2 <- subset(d, d$th_family!="AF")

stop('Data loaded')
#######################
# LT taxonomy model
#######################
# Egg LT & larvae/nymph LT. Since they are different, no need for overall immature LT
LTE1 <- brm(base_temp_Egg ~ Body_Size + Host_Breadth + Experiment + tMode + absLatitude + (1|Order:Family) +(1|th_family), family=gaussian(), control = list(adapt_delta = 0.95),thin=1, iter=4000, save_pars = save_pars(all = TRUE),data=d[d$Length_measured == "Body",])
LTE2 <- brm(base_temp_Egg ~ Body_Size + Host_Breadth + Experiment + tMode + absLatitude + (1|Order:Family) +(1|th_family:th_genus), family=gaussian(), control = list(adapt_delta = 0.95),thin=1, iter=4000, save_pars = save_pars(all = TRUE),data=d2[d2$Length_measured == "Body",])
LTL1 <- brm(base_temp_Larvae ~ Body_Size + Host_Breadth + Experiment + tMode + absLatitude + (1|Order:Family) +(1|th_family), family=gaussian(), control = list(adapt_delta = 0.95),thin=1, iter=4000, save_pars = save_pars(all = TRUE),data=d[d$Length_measured == "Body",])
LTL2 <- brm(base_temp_Larvae ~ Body_Size + Host_Breadth + Experiment + tMode + absLatitude + (1|Order:Family) +(1|th_family:th_genus), family=gaussian(), control = list(adapt_delta = 0.95),thin=1, iter=4000, save_pars = save_pars(all = TRUE),data=d2[d2$Length_measured == "Body",])

# Post hoc test on tMode
LTE1_ph <- emmeans(LTE1, pairwise~tMode)
LTE2_ph <- emmeans(LTE2, pairwise~tMode)

LTL1_ph <- emmeans(LTL1, pairwise~tMode)
LTL2_ph <- emmeans(LTL2, pairwise~tMode)

DD1_ph <- emmeans(DD1, pairwise~tMode)
DD2_ph <- emmeans(DD2, pairwise~tMode)


#######################
# DD taxonomy model
#######################
DD1 <- brm(eggAdDD ~ Body_Size + Host_Breadth  + Experiment + tMode + absLatitude + (1|Order:Family) +(1|th_family), family=gaussian(), control = list(adapt_delta = 0.95),thin=1, iter=4000, save_pars = save_pars(all = TRUE),data=d[d$Length_measured == "Body",])
DD2 <- brm(eggAdDD ~ Body_Size + Host_Breadth  + Experiment + tMode + absLatitude + (1|Order:Family) +(1|th_family:th_genus), family=gaussian(), control = list(adapt_delta = 0.95),thin=1, iter=4000, save_pars = save_pars(all = TRUE),data=d2[d2$Length_measured == "Body",])

########################
# percent error of DD approximation 
########################
# Estimate error introduced from approximating regression of egg to adult with sum of regressions of each developmental stage
d3 <- read.csv('/home/zimo/Desktop/Zimo-PhenModels/phenoDBherbiesDataPack/phenoDBherbies.4.5.csv')
d3$eggDD <- as.numeric(d3$eggDD)
d3$eggAdDD <- as.numeric(d3$eggAdDD)
d3$laDD <- as.numeric(d3$laDD)
d3$puDD <- as.numeric(d3$puDD)
app <- c(which(d3$eggAdDD == d3$eggDD+d3$laDD+d3$puDD),which(d3$eggAdDD == d3$eggDD+d3$laDD))
d3 <- d3[-app,]

d3$LH[d3$Order %in% c("Hemiptera", "Phasmatodea", "Psocoptera", "Orthoptera")] <- "Hemi"
d3$LH[d3$Order %in% c("Coleoptera","Diptera","Hymenoptera","Lepidoptera","Thysanoptera")] <- "Holo"

d31 <- d3[d3$LH=='Holo',]
d31 <- d31[complete.cases(d31[, c('eggAdDD','eggDD','laDD','puDD')]),]
d31$A <- d31$eggDD + d31$laDD + d31$puDD
d31$B <- (d31$A - d31$eggAdDD)/d31$eggAdDD

d32 <- d3[d3$LH=='Hemi',]
d32 <- d32[complete.cases(d32[, c('eggAdDD','eggDD','laDD')]),]
d32$A <- d32$eggDD + d32$laDD
d32$B <- (d32$A - d32$eggAdDD)/d32$eggAdDD

appError <- c(d31$B,d32$B)
mean(appError)
sd(appError)
shapiro.test(appError)
hist(appError)

########################
# LT difference among developmental stages
########################
# Data pre-processing
# Pile LT and DD each into one column, add a column specifying stage

# Get records with full information for model fitting
LTcol <- c("base_temp_Egg","base_temp_Larvae",# No life history yet
	   "Body_Size","Host_Breadth","Experiment","tMode","absLatitude",
	   "Order","Family","th_family","th_genus")

LTc <- d[complete.cases(d[, LTcol]),]
LTc <- LTc[LTc$Length_measured == "Body",]
LTcHe <- LTc[LTc$LH == "Hemi",]
LTcHo <- LTc %>% filter(LH=="Holo") %>% filter(!(is.na(base_temp_Pupae))) # Holo

# Reframe data
sizeLTHe<-length(LTcHe$Order)
rLTcHe <- LTcHe %>% 
reframe(base_temp=c(base_temp_Egg,base_temp_Larvae),
      Body_Size = rep(Body_Size,2),
      Host_Breadth = rep(Host_Breadth,2),
      Experiment = rep(Experiment,2),
      tMode = rep(tMode,2),
      absLatitude = rep(Latitude,2),
      Order = rep(Order,2),
      Family = rep(Family,2),
      th_family = rep(th_family,2),
      th_genus = rep(th_genus,2),
			stage=c(rep("egg",sizeLTHe),rep("la",sizeLTHe)),
			X=rep(X,2))
			
sizeLTHo<-length(LTcHo$Order)			
rLTcHo <- LTcHo %>% 
reframe(base_temp=c(base_temp_Egg,base_temp_Larvae,base_temp_Pupae),
        Body_Size = rep(Body_Size,3),
        Host_Breadth = rep(Host_Breadth,3),
        Experiment = rep(Experiment,3),
        tMode = rep(tMode,3),
        absLatitude = rep(Latitude,3),
        Order = rep(Order,3),
        Family = rep(Family,3),
        th_family = rep(th_family,3),
        th_genus = rep(th_genus,3),
        stage=c(rep("egg",sizeLTHo),rep("la",sizeLTHo),rep("pu",sizeLTHo)),
		  	X=rep(X,3))
			
rLTc <- rbind(rLTcHo,rLTcHe)

rLTc$stage <- as.factor(rLTc$stage)
rLTc$stage <- relevel(rLTc$stage,'egg')


# Analysis the effect of stage under the framework of the full models

# This one cannot stably converge.
#LTstage1<- brm(base_temp ~ stage + Body_Size + Host_Breadth  + Experiment + tMode + absLatitude + (1|Order:Family) +(1|th_family) + (1|X), family=gaussian(), control = list(adapt_delta = 0.95),thin=1, iter=4000, save_pars = save_pars(all = TRUE),
#                data=rLTc) 

LTstage2 <- brm(base_temp ~ stage + Body_Size + Host_Breadth  + Experiment + tMode + absLatitude + (1|Order:Family) +(1|th_family:th_genus) + (1|X), family=gaussian(), control = list(adapt_delta = 0.95),thin=1, iter=4000, save_pars = save_pars(all = TRUE),
                data=rLTc[rLTc$th_family!='AF',]) 

# Post hoc for stage
LTstage2_ph <- emmeans(LTstage2,pairwise~stage)
LTstage2_ph

########################
# Variation in LT variation among stages
# Poor performance
########################
# Take coefficient of variance across stage for each experiment
#LTcHo <- LTcHo %>% 
#	rowwise() %>% 
#	mutate(COV = sd(c(base_temp_Egg,base_temp_Larvae,base_temp_Pupae))/mean(c(base_temp_Egg,base_temp_Larvae,base_temp_Pupae)),
#	       LTavg = mean(c(base_temp_Egg,base_temp_Larvae,base_temp_Pupae)),
#	       LTsd = sd(c(base_temp_Egg,base_temp_Larvae,base_temp_Pupae)))

#LTcHe <- LTcHe %>% 
#	rowwise() %>% 
#	mutate(COV = sd(c(base_temp_Egg,base_temp_Larvae))/mean(c(base_temp_Egg,base_temp_Larvae)),
#	       LTavg = mean(c(base_temp_Egg,base_temp_Larvae)),
#	       LTsd = sd(c(base_temp_Egg,base_temp_Larvae)))

#LTcHo <- data.frame(LTcHo) 
#LTcHe <- data.frame(LTcHe)

#LTcV <- rbind(LTcHo,LTcHe)


#LTv1 <- lmer(COV ~ Body_Size + Host_Breadth + th_form + Experiment + tMode + abs(Latitude) + (1|Order:Family) + (1|th_family)-1, data=LTcV)
#LTv2 <- lmer(COV ~ Body_Size + Host_Breadth + th_form + Experiment + tMode + abs(Latitude) + (1|Order:Family) + (1|th_family:th_genus)-1, data=LTcV[LTcV$th_form != 'AF',]) 


#LTv1 <- brm(COV ~ Body_Size + Host_Breadth  + Experiment + tMode + abs(Latitude) + (1|Order:Family) +(1|th_family), family=gaussian(), control = list(adapt_delta = 0.95),thin=1, iter=4000, save_pars = save_pars(all = TRUE),data=LTcV)
#LTv2 <- brm(COV ~ Body_Size + Host_Breadth  + Experiment + tMode + abs(Latitude) + (1|Order:Family) +(1|th_family:th_genus), family=gaussian(), control = list(adapt_delta = 0.95),thin=1, iter=4000, save_pars = save_pars(all = TRUE),data=LTcV[LTcV$th_family!= 'AF',])

# egg/la ratio
#d$LTEL <- d$base_temp_Egg-d$base_temp_Larvae
#d2$LTEL <- d2$base_temp_Egg-d2$base_temp_Larvae

#EL1 <- brm(LTEL ~ Body_Size + Host_Breadth  + Experiment + abs(Latitude) + tMode + (1|Order:Family) +(1|th_family), family=gaussian(), control = list(adapt_delta = 0.95),thin=1, iter=4000, save_pars = save_pars(all = TRUE),data=d[d$Length_measured == "Body",])
#EL2 <- brm(LTEL ~ Body_Size + Host_Breadth  + Experiment + abs(Latitude) + tMode + (1|Order:Family) +(1|th_family:th_genus), family=gaussian(), control = list(adapt_delta = 0.95),thin=1, iter=4000, save_pars = save_pars(all = TRUE),data=d2[d2$Length_measured == "Body",])
#stop('Model fitted')

########################
# Taxonomy model R2 calculation
########################
getR2 <- function(x){
  a <- bayes_R2(x, re.form = NULL)
  b <- bayes_R2(x, re.form = NA)
  
  print(paste('mar',b[1]))
  print(paste('con',a[1]))
}

lapply(list(LTE1,LTE2,DD1,DD2),getR2)



########################
# Phylogenetic models
########################

# General pre-process
tree <- read.tree('/home/zimo/Desktop/Zimo-PhenModels/Phylonegy/fruityTree2.nwk')
tree <- multi2di(tree)
#reformat scientific name | replace space with underscore
d$Scientific_name <- gsub(' ', '_', d$Scientific_name)
#Add a number in the scientific names of species with more than 1 model. In format Genus_species1. 
pop <- d %>% group_by(Scientific_name) %>% summarise(count=length(Scientific_name)) %>% filter(count!=1) %>% as.data.frame()
for(i in 1:length(pop$count)){
	subset <- d %>% filter(Scientific_name == pop$Scientific_name[i])
	for(j in 1:pop$count[i]){
		subset$Scientific_name[j] = paste(subset$Scientific_name[j],j,sep="")
			for (k in 1:length(subset$X)){
				d$Scientific_name[d$X == subset$X[k]] = subset$Scientific_name[k]
			}
	}
}


######## DD pre-process
info <- c("eggAdDD","Body_Size","Host_Breadth","absLatitude", "tMode", 'Experiment', 'th_family',
	   "Species")
dDD <- d[complete.cases(d[,info]),] %>% filter(Length_measured == "Body")

#Some finicky formating to get geiger's treedata() function to work
#It returns the parts of the phylogeny and dataset that overlap
rownames(dDD) <- dDD$Scientific_name
foo <- dDD
foo$Scientific_name <- NULL
fd <- treedata(tree, foo, warning=FALSE)

d_comp <- data.frame(fd$data) #the overlapping comparative trait data
tree_comp <- fd$phy #the overlapping phylogeny data
d_comp$Scientific_name <- rownames(d_comp)

# Formats are somehow messed up, change back to numeric/factor
d_comp$eggAdDD <- as.numeric(d_comp$eggAdDD)
d_comp$base_temp_eggAd <- as.numeric(d_comp$base_temp_eggAd)
d_comp$Body_Size <- as.numeric(d_comp$Body_Size)
d_comp$Host_Breadth <- as.numeric(d_comp$Host_Breadth)
d_comp$eggAdDD <- as.numeric(d_comp$eggAdDD)
d_comp$tMode <- as.factor(d_comp$tMode)
d_comp$tMode <- relevel(d_comp$tMode,ref="0")
d_comp$absLatitude <- as.numeric(d_comp$absLatitude)
d_comp$Experiment <- as.factor(d_comp$Experiment)
d_comp$Experiment <- relevel(d_comp$Experiment,ref="Field")
d_comp$th_family <- as.factor(d_comp$th_family)
d_comp$th_family <- relevel(d_comp$th_family,ref='Fabaceae')


######### LT Egg pre-process
info2 <- c("base_temp_Egg","Body_Size","Host_Breadth","absLatitude", "tMode", 'Experiment', 'th_family',
           "Species")
dLT <- d[complete.cases(d[,info2]),] %>% filter(Length_measured == "Body")

rownames(dLT) <- dLT$Scientific_name
foo2 <- dLT
foo2$Scientific_name <- NULL
fd2 <- treedata(tree, foo2, warning=FALSE)

d_comp2 <- data.frame(fd2$data) #the overlapping comparative trait data
tree_comp2 <- fd2$phy #the overlapping phylogeny data
d_comp2$Scientific_name <- rownames(d_comp2)

# Formats are somehow messed up, change back to numeric/factor
d_comp2$base_temp_Egg <- as.numeric(d_comp2$base_temp_Egg)
d_comp2$Body_Size <- as.numeric(d_comp2$Body_Size)
d_comp2$Host_Breadth <- as.numeric(d_comp2$Host_Breadth)
d_comp2$eggAdDD <- as.numeric(d_comp2$eggAdDD)
d_comp2$tMode <- as.factor(d_comp2$tMode)
d_comp2$tMode <- relevel(d_comp2$tMode,ref="0")
d_comp2$absLatitude <- as.numeric(d_comp2$absLatitude)
# d_comp2$Experiment <- as.factor(d_comp2$Experiment)
# d_comp2$Experiment <- relevel(d_comp2$Experiment,ref="Field")
d_comp2$th_family <- as.factor(d_comp2$th_family)
d_comp2$th_family <- relevel(d_comp2$th_family,ref='Fabaceae')




########## LT growth stages pre-process
# Remove incomplete rows and separate LT and DD data
info3 <- c("base_temp_Larvae","Body_Size","Host_Breadth","absLatitude", "tMode", 'Experiment', 'th_family',
	   "Species")
dLT2 <- d[complete.cases(d[,info3]),] %>% filter(Length_measured == "Body")

rownames(dLT2) <- dLT2$Scientific_name
foo3 <- dLT2
foo3$Scientific_name <- NULL
fd3 <- treedata(tree, foo3, warning=FALSE)

d_comp3 <- data.frame(fd3$data) #the overlapping comparative trait data
tree_comp3 <- fd3$phy #the overlapping phylogeny data
d_comp3$Scientific_name <- rownames(d_comp3)

# Formats are somehow messed up, change back to numeric/factor
d_comp3$base_temp_Larvae <- as.numeric(d_comp3$base_temp_Larvae)
d_comp3$Body_Size <- as.numeric(d_comp3$Body_Size)
d_comp3$Host_Breadth <- as.numeric(d_comp3$Host_Breadth)
d_comp3$eggAdDD <- as.numeric(d_comp3$eggAdDD)
d_comp3$tMode <- as.factor(d_comp3$tMode)
d_comp3$tMode <- relevel(d_comp3$tMode,ref="0")
d_comp3$absLatitude <- as.numeric(d_comp3$absLatitude)
d_comp3$Experiment <- as.factor(d_comp3$Experiment)
d_comp3$Experiment <- relevel(d_comp3$Experiment,ref="Field")
d_comp3$th_family <- as.factor(d_comp3$th_family)
d_comp3$th_family <- relevel(d_comp3$th_family,ref='Fabaceae')

############# Find max likelihood lambda value with AIC, with accuracy = 0.01, for each model
f <- function(x){
  PDD3 <- phyr::pglmm(eggAdDD ~ Body_Size + Host_Breadth + tMode + absLatitude + 
                        (1|Experiment) + (1|th_family) + (1|Scientific_name__), 
                      data = d_comp, family = "gaussian", cov_ranef = list(Scientific_name=vcv(corPagel(x,tree_comp))))
  print(paste(c(PDD3$logLik,PDD3$AIC)))
}

f2 <- function(x){
  PDD3 <- phyr::pglmm(base_temp_Egg ~ Body_Size + Host_Breadth + tMode + absLatitude + 
                        (1|th_family) + (1|Scientific_name__), 
                      data = d_comp2, family = "gaussian", cov_ranef = list(Scientific_name=vcv(corPagel(x,tree_comp2))))
  print(paste(c(PDD3$logLik,PDD3$AIC)))
}

f3 <- function(x){
  PDD3 <- phyr::pglmm(base_temp_Larvae ~ Body_Size + Host_Breadth + tMode + absLatitude + 
                        (1|Experiment) + (1|th_family) + (1|Scientific_name__), 
                      data = d_comp3, family = "gaussian", cov_ranef = list(Scientific_name=vcv(corPagel(x,tree_comp3))))
  print(paste(c(PDD3$logLik,PDD3$AIC)))
}

# lambda for PDD = 0.98
lambdalist <- seq(0,1,by=0.05)
lapply(lambdalist,f) # around 1.00
lambdalist2 <- seq(0.95,1.00,by=0.01)
lapply(lambdalist2,f) # 0.98

# lambda for PLTE = 0.86
lapply(lambdalist,f2) # around 0.85
lambdalist3 <- seq(0.81,0.89,by=0.01)
lapply(lambdalist3,f2) # 0.86

# lambda for PLTL = 0.86
lapply(lambdalist,f3) # around 0.90
lambdalist3 <- seq(0.86,0.94,by=0.01)
lapply(lambdalist3,f2) # 0.86

############## Fit models
PDD <- phyr::pglmm(eggAdDD ~ Body_Size + Host_Breadth + tMode + absLatitude + 
                      (1|Experiment) + (1|th_family) + (1|Scientific_name__), 
                    data = d_comp, family = "gaussian", cov_ranef = list(Scientific_name=vcv(corPagel(0.98,tree_comp))))

PLTE <- phyr::pglmm(base_temp_Egg ~ Body_Size + Host_Breadth + tMode + absLatitude + 
                      (1|th_family) + (1|Scientific_name__),  # No (1|Experiment) because the data subset only contains one level
                    data = d_comp2, family = "gaussian", cov_ranef = list(Scientific_name=vcv(corPagel(0.86,tree_comp2))))

PLTL <- phyr::pglmm(base_temp_Larvae ~ Body_Size + Host_Breadth + tMode + absLatitude + 
                      (1|Experiment) + (1|th_family) + (1|Scientific_name__), 
                    data = d_comp3, family = "gaussian", cov_ranef = list(Scientific_name=vcv(corPagel(0.86,tree_comp3))))
phyr::pglmm_plot_ranef(eggAdDD ~ Body_Size + Host_Breadth + tMode + absLatitude + 
              (1|Experiment) + (1|th_family) + (1|Scientific_name__), 
            data = d_comp, family = "gaussian", cov_ranef = list(Scientific_name=tree_comp), sp.var="Scientific_name",site.var="Scientific_name")

# Calculate R2 for phylogy models
rr2::R2_lik(PDD)
rr2::R2_lik(PLTE)
rr2::R2_lik(PLTL)

###############################
dm <- read.csv('/home/zimo/Desktop/Zimo-PhenModels/phenoDBherbiesDataPack/phenoDBherbies.4.5.csv')

dm <- dm[dm$lat.range <= 10,]
dm <- dm[dm$Estimated_Location != "E",]
map<-map_data("world")

Lat <- as.numeric(dm$Latitude)
dm$Latitude <- as.numeric(dm$Latitude)

# Chi-square for N vs S hemisphere. 
North <- length(na.omit(dm$Latitude[dm$Latitude>0]))
South <- length(na.omit(dm$Latitude[dm$Latitude<0])) 
ActualNS <- c(North, South)
ExpNS <- c(0.5,0.5)
NSbias <- chisq.test(ActualNS,p=ExpNS)
NSbias
# Chi-square for climate zones
tr <- length(na.omit(dm$Latitude[abs(dm$Latitude)<23.5]))
st <- length(na.omit(dm$Latitude[abs(dm$Latitude)>23.5 & abs(dm$Latitude)<40]))
te <- length(na.omit(dm$Latitude[abs(dm$Latitude)>40 & abs(dm$Latitude)<60]))
co <- length(na.omit(dm$Latitude[abs(dm$Latitude)>60]))
ExpClimate <- c(1/3,1/3,1/3)
ActualClimate <- c(tr,(st+te),co)
Climatebias <- chisq.test(ActualClimate,p=ExpClimate)
Climatebias







