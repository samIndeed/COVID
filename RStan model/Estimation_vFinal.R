rm(list=setdiff(ls(), "filepath"))

################################################################################################################
################### 1 # FRANCE fitting  ########################################################################
################################################################################################################
# POPULATION DISTRIBUTION FROM SALJE et al 2020, SCIENCE
popData<-read.table(file="frPopPyramid.csv", header=TRUE, sep=",")
names(popData)=c("age", "male", "female", "total")

tmpB=rbind(colSums(popData[popData$age==0:19,]),
           colSums(popData[popData$age==20:29,]),
           colSums(popData[popData$age==30:39,]),
           colSums(popData[popData$age==40:49,]),
           colSums(popData[popData$age==50:59,]),
           colSums(popData[popData$age==60:69,]),
           colSums(popData[popData$age==70:79,]),
           ( colSums(popData[popData$age==80:99,])+popData[popData$age==100,]))
tmpB2=tmpB[,c(2:4)]
tmpB3=c("0-19","20-29","30-39","40-49","50-59","60-69","70-79",">=80")
francePopDataB=cbind(tmpB3,tmpB2)
ageDis_frS=francePopDataB[,4]

# mean age in each age group
avrAges=c(sum((popData[popData$age==0:19,4])*(0:19))/sum(popData[popData$age==0:19,4]),
          sum((popData[popData$age==20:29,4])*(20:29))/sum(popData[popData$age==20:29,4]),
          sum((popData[popData$age==30:39,4])*(30:39))/sum(popData[popData$age==30:39,4]),
          sum((popData[popData$age==40:49,4])*(40:49))/sum(popData[popData$age==40:49,4]),
          sum((popData[popData$age==50:59,4])*(50:59))/sum(popData[popData$age==50:59,4]),
          sum((popData[popData$age==60:69,4])*(60:69))/sum(popData[popData$age==60:69,4]),
          sum((popData[popData$age==70:79,4])*(70:79))/sum(popData[popData$age==70:79,4]),
          sum((popData[popData$age==80:99,4])*(80:99)+(popData[popData$age==100,4])*100)/sum(popData[popData$age==80:99,4]+popData[popData$age==100,4]))

###########################################################################################
# CORRECTED DAILY AGE SPECIFIC HOSPITAL ADMISSIONS SALJE et al 2020, SCIENCE
dailyHospCounts_allReg<-read.table(file="dailyHospCounts_allReg.csv", header=TRUE, sep=",")  # Salje's dailyHospCounts.csv
numD=dim(dailyHospCounts_allReg)[1]
dHospData=rbind(dailyHospCounts_allReg[,2],dailyHospCounts_allReg[,3],dailyHospCounts_allReg[,4],dailyHospCounts_allReg[,5],dailyHospCounts_allReg[,6],dailyHospCounts_allReg[,7],dailyHospCounts_allReg[,8],dailyHospCounts_allReg[,9])
cumdHospData=rbind(cumsum(dHospData[1,]),cumsum(dHospData[2,]),cumsum(dHospData[3,]),cumsum(dHospData[4,]),cumsum(dHospData[5,]),cumsum(dHospData[6,]),cumsum(dHospData[7,]),cumsum(dHospData[8,]))
###########################################################################################

###########################################################################################
# CORRECTED HOSPITAL REMOVALS SALJE et al 2020, SCIENCE
dailySIVIC_allReg<-read.table(file="SIVIC/SIVIC_daily_numbers_region_corrected_histo_20200508.csv", header=TRUE, sep=",")
currentSIVIC_allReg<-read.table(file="SIVIC/SIVIC_total_numbers_region_corrected_histo_20200508.csv", header=TRUE, sep=",")
regs=unique(dailySIVIC_allReg$region)
regNewHospSIVIC=numeric(0)
regCurrHospSIVIC=numeric(0)
for (bb in regs) {
  dates=currentSIVIC_allReg[currentSIVIC_allReg$region==bb,2]
  regNewHospSIVIC=rbind(regNewHospSIVIC,dailySIVIC_allReg[dailySIVIC_allReg$region==bb,4])
  regCurrHospSIVIC=rbind(regCurrHospSIVIC,currentSIVIC_allReg[currentSIVIC_allReg$region==bb,4])
}
colSums(regNewHospSIVIC)
CumHospSIVIC=cumsum(colSums(regNewHospSIVIC))                   # SIVIC CUMULATIVE
propRemovSIVIC=1-(colSums(regCurrHospSIVIC)/CumHospSIVIC)       # SIVIC REMOVED
###########################################################################################
propRecover_frS=propRemovSIVIC
CumInfSevere_frS=cumdHospData

############################################################
#library(socialmixr)
#frSurvey=get_survey("https://doi.org/10.5281/zenodo.1157918")
#lower.age.limit=c(0,20,30,40,50,60,70,80)
#popdf=data.frame(lower.age.limit,population2)    # NOTE this is same as ageDis_frS
#cm3=contact_matrix(frSurvey, countries = "France",survey.pop=popdf, age.limits = c(0, 20, 30, 40, 50, 60, 70, 80))
# HARDCODE CM in case of zenodo failure:
# $matrix
# contact.age.group
# age.group    [0,20)  [20,30)  [30,40)  [40,50)  [50,60)  [60,70)   [70,80)       80+
#   [0,20)  13.147714 1.966002 3.992966 4.021102 1.711606 1.029308 0.4701055 0.1606096
# [20,30)  2.988095 8.410714 3.571429 3.380952 2.952381 0.922619 0.4226190 0.3630952
# [30,40)  3.188034 3.487179 6.000000 5.094017 2.692308 1.623932 0.5470085 0.4273504
# [40,50)  3.130435 3.539130 4.869565 6.460870 5.139130 2.000000 0.8434783 0.1739130
# [50,60)  2.118280 3.870968 4.758065 5.333333 6.360215 2.483871 0.7258065 0.8172043
# [60,70)  2.172805 1.600567 3.303116 3.124646 3.121813 4.983003 1.6402266 0.9943343
# [70,80)  1.182927 1.176829 2.158537 3.000000 2.060976 3.201220 2.9329268 1.1280488
# 80+      1.072727 1.127273 1.672727 2.309091 2.909091 1.654545 1.7454545 1.5272727
############################################################
# Entering Hard-coded socialmixr numbers
tmpMatrix = matrix( 
  +   c(13.147714, 2.988095, 3.188034, 3.130435, 2.118280, 2.172805, 1.182927, 1.072727,1.966002, 8.410714, 3.487179, 3.539130, 3.870968, 1.600567, 1.176829, 1.127273,3.992966, 3.571429, 6.000000, 4.869565, 4.758065, 3.303116, 2.158537, 1.672727,4.021102, 3.380952, 5.094017, 6.460870, 5.333333, 3.124646, 3.000000, 2.309091,1.711606, 2.952381, 2.692308, 5.139130, 6.360215, 3.121813, 2.060976, 2.909091,1.029308, 0.922619, 1.623932, 2.000000, 2.483871, 4.983003, 3.201220, 1.654545,0.4701055, 0.4226190, 0.5470085, 0.8434783, 0.7258065, 1.6402266, 2.9329268, 1.7454545,0.1606096, 0.3630952, 0.4273504, 0.1739130, 0.8172043, 0.9943343, 1.1280488, 1.5272727), 
  +   8, 
  +   8) 
cmIn=tmpMatrix
#cmIn=cm3$matrix

library(rstan)
options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
ages=c(9.752804,24.495934,34.583572,44.634635,54.469143,64.418965,73.794358,87.23397)
baseNdays=24;
data <- list( D_NDays=baseNdays,
              D_NAgeGroups=8,
              D_cumulHospByAge=CumInfSevere_frS,
              D_ageDisIn=ageDis_frS,
              D_meanAges=ages,
              D_removedsByAge_frac=propRecover_frS,
              D_contactMPre=cmIn,
              D_contactMPost=cmIn,
              D_runType=0,
              D_numExtraDays=0,
              D_negTau=0,
              D_posTau=0,
              D_meancontact=mean(rowSums(cmIn))
)
warmUpT=1000
longRunT=10000; 
adaptVal=0.99999999999
sink(file='newRecords.txt')
########################################################################################
# There are 4 sets of results representing Bayesian fitting of 4 variants of a single model
########################################################################################
# Each fit uses the model file ageDepNB_vFinal.stan
# Move between the 4 variants in the .stan file using the data$D_runType option \in [1,2,3,4]
########################################################################################

########################################################################################
# A: if plotting Psevere fitted individually with no exp assumption: INITIAL MODEL
########################################################################################
#includeA=0;
#if (includeA) {

data$D_runType=0
print('A')
initf1 <- function() list(relSev_inv=rep(2,8),transMC = 0.999, transSC = 0.999, transMU = 0.999, transSU = 0.999, transMU2 = 0.999, transSU2 = 0.999, phi_inv = 1.001) #
fit_A <- stan(file = "ageDepNB_vFinal.stan",
             data=data,
             init=initf1,
             iter=longRunT,
             warmup = warmUpT,
             chains=4,
             control = list(max_treedepth = 15,adapt_delta=adaptVal) 
) 
smA=summary(fit_A)
smmA=smA$summary
chainyA_NoPerm<- extract(fit_A, permuted = FALSE,inc_warmup=FALSE)
chainyA_YesPerm<- extract(fit_A, permuted = TRUE)
save(chainyA_NoPerm, file = "vFinal_stanA_NP.RData")
save(chainyA_YesPerm, file = "vFinal_stanA_YP.RData")
save(smmA, file = "vFinal_stanASummary.RData")

################# DIAGNOSTICS #######################################################
namesList=colnames(head(chainyA_NoPerm[1,,]))
tmpny=which(is.na(smmA[,10])) 
invalidnamesList=namesList[tmpny]   # CHECK these NA's from the summary - were they uniform? Then discard uniforms for diagnostics...
# Add to newnamesList: contacts_vec[1..8] SevereInc[1..8] both of which are uniform in the sense of not involving fitted parameters...
invalidnamesList=c(invalidnamesList,"contacts_vec[1]","contacts_vec[2]","contacts_vec[3]","contacts_vec[4]","contacts_vec[5]","contacts_vec[6]","contacts_vec[7]","contacts_vec[8]")
invalidnamesList=c(invalidnamesList,"SevereInc[1]","SevereInc[2]","SevereInc[3]","SevereInc[4]","SevereInc[5]","SevereInc[6]","SevereInc[7]","SevereInc[8]")
invalidnamesList=c(invalidnamesList,"lp__")
Rhatvec=numeric(0)
Bulkvec=numeric(0)
Tailvec=numeric(0)
newnamesList=colnames(chainyA_NoPerm[,1,])
numNames=length(colnames(chainyA_NoPerm[,1,]))
for (tt in 1:(length(numNames))) {
  if (sum(newnamesList[tt]==invalidnamesList)==0) {
      #print(newnamesList[tt])
      tmpVal_Rhat=Rhat(chainyA_NoPerm[,,colnames(chainyA_NoPerm[,1,])==newnamesList[tt]])
      #print(tmpVal_Rhat)
      tmpVal_bulk=ess_bulk(chainyA_NoPerm[,,colnames(chainyA_NoPerm[,1,])==newnamesList[tt]])
      #print(tmpVal_bulk)
      tmpVal_tail=ess_tail(chainyA_NoPerm[,,colnames(chainyA_NoPerm[,1,])==newnamesList[tt]])
      #print(tmpVal_tail)
      Rhatvec=c(Rhatvec,tmpVal_Rhat)
      Bulkvec=c(Bulkvec,tmpVal_bulk)
      Tailvec=c(Tailvec,tmpVal_tail)
  }
}
print(max(Rhatvec))  # Check this is <1.05
print(min(Bulkvec))  # Check this is >100
print(min(Tailvec))  # Check this is >100
##############################################################################################################
#}

########################################################################################
# B: if plotting Psevere under exp assumption with individual cohort deviations : MAIN MODEL
# stan file is set up so that epidemic prediction is only done in B
########################################################################################
data <- list( D_NDays=baseNdays,
              D_NAgeGroups=8,
              D_cumulHospByAge=CumInfSevere_frS,
              D_ageDisIn=ageDis_frS,
              D_meanAges=ages,
              D_removedsByAge_frac=propRecover_frS,
              D_contactMPre=cmIn,
              D_contactMPost=cmIn,
              D_runType=1,
              D_numExtraDays=0,
              D_negTau=0,
              D_posTau=0,
              D_meancontact=mean(rowMeans(cmIn))
)

#includeB=1;
#if (includeB) {
  
data$D_runType=1   # This will trigger conditional declaration of parameters in the stan file for the second run variant (i.e. slope with deviations)
initB_f1 <- function() list(ageExp=array(0.001,dim = 1),coefExp_inv=array(2,dim = 1),deviationFromExp=array(0,dim = 8),transMC = 0.999, transSC = 0.999, transMU = 0.999, transSU = 0.999, phi_inv = 1.001) #
print('B')
fit_B <- stan(file="ageDepNB_vFinal.stan", 
             data=data,
             init=initB_f1,
             iter=longRunT,
             chains=4,
             warmup = warmUpT,
             control = list(max_treedepth = 15,adapt_delta=adaptVal)
) 
smB=summary(fit_B)
smmB=smB$summary
chainyB_NoPerm<- extract(fit_B, permuted = FALSE,inc_warmup=FALSE)
chainyB_YesPerm<- extract(fit_B, permuted = TRUE)
save(chainyB_NoPerm, file = "vFinal_stanB_NP.RData")
save(chainyB_YesPerm, file = "vFinal_stanB_YP.RData")
save(smmB, file = "vFinal_stanBSummary.RData")

p=stan_hist(fit_B, pars =c("ageExp[1]"),bins=50)
p+geom_vline(xintercept=c(summary(fit_B, pars = "ageExp[1]")$summary[4],summary(fit_B, pars = "ageExp[1]")$summary[8]),color = "black", size=1.5) + 
  geom_point(aes(x=0.044,y=0),color = "black", size=4,shape=23,fill="yellow",stroke = 2) +
  labs(x = "Posterior m (exponential function)")+
  theme(axis.text.x = element_text(face="bold", size=16))+
  xlim(0.0375,0.06)


################# DIAGNOSTICS #######################################################
namesList=colnames(head(chainyB_NoPerm[1,,]))
tmpny=which(is.na(smmB[,10])) 
invalidnamesList=namesList[tmpny]   # CHECK these NA's from the summary - were they uniform? Then discard uniforms for diagnostics...
# Add to newnamesList: contacts_vec[1..8] SevereInc[1..8] both of which are uniform in the sense of not involving fitted parameters...
invalidnamesList=c(invalidnamesList,"contacts_vec[1]","contacts_vec[2]","contacts_vec[3]","contacts_vec[4]","contacts_vec[5]","contacts_vec[6]","contacts_vec[7]","contacts_vec[8]")
invalidnamesList=c(invalidnamesList,"SevereInc[1]","SevereInc[2]","SevereInc[3]","SevereInc[4]","SevereInc[5]","SevereInc[6]","SevereInc[7]","SevereInc[8]")
invalidnamesList=c(invalidnamesList,"lp__")
Rhatvec=numeric(0)
Bulkvec=numeric(0)
Tailvec=numeric(0)
newnamesList=colnames(chainyB_NoPerm[,1,])
numNames=length(colnames(chainyB_NoPerm[,1,]))
for (tt in 1:(length(numNames))) {
  if (sum(newnamesList[tt]==invalidnamesList)==0) {
    #print(newnamesList[tt])
    tmpVal_Rhat=Rhat(chainyB_NoPerm[,,colnames(chainyB_NoPerm[,1,])==newnamesList[tt]])
    #print(tmpVal_Rhat)
    tmpVal_bulk=ess_bulk(chainyB_NoPerm[,,colnames(chainyB_NoPerm[,1,])==newnamesList[tt]])
    #print(tmpVal_bulk)
    tmpVal_tail=ess_tail(chainyB_NoPerm[,,colnames(chainyB_NoPerm[,1,])==newnamesList[tt]])
    #print(tmpVal_tail)
    Rhatvec=c(Rhatvec,tmpVal_Rhat)
    Bulkvec=c(Bulkvec,tmpVal_bulk)
    Tailvec=c(Tailvec,tmpVal_tail)
  }
}
print(max(Rhatvec))  # Check this is <1.05
print(min(Bulkvec))  # Check this is >100
print(min(Tailvec))  # Check this is >100
##############################################################################################################
#}

#######################################################################################################################
#######################################################################################################################
# MODELS WITH DEVIATIONS (for sensitivity analysis/exploration of factors underlying low juvenile hospitalisations) ###
#######################################################################################################################
#######################################################################################################################


########################################################################################
# D: if extending exp assumption with individual cohort deviations to  #################
# also include deviations in infection|exposure                        #################
########################################################################################
data$D_runType=1   # triggers conditional declaration of parameters in stan file 
initB_f1 <- function() list(ageExp=array(0.001,dim = 1),coefExp_inv=array(2,dim = 1),deviationFromExp=array(0,dim = 8),deviationFromUE=array(0,dim = 8),transMC = 0.999, transSC = 0.999, transMU = 0.999, transSU = 0.999, phi_inv = 1.001) #

print('D1')
fitD_dev2 <- stan(file = "ageDepNB_vFinal.stan",
              data=data,
              init=initB_f1,
              iter=longRunT,
              warmup = warmUpT,
              chains=4,
              control = list(max_treedepth = 15,adapt_delta=adaptVal) 
) 
smD2=summary(fitD_dev2)
smmD2=smD2$summary
chainyD2_NoPerm<- extract(fitD_dev2, permuted = FALSE,inc_warmup=FALSE)
chainyD2_YesPerm<- extract(fitD_dev2, permuted = TRUE)
save(chainyD2_NoPerm, file = "vFinal_stanD2_NP.RData")
save(chainyD2_YesPerm, file = "vFinal_stanD2_YP.RData")
save(smmD2, file = "vFinal_stanD2Summary.RData")

################# DIAGNOSTICS #######################################################
namesList=colnames(head(chainyD2_NoPerm[1,,]))
tmpny=which(is.na(smmD2[,10])) 
invalidnamesList=namesList[tmpny]   # CHECK these NA's from the summary - were they uniform? Then discard uniforms for diagnostics...
# Add to newnamesList: contacts_vec[1..8] SevereInc[1..8] both of which are uniform in the sense of not involving fitted parameters...
invalidnamesList=c(invalidnamesList,"contacts_vec[1]","contacts_vec[2]","contacts_vec[3]","contacts_vec[4]","contacts_vec[5]","contacts_vec[6]","contacts_vec[7]","contacts_vec[8]")
invalidnamesList=c(invalidnamesList,"SevereInc[1]","SevereInc[2]","SevereInc[3]","SevereInc[4]","SevereInc[5]","SevereInc[6]","SevereInc[7]","SevereInc[8]")
invalidnamesList=c(invalidnamesList,"lp__")
Rhatvec=numeric(0)
Bulkvec=numeric(0)
Tailvec=numeric(0)
newnamesList=colnames(chainyD2_NoPerm[,1,])
numNames=length(colnames(chainyD2_NoPerm[,1,]))
for (tt in 1:(length(numNames))) {
  if (sum(newnamesList[tt]==invalidnamesList)==0) {
    #print(newnamesList[tt])
    tmpVal_Rhat=Rhat(chainyD2_NoPerm[,,colnames(chainyD2_NoPerm[,1,])==newnamesList[tt]])
    #print(tmpVal_Rhat)
    tmpVal_bulk=ess_bulk(chainyD2_NoPerm[,,colnames(chainyD2_NoPerm[,1,])==newnamesList[tt]])
    #print(tmpVal_bulk)
    tmpVal_tail=ess_tail(chainyD2_NoPerm[,,colnames(chainyD2_NoPerm[,1,])==newnamesList[tt]])
    #print(tmpVal_tail)
    Rhatvec=c(Rhatvec,tmpVal_Rhat)
    Bulkvec=c(Bulkvec,tmpVal_bulk)
    Tailvec=c(Tailvec,tmpVal_tail)
  }
}
print(max(Rhatvec))  # Check this is <1.05
print(min(Bulkvec))  # Check this is >100
print(min(Tailvec))  # Check this is >100
##############################################################################################################





########################################################################################
# D: if extending exp assumption with individual cohort deviations to  #################
# also include deviations in infectiousness of infecteds, first of 2 types!     ########
########################################################################################
data$D_runType=3   # triggers conditional declaration of parameters in stan file 
initB_f1 <- function() list(ageExp=array(0.001,dim = 1),coefExp_inv=array(2,dim = 1),deviationFromUT=array(0,dim = 8),deviationFromExp=array(0,dim = 8),transMC = 0.999, transSC = 0.999, transMU = 0.999, transSU = 0.999, phi_inv = 1.001) #
print('D2')
fitD_dev3 <- stan(file = "ageDepNB_vFinal.stan",
              data=data,
              init=initB_f1,
              iter=longRunT, 
              warmup = warmUpT,
              chains=4,
              control = list(max_treedepth = 15,adapt_delta=adaptVal)  
) 
smD3=summary(fitD_dev3)
smmD3=smD3$summary
chainyD3_NoPerm<- extract(fitD_dev3, permuted = FALSE,inc_warmup=FALSE)
chainyD3_YesPerm<- extract(fitD_dev3, permuted = TRUE)
save(chainyD3_NoPerm, file = "vFinal_stanD3_NP.RData")
save(chainyD3_YesPerm, file = "vFinal_stanD3_YP.RData")
save(smmD3, file = "vFinal_stanD3Summary.RData")

################# DIAGNOSTICS #######################################################
namesList=colnames(head(chainyD3_NoPerm[1,,]))
tmpny=which(is.na(smmD3[,10])) 
invalidnamesList=namesList[tmpny]   # CHECK these NA's from the summary - were they uniform? Then discard uniforms for diagnostics...
# Add to newnamesList: contacts_vec[1..8] SevereInc[1..8] both of which are uniform in the sense of not involving fitted parameters...
invalidnamesList=c(invalidnamesList,"contacts_vec[1]","contacts_vec[2]","contacts_vec[3]","contacts_vec[4]","contacts_vec[5]","contacts_vec[6]","contacts_vec[7]","contacts_vec[8]")
invalidnamesList=c(invalidnamesList,"SevereInc[1]","SevereInc[2]","SevereInc[3]","SevereInc[4]","SevereInc[5]","SevereInc[6]","SevereInc[7]","SevereInc[8]")
invalidnamesList=c(invalidnamesList,"lp__")
Rhatvec=numeric(0)
Bulkvec=numeric(0)
Tailvec=numeric(0)
newnamesList=colnames(chainyD3_NoPerm[,1,])
numNames=length(colnames(chainyD3_NoPerm[,1,]))
for (tt in 1:(length(numNames))) {
  if (sum(newnamesList[tt]==invalidnamesList)==0) {
    #print(newnamesList[tt])
    tmpVal_Rhat=Rhat(chainyD3_NoPerm[,,colnames(chainyD3_NoPerm[,1,])==newnamesList[tt]])
    #print(tmpVal_Rhat)
    tmpVal_bulk=ess_bulk(chainyD3_NoPerm[,,colnames(chainyD3_NoPerm[,1,])==newnamesList[tt]])
    #print(tmpVal_bulk)
    tmpVal_tail=ess_tail(chainyD3_NoPerm[,,colnames(chainyD3_NoPerm[,1,])==newnamesList[tt]])
    #print(tmpVal_tail)
    Rhatvec=c(Rhatvec,tmpVal_Rhat)
    Bulkvec=c(Bulkvec,tmpVal_bulk)
    Tailvec=c(Tailvec,tmpVal_tail)
  }
}
print(max(Rhatvec))  # Check this is <1.05
print(min(Bulkvec))  # Check this is >100
print(min(Tailvec))  # Check this is >100
##############################################################################################################



########################################################################################
# D: if extending exp assumption with individual cohort deviations to  #################
# also include deviations in infectiousness of infecteds, second of 2 types!     #######
########################################################################################
data$D_runType=4   # trigger conditional declaration of parameters in stan file 
initB_f1 <- function() list(ageExp=array(0.001,dim = 1),coefExp_inv=array(2,dim = 1),deviationFromUT=array(0,dim = 8),deviationFromUT_S=array(0,dim = 8),deviationFromExp=array(0,dim = 8),transMC = 0.999, transSC = 0.999, transMU = 0.999, transSU = 0.999, phi_inv = 1.001) #
print('D3')
fitD_dev4 <- stan(file = "ageDepNB_vFinal.stan",
              data=data,
              init=initB_f1,
              iter=longRunT, 
              warmup = warmUpT,
              chains=4,
              control = list(max_treedepth = 15,adapt_delta=adaptVal) 
) 
smD4=summary(fitD_dev4)
smmD4=smD4$summary
chainyD4_NoPerm<- extract(fitD_dev4, permuted = FALSE,inc_warmup=FALSE)
chainyD4_YesPerm<- extract(fitD_dev4, permuted = TRUE)
save(chainyD4_NoPerm, file = "vFinal_stanD4_NP.RData")
save(chainyD4_YesPerm, file = "vFinal_stanD4_YP.RData")
save(smmD4, file = "vFinal_stanD4Summary.RData")

################# DIAGNOSTICS #######################################################
namesList=colnames(head(chainyD4_NoPerm[1,,]))
tmpny=which(is.na(smmD4[,10])) 
invalidnamesList=namesList[tmpny]   # CHECK these NA's from the summary - were they uniform? Then discard uniforms for diagnostics...
# Add to newnamesList: contacts_vec[1..8] SevereInc[1..8] both of which are uniform in the sense of not involving fitted parameters...
invalidnamesList=c(invalidnamesList,"contacts_vec[1]","contacts_vec[2]","contacts_vec[3]","contacts_vec[4]","contacts_vec[5]","contacts_vec[6]","contacts_vec[7]","contacts_vec[8]")
invalidnamesList=c(invalidnamesList,"SevereInc[1]","SevereInc[2]","SevereInc[3]","SevereInc[4]","SevereInc[5]","SevereInc[6]","SevereInc[7]","SevereInc[8]")
invalidnamesList=c(invalidnamesList,"lp__")
Rhatvec=numeric(0)
Bulkvec=numeric(0)
Tailvec=numeric(0)
newnamesList=colnames(chainyD4_NoPerm[,1,])
numNames=length(colnames(chainyD4_NoPerm[,1,]))
for (tt in 1:(length(numNames))) {
  if (sum(newnamesList[tt]==invalidnamesList)==0) {
    #print(newnamesList[tt])
    tmpVal_Rhat=Rhat(chainyD4_NoPerm[,,colnames(chainyD4_NoPerm[,1,])==newnamesList[tt]])
    #print(tmpVal_Rhat)
    tmpVal_bulk=ess_bulk(chainyD4_NoPerm[,,colnames(chainyD4_NoPerm[,1,])==newnamesList[tt]])
    #print(tmpVal_bulk)
    tmpVal_tail=ess_tail(chainyD4_NoPerm[,,colnames(chainyD4_NoPerm[,1,])==newnamesList[tt]])
    #print(tmpVal_tail)
    Rhatvec=c(Rhatvec,tmpVal_Rhat)
    Bulkvec=c(Bulkvec,tmpVal_bulk)
    Tailvec=c(Tailvec,tmpVal_tail)
  }
}
print(max(Rhatvec))  # Check this is <1.05
print(min(Bulkvec))  # Check this is >100
print(min(Tailvec))  # Check this is >100
##############################################################################################################


########################################################################################
# E: another sensitivity a. of MAIN MODEL (i.e. B above again), this time with total#days data +-delta
########################################################################################
deltaNdays=1;
data <- list( D_NDays=baseNdays+deltaNdays,
              D_NAgeGroups=8,
              D_cumulHospByAge=CumInfSevere_frS,
              D_ageDisIn=ageDis_frS,
              D_meanAges=ages,
              D_removedsByAge_frac=propRecover_frS,
              D_contactMPre=cmIn,
              D_contactMPost=cmIn,
              D_runType=1,
              D_numExtraDays=0,
              D_negTau=0,
              D_posTau=0,
              D_meancontact=mean(rowSums(cmIn))
)
initB_f1 <- function() list(ageExp=array(0.001,dim = 1),coefExp_inv=array(2,dim = 1),deviationFromExp=array(0,dim = 8),transMC = 0.999, transSC = 0.999, transMU = 0.999, transSU = 0.999, phi_inv = 1.001) #
print('E')
fitE_extra1 <- stan(file = "ageDepNB_vFinal.stan",
              data=data,
              init=initB_f1,
              iter=longRunT, 
              warmup = warmUpT,
              chains=4,
              control = list(max_treedepth = 15,adapt_delta=adaptVal) 
)
summary(fitE_extra1, pars = "ageExp[1]")$summary[c(4,8)]
summary(fitE_extra1, pars = "coefExp[1]")$summary[c(4,8)]
smE1=summary(fitE_extra1)
smmE1=smE1$summary
chainyE1_NoPerm<- extract(fitE_extra1, permuted = FALSE,inc_warmup=FALSE)
chainyE1_YesPerm<- extract(fitE_extra1, permuted = TRUE)
save(chainyE1_NoPerm, file = "vFinal_stanE1_NP.RData")
save(chainyE1_YesPerm, file = "vFinal_stanE1_YP.RData")
save(smmE1, file = "vFinal_stanE1Summary.RData")
################################

################# DIAGNOSTICS #######################################################
namesList=colnames(head(chainyE1_NoPerm[1,,]))
tmpny=which(is.na(smmE1[,10])) 
invalidnamesList=namesList[tmpny]   # CHECK these NA's from the summary - were they uniform? Then discard uniforms for diagnostics...
# Add to newnamesList: contacts_vec[1..8] SevereInc[1..8] both of which are uniform in the sense of not involving fitted parameters...
invalidnamesList=c(invalidnamesList,"contacts_vec[1]","contacts_vec[2]","contacts_vec[3]","contacts_vec[4]","contacts_vec[5]","contacts_vec[6]","contacts_vec[7]","contacts_vec[8]")
invalidnamesList=c(invalidnamesList,"SevereInc[1]","SevereInc[2]","SevereInc[3]","SevereInc[4]","SevereInc[5]","SevereInc[6]","SevereInc[7]","SevereInc[8]")
invalidnamesList=c(invalidnamesList,"lp__")
Rhatvec=numeric(0)
Bulkvec=numeric(0)
Tailvec=numeric(0)
newnamesList=colnames(chainyE1_NoPerm[,1,])
numNames=length(colnames(chainyE1_NoPerm[,1,]))
for (tt in 1:(length(numNames))) {
  if (sum(newnamesList[tt]==invalidnamesList)==0) {
    #print(newnamesList[tt])
    tmpVal_Rhat=Rhat(chainyE1_NoPerm[,,colnames(chainyE1_NoPerm[,1,])==newnamesList[tt]])
    #print(tmpVal_Rhat)
    tmpVal_bulk=ess_bulk(chainyE1_NoPerm[,,colnames(chainyE1_NoPerm[,1,])==newnamesList[tt]])
    #print(tmpVal_bulk)
    tmpVal_tail=ess_tail(chainyE1_NoPerm[,,colnames(chainyE1_NoPerm[,1,])==newnamesList[tt]])
    #print(tmpVal_tail)
    Rhatvec=c(Rhatvec,tmpVal_Rhat)
    Bulkvec=c(Bulkvec,tmpVal_bulk)
    Tailvec=c(Tailvec,tmpVal_tail)
  }
}
print(max(Rhatvec))  # Check this is <1.05
print(min(Bulkvec))  # Check this is >100
print(min(Tailvec))  # Check this is >100
##############################################################################################################




########################################################################################
# F: another sensitivity a. of MAIN MODEL (i.e. B above again), this time with total#days data +-delta
########################################################################################
deltaNdays=-1;
data <- list( D_NDays=baseNdays+deltaNdays,
              D_NAgeGroups=8,
              D_cumulHospByAge=CumInfSevere_frS,
              D_ageDisIn=ageDis_frS,
              D_meanAges=ages,
              D_removedsByAge_frac=propRecover_frS,
              D_contactMPre=cmIn,
              D_contactMPost=cmIn,
              D_runType=1,
              D_numExtraDays=0,
              D_negTau=0,
              D_posTau=0,
              D_meancontact=mean(rowSums(cmIn))
)
initB_f1 <- function() list(ageExp=array(0.001,dim = 1),coefExp_inv=array(2,dim = 1),deviationFromExp=array(0,dim = 8),transMC = 0.999, transSC = 0.999, transMU = 0.999, transSU = 0.999, phi_inv = 1.001) #
print('F')
fitE_extra2 <- stan(file = "ageDepNB_vFinal.stan",
              data=data,
              init=initB_f1,
              iter=longRunT, 
              chains=4,
              warmup = warmUpT,
              control = list(max_treedepth = 15,adapt_delta=adaptVal) 
) 
summary(fitE_extra2, pars = "ageExp[1]")$summary[c(4,8)]
summary(fitE_extra2, pars = "coefExp[1]")$summary[c(4,8)]
smE2=summary(fitE_extra2)
smmE2=smE2$summary
chainyE2_NoPerm<- extract(fitE_extra2, permuted = FALSE,inc_warmup=FALSE)
chainyE2_YesPerm<- extract(fitE_extra2, permuted = TRUE)
save(chainyE2_NoPerm, file = "vFinal_stanE2_NP.RData")
save(chainyE2_YesPerm, file = "vFinal_stanE2_YP.RData")
save(smmE2, file = "vFinal_stanE2Summary.RData")
################################

################# DIAGNOSTICS #######################################################
namesList=colnames(head(chainyE2_NoPerm[1,,]))
tmpny=which(is.na(smmE2[,10])) 
invalidnamesList=namesList[tmpny]   # CHECK these NA's from the summary - were they uniform? Then discard uniforms for diagnostics...
# Add to newnamesList: contacts_vec[1..8] SevereInc[1..8] both of which are uniform in the sense of not involving fitted parameters...
invalidnamesList=c(invalidnamesList,"contacts_vec[1]","contacts_vec[2]","contacts_vec[3]","contacts_vec[4]","contacts_vec[5]","contacts_vec[6]","contacts_vec[7]","contacts_vec[8]")
invalidnamesList=c(invalidnamesList,"SevereInc[1]","SevereInc[2]","SevereInc[3]","SevereInc[4]","SevereInc[5]","SevereInc[6]","SevereInc[7]","SevereInc[8]")
invalidnamesList=c(invalidnamesList,"lp__")
Rhatvec=numeric(0)
Bulkvec=numeric(0)
Tailvec=numeric(0)
newnamesList=colnames(chainyE2_NoPerm[,1,])
numNames=length(colnames(chainyE2_NoPerm[,1,]))
for (tt in 1:(length(numNames))) {
  if (sum(newnamesList[tt]==invalidnamesList)==0) {
    #print(newnamesList[tt])
    tmpVal_Rhat=Rhat(chainyE2_NoPerm[,,colnames(chainyE2_NoPerm[,1,])==newnamesList[tt]])
    #print(tmpVal_Rhat)
    tmpVal_bulk=ess_bulk(chainyE2_NoPerm[,,colnames(chainyE2_NoPerm[,1,])==newnamesList[tt]])
    #print(tmpVal_bulk)
    tmpVal_tail=ess_tail(chainyE2_NoPerm[,,colnames(chainyE2_NoPerm[,1,])==newnamesList[tt]])
    #print(tmpVal_tail)
    Rhatvec=c(Rhatvec,tmpVal_Rhat)
    Bulkvec=c(Bulkvec,tmpVal_bulk)
    Tailvec=c(Tailvec,tmpVal_tail)
  }
}
print(max(Rhatvec))  # Check this is <1.05
print(min(Bulkvec))  # Check this is >100
print(min(Tailvec))  # Check this is >100
##############################################################################################################



########################################################################################
# G: another sensitivity a. of MAIN MODEL  on DELAYS delay=tau +tau shifts the new admissions day forward
########################################################################################
deltaNdays=0;
data <- list( D_NDays=baseNdays+deltaNdays,
              D_NAgeGroups=8,
              D_cumulHospByAge=CumInfSevere_frS,
              D_ageDisIn=ageDis_frS, 
              D_meanAges=ages,
              D_removedsByAge_frac=propRecover_frS,
              D_contactMPre=cmIn,
              D_contactMPost=cmIn,
              D_runType=1,
              D_negTau=3,
              D_posTau=0,
              D_meancontact=mean(rowSums(cmIn))
)
initB_f1 <- function() list(ageExp=array(0.001,dim = 1),coefExp_inv=array(2,dim = 1),deviationFromExp=array(0,dim = 8),transMC = 0.999, transSC = 0.999, transMU = 0.999, transSU = 0.999, phi_inv = 1.001) #
print('G')
fitE_extra3 <- stan(file = "ageDepNB_vFinal.stan",
                     data=data,
                     init=initB_f1,
                     iter=longRunT, 
                     chains=4,
                     warmup = warmUpT,
                     control = list(max_treedepth = 15,adapt_delta=adaptVal)  
) 
summary(fitE_extra3, pars = "ageExp[1]")$summary[c(4,8)]
summary(fitE_extra3, pars = "coefExp[1]")$summary[c(4,8)]
smE3=summary(fitE_extra3)
smmE3=smE3$summary
chainyE3_NoPerm<- extract(fitE_extra3, permuted = FALSE,inc_warmup=FALSE)
chainyE3_YesPerm<- extract(fitE_extra3, permuted = TRUE)
save(chainyE3_NoPerm, file = "vFinal_stanE3_NP.RData")
save(chainyE3_YesPerm, file = "vFinal_stanE3_YP.RData")
save(smmE3, file = "vFinal_stanE3Summary.RData")
################# DIAGNOSTICS #######################################################
namesList=colnames(head(chainyE3_NoPerm[1,,]))
tmpny=which(is.na(smmE3[,10])) 
invalidnamesList=namesList[tmpny]   # CHECK these NA's from the summary - were they uniform? Then discard uniforms for diagnostics...
# Add to newnamesList: contacts_vec[1..8] SevereInc[1..8] both of which are uniform in the sense of not involving fitted parameters...
invalidnamesList=c(invalidnamesList,"contacts_vec[1]","contacts_vec[2]","contacts_vec[3]","contacts_vec[4]","contacts_vec[5]","contacts_vec[6]","contacts_vec[7]","contacts_vec[8]")
invalidnamesList=c(invalidnamesList,"SevereInc[1]","SevereInc[2]","SevereInc[3]","SevereInc[4]","SevereInc[5]","SevereInc[6]","SevereInc[7]","SevereInc[8]")
invalidnamesList=c(invalidnamesList,"lp__")
Rhatvec=numeric(0)
Bulkvec=numeric(0)
Tailvec=numeric(0)
newnamesList=colnames(chainyE3_NoPerm[,1,])
numNames=length(colnames(chainyE3_NoPerm[,1,]))
for (tt in 1:(length(numNames))) {
  if (sum(newnamesList[tt]==invalidnamesList)==0) {
    #print(newnamesList[tt])
    tmpVal_Rhat=Rhat(chainyE3_NoPerm[,,colnames(chainyE3_NoPerm[,1,])==newnamesList[tt]])
    #print(tmpVal_Rhat)
    tmpVal_bulk=ess_bulk(chainyE3_NoPerm[,,colnames(chainyE3_NoPerm[,1,])==newnamesList[tt]])
    #print(tmpVal_bulk)
    tmpVal_tail=ess_tail(chainyE3_NoPerm[,,colnames(chainyE3_NoPerm[,1,])==newnamesList[tt]])
    #print(tmpVal_tail)
    Rhatvec=c(Rhatvec,tmpVal_Rhat)
    Bulkvec=c(Bulkvec,tmpVal_bulk)
    Tailvec=c(Tailvec,tmpVal_tail)
  }
}
print(max(Rhatvec))  # Check this is <1.05
print(min(Bulkvec))  # Check this is >100
print(min(Tailvec))  # Check this is >100
##############################################################################################################






########################################################################################
# H: another sensitivity a. of MAIN MODEL  on DELAYS delay=-tau, -tau shifts the current hospitalisastions days forwards
########################################################################################
deltaNdays=0;
data <- list( D_NDays=baseNdays+deltaNdays,
              D_NAgeGroups=8,
              D_cumulHospByAge=CumInfSevere_frS,
              D_ageDisIn=ageDis_frS, 
              D_meanAges=ages,
              D_removedsByAge_frac=propRecover_frS,
              D_contactMPre=cmIn,
              D_contactMPost=cmIn,
              D_runType=1,
              D_negTau=6, # D_negTau>0 reflects the added delay of mean hospitalisation time
              D_posTau=0,
              D_meancontact=mean(rowSums(cmIn))
)
initB_f1 <- function() list(ageExp=array(0.001,dim = 1),coefExp_inv=array(2,dim = 1),deviationFromExp=array(0,dim = 8),transMC = 0.999, transSC = 0.999, transMU = 0.999, transSU = 0.999, phi_inv = 1.001) #
print('H')
fitE_extra4 <- stan(file = "ageDepNB_vFinal.stan",
                     data=data,
                     init=initB_f1,
                     iter=longRunT, 
                     chains=4,
                     warmup = warmUpT,
                     control = list(max_treedepth = 15,adapt_delta=adaptVal) 
) 
smE4=summary(fitE_extra4)
smmE4=smE4$summary
chainyE4_NoPerm<- extract(fitE_extra4, permuted = FALSE,inc_warmup=FALSE)
chainyE4_YesPerm<- extract(fitE_extra4, permuted = TRUE)
save(chainyE4_NoPerm, file = "vFinal_stanE4_NP.RData")
save(chainyE4_YesPerm, file = "vFinal_stanE4_YP.RData")
save(smmE4, file = "vFinal_stanE4Summary.RData")

################# DIAGNOSTICS #######################################################
namesList=colnames(head(chainyE4_NoPerm[1,,]))
tmpny=which(is.na(smmE4[,10])) 
invalidnamesList=namesList[tmpny]   # CHECK these NA's from the summary - were they uniform? Then discard uniforms for diagnostics...
# Add to newnamesList: contacts_vec[1..8] SevereInc[1..8] both of which are uniform in the sense of not involving fitted parameters...
invalidnamesList=c(invalidnamesList,"contacts_vec[1]","contacts_vec[2]","contacts_vec[3]","contacts_vec[4]","contacts_vec[5]","contacts_vec[6]","contacts_vec[7]","contacts_vec[8]")
invalidnamesList=c(invalidnamesList,"SevereInc[1]","SevereInc[2]","SevereInc[3]","SevereInc[4]","SevereInc[5]","SevereInc[6]","SevereInc[7]","SevereInc[8]")
invalidnamesList=c(invalidnamesList,"lp__")
Rhatvec=numeric(0)
Bulkvec=numeric(0)
Tailvec=numeric(0)
newnamesList=colnames(chainyE4_NoPerm[,1,])
numNames=length(colnames(chainyE4_NoPerm[,1,]))
for (tt in 1:(length(numNames))) {
  if (sum(newnamesList[tt]==invalidnamesList)==0) {
    #print(newnamesList[tt])
    tmpVal_Rhat=Rhat(chainyE4_NoPerm[,,colnames(chainyE4_NoPerm[,1,])==newnamesList[tt]])
    #print(tmpVal_Rhat)
    tmpVal_bulk=ess_bulk(chainyE4_NoPerm[,,colnames(chainyE4_NoPerm[,1,])==newnamesList[tt]])
    #print(tmpVal_bulk)
    tmpVal_tail=ess_tail(chainyE4_NoPerm[,,colnames(chainyE4_NoPerm[,1,])==newnamesList[tt]])
    #print(tmpVal_tail)
    Rhatvec=c(Rhatvec,tmpVal_Rhat)
    Bulkvec=c(Bulkvec,tmpVal_bulk)
    Tailvec=c(Tailvec,tmpVal_tail)
  }
}
print(max(Rhatvec))  # Check this is <1.05
print(min(Bulkvec))  # Check this is >100
print(min(Tailvec))  # Check this is >100
##############################################################################################################

sink()




