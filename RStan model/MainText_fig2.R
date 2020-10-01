rm(list=setdiff(ls(), "filepath"))

library(ggpubr)
library(rstan)

###########################################################################################
# CORRECTED DAILY AGE SPECIFIC HOSPITAL ADMISSIONS SALJE et al 2020, SCIENCE
dailyHospCounts_allReg<-read.table(file="dailyHospCounts_allReg.csv", header=TRUE, sep=",")  # Salje's dailyHospCounts.csv
numD=dim(dailyHospCounts_allReg)[1]
dHospData=rbind(dailyHospCounts_allReg[,2],dailyHospCounts_allReg[,3],dailyHospCounts_allReg[,4],dailyHospCounts_allReg[,5],dailyHospCounts_allReg[,6],dailyHospCounts_allReg[,7],dailyHospCounts_allReg[,8],dailyHospCounts_allReg[,9])
CumInfSevere_frS=rbind(cumsum(dHospData[1,]),cumsum(dHospData[2,]),cumsum(dHospData[3,]),cumsum(dHospData[4,]),cumsum(dHospData[5,]),cumsum(dHospData[6,]),cumsum(dHospData[7,]),cumsum(dHospData[8,]))
###########################################################################################


# This code reads in the chains that resulted from r & stan for the various model formulations.
# The data is then use to generate Fig. 2 main text (Palmer et al 2020)

print(load(file = "vFinal_Data/vFinal_stanB_NP.RData"))
print(load(file = "vFinal_Data/vFinal_stanBSummary.RData"))
print(load(file = "vFinal_Data/vFinal_stanB_YP.RData"))

chainy=chainyB_NoPerm

#s1=smmB[rownames(smmB)=="ageExp[1]",][c(1,4,6,8)]
s2=smmB[rownames(smmB)=="coefExp[1]",][c(1,4,6,8)]
s3=smmB[rownames(smmB)=="relTransCvU_S",][c(1,4,6,8)]
s4=smmB[rownames(smmB)=="relTransCvU_M",][c(1,4,6,8)]
s5=smmB[rownames(smmB)=="relTransC",][c(1,4,6,8)]
s6=smmB[rownames(smmB)=="relTransU",][c(1,4,6,8)]

ss1=smmB[rownames(smmB)=="additProtect[1]",][c(1,4,6,8)]
ss2=smmB[rownames(smmB)=="additProtect[2]",][c(1,4,6,8)]
ss3=smmB[rownames(smmB)=="additProtect[3]",][c(1,4,6,8)]
ss4=smmB[rownames(smmB)=="additProtect[4]",][c(1,4,6,8)]
ss5=smmB[rownames(smmB)=="additProtect[5]",][c(1,4,6,8)]
ss6=smmB[rownames(smmB)=="additProtect[6]",][c(1,4,6,8)]
ss7=smmB[rownames(smmB)=="additProtect[7]",][c(1,4,6,8)]
ss8=smmB[rownames(smmB)=="additProtect[8]",][c(1,4,6,8)]

sBig=rbind(s5,s6,s3,s4,s2,ss1)
rownames(sBig)=c("relTransC","relTransU","relTransCvU_S","relTransCvU_M","coefExp[1]","additProtect[1]")
print('Table S1, Palmer et al 2020')
print(sBig)

### EPIDEMIC SIMULATION ########
days=seq(1,24)
simRange2=cbind(rbind(smmB[rownames(smmB)=="epi_pred[1]",][c(1,4,6,8)],
               smmB[rownames(smmB)=="epi_pred[2]",][c(1,4,6,8)],
               smmB[rownames(smmB)=="epi_pred[3]",][c(1,4,6,8)],
               smmB[rownames(smmB)=="epi_pred[4]",][c(1,4,6,8)],
               smmB[rownames(smmB)=="epi_pred[5]",][c(1,4,6,8)],
               smmB[rownames(smmB)=="epi_pred[6]",][c(1,4,6,8)],
               smmB[rownames(smmB)=="epi_pred[7]",][c(1,4,6,8)],
               smmB[rownames(smmB)=="epi_pred[8]",][c(1,4,6,8)],
               smmB[rownames(smmB)=="epi_pred[9]",][c(1,4,6,8)],
               smmB[rownames(smmB)=="epi_pred[10]",][c(1,4,6,8)],
               smmB[rownames(smmB)=="epi_pred[11]",][c(1,4,6,8)],
               smmB[rownames(smmB)=="epi_pred[12]",][c(1,4,6,8)],
               smmB[rownames(smmB)=="epi_pred[13]",][c(1,4,6,8)],
               smmB[rownames(smmB)=="epi_pred[14]",][c(1,4,6,8)],
               smmB[rownames(smmB)=="epi_pred[15]",][c(1,4,6,8)],
               smmB[rownames(smmB)=="epi_pred[16]",][c(1,4,6,8)],
               smmB[rownames(smmB)=="epi_pred[17]",][c(1,4,6,8)],
               smmB[rownames(smmB)=="epi_pred[18]",][c(1,4,6,8)],
               smmB[rownames(smmB)=="epi_pred[19]",][c(1,4,6,8)],
               smmB[rownames(smmB)=="epi_pred[20]",][c(1,4,6,8)],
               smmB[rownames(smmB)=="epi_pred[21]",][c(1,4,6,8)],
               smmB[rownames(smmB)=="epi_pred[22]",][c(1,4,6,8)],
               smmB[rownames(smmB)=="epi_pred[23]",][c(1,4,6,8)],
               smmB[rownames(smmB)=="epi_pred[24]",][c(1,4,6,8)]),days)


################################
val=c(simRange2[,3],simRange2[,2],simRange2[,4])
pcile=c(rep('a',length(simRange2[,1])),rep('b',length(simRange2[,1])),rep('c',length(simRange2[,1])))
dayIs=c(seq(1,24),seq(1,24),seq(1,24))
dataIs=data.frame(val,pcile,dayIs)

nums=c(colSums(CumInfSevere_frS[,2:25]-CumInfSevere_frS[,1:24]),2500)
dummy=rep('3',25)
days=c(seq(1,24),2.5)
data2=data.frame(days,nums,dummy)

seg_df1 <- data.frame(x=c(17),
                     y=c(0),
                     xend=c(17),
                     yend=c(2000))

plot1=ggplot(dataIs, aes(x = dayIs,y = val))+
  geom_line(aes(color = pcile,linetype= pcile), size=1.25)+
  geom_point(data = data2, mapping = aes(x = days, y = nums), color = 'black', size=2)+
  scale_colour_manual(name="",labels=c('Median fit','95% credible interval','95% credible interval'),values = c("blue", "lightblue", "lightblue")) +  
  scale_linetype_manual(name="",labels=c('Median fit','95% credible interval','95% credible interval'),values = c('dashed', 'solid', 'solid')) +  
  xlab("Time, days (from 01/03/2020)") +
  ylab("#New severe cases (all ages)") +
  expand_limits(y = c(0, 4000))+
  scale_y_continuous(breaks=seq(0,4000,1000)) +         # Set tick every 4
  theme_bw() +
  annotate("text", x = c(7,21.5), y = c(2500,80), label = c("France data","National lockdown"),colour=c('black','darkred'), size = c(3.5,3))+
  geom_segment(data=seg_df1, aes(x, y, xend=xend, yend=yend),size=1,col="darkred")+
  theme(legend.justification=c(1,0),
         legend.position=c(0.6,0.6),
         axis.title.x=element_text(size=12),
         axis.title.y=element_text(size=12),
         axis.text.x=element_text(size=12),
         axis.text.y=element_text(size=12),
         #panel.background = element_rect(fill = "transparent",colour =NA), # or theme_blank()
         #panel.grid.minor = element_blank(),
         #panel.grid.major = element_blank(),
         legend.text=element_text(size=10),
         legend.key=element_rect(fill = NA,colour = NA, size = 0.5),
         legend.background = element_rect(fill = "transparent",colour = NA),
         legend.key.size = unit(0.4, "cm"))+
  xlim(0,25)

################################

tmp=chainy[,,colnames(chainy[1,,])=="ageExp[1]"]
s1=smmB[rownames(smmB)=="ageExp[1]",][c(4,8)]
dim(tmp)
tmp2=as.vector(tmp) 
length(tmp2)
tmp3=as.data.frame(tmp2)

seg_df <- data.frame(x=c(s1[1],s1[2]),
                     y=c(-150,-150),
                     xend=c(s1[1],s1[2]),
                     yend=c(150,150))
seg_df2 <- data.frame(x=c(s1[1],s1[2]),
                     y=c(-120,-120),
                     xend=c(s1[1],s1[2]),
                     yend=c(120,120))
op <- par(cex = 0.75)   
p=ggplot(data=tmp3, aes(tmp3$tmp2))
plot2=p+geom_histogram(bins=40,col="black",fill="red3")+
  xlab("Posterior distribution, m (exp. rate)") +
  ylab(" ") +
  expand_limits(y = c(-100, 100))+
  scale_y_continuous(breaks=NULL) +         # Set tick every 4
  theme_bw() +
  theme(legend.justification=c(1,0),
        legend.position=c(1,0),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        #panel.background = element_rect(fill = "transparent",colour =NA), # or theme_blank()
        #panel.grid.minor = element_blank(),
        #panel.grid.major = element_blank(),
        legend.text=element_text(size=10),
        legend.key=element_rect(fill = NA,colour = NA, size = 0.25),
        legend.background = element_rect(fill = "transparent",colour = NA),
        legend.key.size = unit(0.4, "cm"))+
  geom_segment(data=seg_df, aes(x, y, xend=xend, yend=yend),size=2,col="black")+
  geom_segment(data=seg_df2, aes(x, y, xend=xend, yend=yend),size=1.5,col="black")+
  geom_point(aes(x=0.044,y=90),color = "black", size=3.5,shape=23,fill="yellow",stroke = 1)+
  xlim(0.04,0.0565)

################################

additProtectExpOnly=rbind(smmB[rownames(smmB)=="additProtect[1]",][c(4,8)],smmB[rownames(smmB)=="additProtect[2]",][c(4,8)],smmB[rownames(smmB)=="additProtect[3]",][c(4,8)],smmB[rownames(smmB)=="additProtect[4]",][c(4,8)],smmB[rownames(smmB)=="additProtect[5]",][c(4,8)],smmB[rownames(smmB)=="additProtect[6]",][c(4,8)],smmB[rownames(smmB)=="additProtect[7]",][c(4,8)],smmB[rownames(smmB)=="additProtect[8]",][c(4,8)])

devTit1=character(0)
for (kk in 1:8) {
  devTit1=c(devTit1,"transExp")
}
devType=c(devTit1)
cAge=seq(1,8)
devL=additProtectExpOnly[,1]*100
devU=additProtectExpOnly[,2]*100
devFromTrans_MS=data.frame(devType,cAge,devL,devU)
library(ggplot2)
boolColors <- as.character(c("transProt" = "#000000"))
pd <- position_dodge(4.5) # move them .05 to the left and right
plot3=ggplot(devFromTrans_MS, aes(x=cAge, y=devL, colour=devType, group=devType)) + 
  geom_errorbar(aes(ymin=devL, ymax=devU, colour=devType), width=0.5, position=pd) +
  xlab("Age group") +
  ylab("Age-specific deviation (%)") +
  scale_x_continuous(breaks = seq(1,8))+
  scale_colour_manual(name="", 
                      labels=c("Deviation from exponential P(severe|infection)"),
                      values=boolColors) +                    # Use darker colors, lightness=40
  geom_hline(yintercept=0, colour='darkred')+
  annotate("text", x = 7, y = 15, label = "line of no deviation", fontface = 'bold',colour='darkred', size = 3)+
  #annotate("text", x = 10, y = 85, label = "B", size = 9)+
  #ggtitle("Fig. 2b, uniform infectiousness, Prob. infection | exposure") +
  expand_limits(y = c(-100, 100))+
  scale_y_continuous(breaks=seq(-100,100,20)) +         # Set tick every 4
  expand_limits(x = c(0.5, 8.5))+
  #scale_x_continuous(breaks=seq(1,8)) +         # Set tick every 4
  theme_bw() +
  theme(legend.justification=c(1,0),
        legend.position=c(1,0),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        #panel.background = element_rect(fill = "transparent",colour =NA), # or theme_blank()
        #panel.grid.minor = element_blank(),
        #panel.grid.major = element_blank(),
        legend.text=element_text(size=10),
        legend.key=element_rect(fill = NA,colour = NA, size = 0.25),
        legend.background = element_rect(fill = "transparent",colour = NA),
        legend.key.size = unit(0.4, "cm")
  )
#################################################################################
ggarrange(plot1,plot2,plot3,labels = c("A", "B", "C"),ncol = 3, nrow = 1)
#################################################################################
ggsave("Figure2.pdf", width = 32, height = 10, units = "cm")










