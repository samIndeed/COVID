rm(list=setdiff(ls(), "filepath"))

library(ggplot2)
library(Rmisc)
library(grid)
library(gridExtra)
library(latex2exp)
# mean age for france in each age cohort (calculated elsewhere i.e. in main script Estimation_vFinal.r)
ages=c(9.752804,24.495934,34.583572,44.634635,54.469143,64.418965,73.794358,87.23397)

################ INITIAL A. ############################
################ INITIAL A. ############################
################ INITIAL A. ############################
################ INITIAL A. ############################
################ INITIAL A. ############################
################ INITIAL A. ############################
print(load(file = "vFinal_Data/vFinal_stanA_NP.RData"))
print(load(file = "vFinal_Data/vFinal_stanASummary.RData"))
chainy=chainyA_NoPerm
basel=0.018*exp(0.044*ages)
smMean=smmA[16:23-1,1]
sm25=smmA[16:23-1,4]
sm50=smmA[16:23-1,6]
sm975=smmA[16:23-1,8]
###############
df=data.frame(ages,basel)
p<- ggplot(df, aes(x=ages, y=basel))+
  geom_line(size=0.5, aes(colour = c("Doubling every 16 years")),linetype='dashed',show.legend = TRUE) +
  geom_errorbar(aes(ymin=sm25, ymax=sm975), width=1.5,
                position=position_dodge(.9))+
  geom_point(aes(x=ages, y=sm50), width=1.5,
             position=position_dodge(.9))+
  xlab("Age group")+
  ylab("Prob(Severe disease | infection)")+
  scale_colour_manual("Legend", values = c("green"))+
  scale_x_continuous(breaks = ages,labels=seq(1,8))+
  #scale_x_continuous(breaks = round(ages))+
  theme_bw() +
  theme(legend.justification=c(1,0),
        legend.position=c(0.55,0.875),
        legend.title = element_blank(),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        panel.background = element_rect(fill = "transparent",colour =NA), # or theme_blank()
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.text=element_text(size=10),
        legend.key=element_rect(fill = NA,colour = NA, size = 1),
        legend.background = element_rect(fill = "transparent",colour = NA),
        legend.key.size = unit(1.5, "cm"))
p+scale_y_log10()
ggsave("FigureS4.pdf", width = 20, height = 15, units = "cm")
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################





################ DEVIATION A. ##############################
################ DEVIATION A. ##############################
################ DEVIATION A. ##############################
################ DEVIATION A. ##############################
################ DEVIATION A. ##############################
################ DEVIATION A. ##############################
################ DEVIATION A. ##############################

######################### ORIGINAL DEVIATION EXP ONLY #################################
print(load(file = "vFinal_Data/vFinal_stanB_NP.RData"))
print(load(file = "vFinal_Data/vFinal_stanBSummary.RData"))
chainy=chainyB_NoPerm
additProtectExpOnly=rbind(smmB[rownames(smmB)=="additProtect[1]",][c(4,8)],smmB[rownames(smmB)=="additProtect[2]",][c(4,8)],smmB[rownames(smmB)=="additProtect[3]",][c(4,8)],smmB[rownames(smmB)=="additProtect[4]",][c(4,8)],smmB[rownames(smmB)=="additProtect[5]",][c(4,8)],smmB[rownames(smmB)=="additProtect[6]",][c(4,8)],smmB[rownames(smmB)=="additProtect[7]",][c(4,8)],smmB[rownames(smmB)=="additProtect[8]",][c(4,8)])
###

devTit2=character(0)
for (kk in 1:8) {
  devTit2=c(devTit2,"transProt")
}
devType=devTit2
cAge=ages
devL=additProtectExpOnly[,1]*100
devU=additProtectExpOnly[,2]*100
devFromTrans_MS=data.frame(devType,cAge,devL,devU)

boolColors <- as.character(c("transProt" = "#000000"))

pd <- position_dodge(5.5) # move them .05 to the left and right
p4=ggplot(devFromTrans_MS, aes(x=cAge, y=devL, colour=devType, group=devType)) + 
  geom_errorbar(aes(ymin=devL, ymax=devU, colour=devType), width=2, position=pd) +
  xlab("Cohort age (yrs)") +
  ylab("Age-specific deviation (%)") +
  scale_colour_manual(name="", 
                      labels=c("Deviation from exponential P(severe|infection)"),
                      values=boolColors) +                    # Use darker colors, lightness=40
  geom_hline(yintercept=0)+
  annotate("text", x = 20, y = 25, label = "line of no deviation", fontface = 'bold', size = 2.5)+
  annotate("text", x = 10, y = 85, label = "A", size = 5)+
  expand_limits(y = c(-100, 100))+
  scale_y_continuous(breaks=seq(-100,100,20)) +         # Set tick every 4
  theme_bw() +
  theme(legend.justification=c(1,0),
        legend.position=c(1,0),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        legend.text=element_text(size=7),
        legend.key=element_rect(fill = NA,colour = NA, size = 0.25),
        legend.background = element_rect(fill = "transparent",colour = NA),
        legend.key.size = unit(0.4, "cm")
  )               # Position legend in bottom right
######################### DEVIATION FROM INFECTION GIVEN EXPOSURE ##############
print(load(file = "vFinal_Data/vFinal2_stanD2_NP.RData"))
print(load(file = "vFinal_Data/vFinal2_stanD2Summary.RData"))
chainy=chainyD2_NoPerm
namesList=colnames(head(chainy[1,,]))
additProtectExpOnly2=rbind(smmD2[rownames(smmD2)=="additProtect[1]",][c(4,8)],smmD2[rownames(smmD2)=="additProtect[2]",][c(4,8)],smmD2[rownames(smmD2)=="additProtect[3]",][c(4,8)],smmD2[rownames(smmD2)=="additProtect[4]",][c(4,8)],smmD2[rownames(smmD2)=="additProtect[5]",][c(4,8)],smmD2[rownames(smmD2)=="additProtect[6]",][c(4,8)],smmD2[rownames(smmD2)=="additProtect[7]",][c(4,8)],smmD2[rownames(smmD2)=="additProtect[8]",][c(4,8)])
deviationFromUE=rbind(smmD2[rownames(smmD2)=="deviationFromUE[1]",][c(4,8)],smmD2[rownames(smmD2)=="deviationFromUE[2]",][c(4,8)],smmD2[rownames(smmD2)=="deviationFromUE[3]",][c(4,8)],smmD2[rownames(smmD2)=="deviationFromUE[4]",][c(4,8)],smmD2[rownames(smmD2)=="deviationFromUE[5]",][c(4,8)],smmD2[rownames(smmD2)=="deviationFromUE[6]",][c(4,8)],smmD2[rownames(smmD2)=="deviationFromUE[7]",][c(4,8)],smmD2[rownames(smmD2)=="deviationFromUE[8]",][c(4,8)])
###
devTit1=character(0)
devTit2=character(0)
for (kk in 1:8) {
  devTit1=c(devTit1,"transDev")
  devTit2=c(devTit2,"transProt")
}
devType=c(devTit1,devTit2)
cAge=c(ages,ages)
devL=c(deviationFromUE[,1],additProtectExpOnly2[,1])*100
devU=c(deviationFromUE[,2],additProtectExpOnly2[,2])*100
devFromTrans_MS=data.frame(devType,cAge,devL,devU)

boolColors <- as.character(c("transProt" = "#0033FF","transDev"="#000000"))

pd <- position_dodge(4.5) # move them .05 to the left and right
p3=ggplot(devFromTrans_MS, aes(x=cAge, y=devL, colour=devType, group=devType)) + 
  geom_errorbar(aes(ymin=devL, ymax=devU, colour=devType), width=3.5, position=pd) +
  xlab("Cohort age (yrs)") +
  ylab("Age-specific deviation (%)") +
  scale_colour_manual(name="", 
                      labels=c("Deviation from uniform P(infection|exposure)", "Deviation from exponential P(severe|infection)"),
                      values=boolColors) +                    # Use darker colors, lightness=40
  geom_hline(yintercept=0)+
  annotate("text", x = 10, y = 85, label = "B",size = 5)+
  expand_limits(y = c(-100, 100))+
  scale_y_continuous(breaks=seq(-100,100,20)) +         # Set tick every 4
  theme_bw() +
  theme(legend.justification=c(1,0),
        legend.position=c(1,0),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        legend.text=element_text(size=7),
        legend.key=element_rect(fill = NA,colour = NA, size = 0.25),
        legend.background = element_rect(fill = "transparent",colour = NA),
        legend.key.size = unit(0.4, "cm"))               # Position legend in bottom right

######################### DEVIATION FROM TRANSMISSION MILD ONLY ##############
print(load(file = "vFinal_Data/vFinal3_stanD3_NP.RData"))
print(load(file = "vFinal_Data/vFinal3_stanD3Summary.RData"))
chainy=chainyD3_NoPerm
additProtectExpOnly3=rbind(smmD3[rownames(smmD3)=="additProtect[1]",][c(4,8)],smmD3[rownames(smmD3)=="additProtect[2]",][c(4,8)],smmD3[rownames(smmD3)=="additProtect[3]",][c(4,8)],smmD3[rownames(smmD3)=="additProtect[4]",][c(4,8)],smmD3[rownames(smmD3)=="additProtect[5]",][c(4,8)],smmD3[rownames(smmD3)=="additProtect[6]",][c(4,8)],smmD3[rownames(smmD3)=="additProtect[7]",][c(4,8)],smmD3[rownames(smmD3)=="additProtect[8]",][c(4,8)])
deviationFromUT3=rbind(smmD3[rownames(smmD3)=="deviationFromUT[1]",][c(4,8)],smmD3[rownames(smmD3)=="deviationFromUT[2]",][c(4,8)],smmD3[rownames(smmD3)=="deviationFromUT[3]",][c(4,8)],smmD3[rownames(smmD3)=="deviationFromUT[4]",][c(4,8)],smmD3[rownames(smmD3)=="deviationFromUT[5]",][c(4,8)],smmD3[rownames(smmD3)=="deviationFromUT[6]",][c(4,8)],smmD3[rownames(smmD3)=="deviationFromUT[7]",][c(4,8)],smmD3[rownames(smmD3)=="deviationFromUT[8]",][c(4,8)])
namesList=colnames(head(chainy[1,,]))

devType=c(devTit1,devTit2)
cAge=c(ages,ages)
devL=c(deviationFromUT3[,1],additProtectExpOnly3[,1])*100
devU=c(deviationFromUT3[,2],additProtectExpOnly3[,2])*100
devFromTrans_MS=data.frame(devType,cAge,devL,devU)

boolColors <- as.character(c("transProt" = "#990000","transDevM"="#000000"))

pd <- position_dodge(4.5) # move them .05 to the left and right
p2=ggplot(devFromTrans_MS, aes(x=cAge, y=devL, colour=devType, group=devType)) + 
  geom_errorbar(aes(ymin=devL, ymax=devU, colour=devType), width=4.5, position=pd) +
  xlab("Cohort age (yrs)") +
  ylab("Age-specific deviation (%)") +
  scale_colour_manual(name="", 
                      labels=c("Deviation from uniform infectiousness (mild)", "Deviation from exponential P(severe|infection)"),
                      values=boolColors) +                    # Use darker colors, lightness=40
  geom_hline(yintercept=0)+
  annotate("text", x = 10, y = 85, label = "C", size = 5)+
  expand_limits(y = c(-100, 100)) +                        # Expand y range
  scale_y_continuous(breaks=seq(-100,100,20)) +         # Set tick every 4
  theme_bw() +
  theme(legend.justification=c(1,0),
        legend.position=c(1,0),
        legend.text=element_text(size=7),
        legend.key=element_rect(fill = NA,colour = NA, size = 0.25),
        legend.background = element_rect(fill = "transparent",colour = NA),
        legend.key.size = unit(0.4, "cm"))               # Position legend in bottom right

######################### DEVIATION FROM TRANSMISSION SEVERE AND MILD ###################
print(load(file = "vFinal_Data/vFinal2_stanD4_NP.RData"))
print(load(file = "vFinal_Data/vFinal2_stanD4Summary.RData"))
chainy4=chainyD4_NoPerm
additProtectExpOnly4=rbind(smmD4[rownames(smmD4)=="additProtect[1]",][c(4,8)],smmD4[rownames(smmD4)=="additProtect[2]",][c(4,8)],smmD4[rownames(smmD4)=="additProtect[3]",][c(4,8)],smmD4[rownames(smmD4)=="additProtect[4]",][c(4,8)],smmD4[rownames(smmD4)=="additProtect[5]",][c(4,8)],smmD4[rownames(smmD4)=="additProtect[6]",][c(4,8)],smmD4[rownames(smmD4)=="additProtect[7]",][c(4,8)],smmD4[rownames(smmD4)=="additProtect[8]",][c(4,8)])
deviationFromUT=rbind(smmD4[rownames(smmD4)=="deviationFromUT[1]",][c(4,8)],smmD4[rownames(smmD4)=="deviationFromUT[2]",][c(4,8)],smmD4[rownames(smmD4)=="deviationFromUT[3]",][c(4,8)],smmD4[rownames(smmD4)=="deviationFromUT[4]",][c(4,8)],smmD4[rownames(smmD4)=="deviationFromUT[5]",][c(4,8)],smmD4[rownames(smmD4)=="deviationFromUT[6]",][c(4,8)],smmD4[rownames(smmD4)=="deviationFromUT[7]",][c(4,8)],smmD4[rownames(smmD4)=="deviationFromUT[8]",][c(4,8)])
deviationFromUT_S=rbind(smmD4[rownames(smmD4)=="deviationFromUT_S[1]",][c(4,8)],smmD4[rownames(smmD4)=="deviationFromUT_S[2]",][c(4,8)],smmD4[rownames(smmD4)=="deviationFromUT_S[3]",][c(4,8)],smmD4[rownames(smmD4)=="deviationFromUT_S[4]",][c(4,8)],smmD4[rownames(smmD4)=="deviationFromUT_S[5]",][c(4,8)],smmD4[rownames(smmD4)=="deviationFromUT_S[6]",][c(4,8)],smmD4[rownames(smmD4)=="deviationFromUT_S[7]",][c(4,8)],smmD4[rownames(smmD4)=="deviationFromUT_S[8]",][c(4,8)])
namesList=colnames(head(chainy[1,,]))

devTit1=character(0)
devTit1b=character(0)
devTit2=character(0)
for (kk in 1:8) {
  devTit1=c(devTit1,"transDevM")
  devTit1b=c(devTit1b,"transDevS")
  devTit2=c(devTit2,"transProt")
}

devType=c(devTit1,devTit1b,devTit2)
cAge=c(ages,ages,ages)
devL=c(deviationFromUT[,1],deviationFromUT_S[,1],additProtectExpOnly4[,1])*100
devU=c(deviationFromUT[,2],deviationFromUT_S[,2],additProtectExpOnly4[,2])*100
devFromTrans_MS=data.frame(devType,cAge,devL,devU)

boolColors <- as.character(c("transProt" = "#990000","transDevS" = "#9933FF","transDevM"="#000000"))

pd <- position_dodge(7.5) # move them .05 to the left and right
p1=ggplot(devFromTrans_MS, aes(x=cAge, y=devL, colour=devType, group=devType)) + 
  geom_errorbar(aes(ymin=devL, ymax=devU, colour=devType), width=5.5, position=pd) +
  xlab("Cohort age (yrs)") +
  ylab("Age-specific deviation (%)") +
  scale_colour_manual(name="", 
                      labels=c("Deviation from uniform infectiousness (mild)","Deviation from uniform infectiousness (severe)", "Deviation from exponential P(severe|infection)"),
                   values=boolColors) +                    # Use darker colors, lightness=40
  geom_hline(yintercept=0)+
  annotate("text", x = 10, y = 85, label = "D", size = 5)+
  expand_limits(y = c(-100, 100)) +                        # Expand y range
  scale_y_continuous(breaks=seq(-100,100,20)) +         # Set tick every 4
  theme_bw() +
  theme(legend.justification=c(1,0),
        legend.position=c(1,0),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        legend.text=element_text(size=7),
        legend.key=element_rect(fill = NA,colour = NA, size = 0.25),
        legend.background = element_rect(fill = "transparent",colour = NA),
        legend.key.size = unit(0.4, "cm"))               # Position legend in bottom right

################### MAKE FIGURE #############################
#multiplot(p1, p2, p3, p4, cols=2)
grid.newpage()

pdf(file="FigureS5.pdf", width = (1/2.54)*20, height = (1/2.54)*15)
pushViewport(viewport(layout=grid.layout(2,2,widths=c(0.5425,0.4575),heights=c(0.4625,0.5375))))
print(p4, vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(p2, vp=viewport(layout.pos.row=2,layout.pos.col=1))
print(p3, vp=viewport(layout.pos.row=1,layout.pos.col=2))
print(p1, vp=viewport(layout.pos.row=2,layout.pos.col=2))
dev.off()
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################










################ SENSITIVITY A. ############################
################ SENSITIVITY A. ############################
################ SENSITIVITY A. ############################
################ SENSITIVITY A. ############################
################ SENSITIVITY A. ############################
################ SENSITIVITY A. ############################
ageExp1=smmB[rownames(smmB)=="ageExp[1]",][c(4,8)]
ageExp2=smmD2[rownames(smmD2)=="ageExp[1]",][c(4,8)]
ageExp3=smmD3[rownames(smmD3)=="ageExp[1]",][c(4,8)]
ageExp4=smmD4[rownames(smmD4)=="ageExp[1]",][c(4,8)]

print(load(file = "vFinal_Data/vFinal4_stanE1Summary.RData"))
plus1=smmE1[rownames(smmE1)=="ageExp[1]",][c(4,8)]

print(load(file = "vFinal_Data/vFinal_stanE2Summary.RData"))
minus1=smmE2[rownames(smmE2)=="ageExp[1]",][c(4,8)]

print(load(file = "vFinal_Data/vFinal_stanE3Summary.RData"))
delay5=smmE3[rownames(smmE3)=="ageExp[1]",][c(4,8)]

print(load(file = "vFinal_Data/vFinal_stanE4Summary.RData"))
delay10=smmE4[rownames(smmE4)=="ageExp[1]",][c(4,8)]


test=c(1,2,3,4,5,6,7,8)
pin=c('one','one','one','one','one','one','one','one')
ageL=c(ageExp1[1],ageExp2[1],ageExp3[1],ageExp4[1],plus1[1],minus1[1],delay5[1],delay10[1])
ageU=c(ageExp1[2],ageExp2[2],ageExp3[2],ageExp4[2],plus1[2],minus1[2],delay5[2],delay10[2])
sensVec=data.frame(pin,test,ageL,ageU)

ytexts=c("baseline (deviation#1)","baseline + deviation#2","baseline + deviation#3","baseline + deviation#4","baseline + 1 days data","baseline - 1 days data","baseline with 5 day fwd delay","baseline with 10 day fwd delay")

  ggplot(sensVec, aes(y=test, x=ageL, colour=pin, group=pin)) + 
  geom_errorbar(aes(xmin=ageL, xmax=ageU, colour=pin), width=0.25)+
  xlab("Posterior 95% credible interval, m (exp. rate)") +
  ylab(" ")+
  expand_limits(x = c(0, 0.1))+
  expand_limits(y = c(1, 7))+
  scale_x_continuous(breaks=seq(0,0.1,0.01))+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8),labels=ytexts)+
  theme_bw() +
  theme(legend.justification=c(1,0),
          legend.position=c(0.6,0.6),
          axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          legend.text=element_text(size=10),
          legend.key=element_rect(fill = NA,colour = NA, size = 0.5),
          legend.background = element_rect(fill = "transparent",colour = NA),
          legend.key.size = unit(0.4, "cm"))+
  theme(axis.text.y = element_text(face="italic", color="#993333", size=10, angle=35),legend.position = "none")
  ggsave("FigureS6.pdf", width = 20, height = 15, units = "cm")
  #############################################################
  #############################################################
  #############################################################
  #############################################################
  #############################################################
  #############################################################
  
    


  
  
  