##Call necessary libraries
library(ggplot2)
library(cowplot)
library(reshape2)
library(ggpubr)
library(dplyr)
library(scales)
library(ggsignif)
library(data.table)
library(Hmisc)


###set working directory to Github Folder
#setwd("<Path to GitHub Folder>")

Biofilm_df<-read.table("Figure_2/Biofilm_DataSet.txt",sep="\t",header=TRUE)

Biofilm_df["Treatment"][Biofilm_df["Treatment"] == "WT_Ans"] <- "Ancestor, WT"
Biofilm_df["Treatment"][Biofilm_df["Treatment"] == "WCT"] <- "Tube, WT"
Biofilm_df["Treatment"][Biofilm_df["Treatment"] == "WFL"] <- "Flask, WT"
Biofilm_df["Treatment"][Biofilm_df["Treatment"] == "FimA_Ans"] <- "Ancestor, FimA-"
Biofilm_df["Treatment"][Biofilm_df["Treatment"] == "FCT"] <- "Tube, FimA-"
Biofilm_df["Treatment"][Biofilm_df["Treatment"] == "FFL"] <- "Flask, FimA-"

Biofilm_df$Treatment <- factor(Biofilm_df$Treatment, levels = c("Ancestor, FimA-", "Tube, FimA-", "Flask, FimA-", "Ancestor, WT", "Tube, WT", "Flask, WT"))


##Figure 2A

BF_Ancestors_Only<-Biofilm_df[which(Biofilm_df$Treatment=="Ancestor, WT"|Biofilm_df$Treatment=="Ancestor, FimA-"),]
BF_Ancestors_Only

###Anova for Ancestors, linear model OD550~Temperature*Treatment
Ancestor.aov <-aov(FinalOD550~Temp*Treatment,data=BF_Ancestors_Only)
summary(Ancestor.aov) 
TukeyHSD(Ancestor.aov)

BF_Ancestor_Plot_Fig_2A<-ggplot(data=BF_Ancestors_Only, aes(x=Temp, y=FinalOD550, col=Background))+
  #geom_jitter(color="#E3E4E8")+
  geom_jitter(alpha=0.3)+
  scale_color_manual(values=c("#2A6DD2","#99282E"))+
  stat_summary(fun=mean, geom="point", pch=20,size=4)+
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar",size=1)+
  ylab("Ancestral Biofilm (Abs: 550 nm)")+
  xlab("Temperature (°C)")+
  scale_y_continuous(limits=c(0,1.5))+
  scale_x_continuous(breaks=c(25,37))+
  #geom_hline(yintercept=0.132, linetype="dashed", color="#6AA0EF")+ 
  #geom_hline(yintercept=0.612, linetype="dashed", color="#c8686d")+
  theme_bw()+
  theme(text = element_text(size=10),axis.text.y = element_text(size=10),panel.grid = element_blank(),axis.text.x = element_text(size=10),axis.title = element_text(size=8,face="bold"), legend.position = "none", legend.title = element_text(size=10,face="bold"), legend.background = element_rect(color="black"))+
  facet_wrap(.~Background)

BF_Ancestor_Plot_Fig_2A

###Figure 2B
Biofilm_Ancestor_and_Final<-Biofilm_df[which(Biofilm_df$Time!=7),]

Biofilm_Ancestor_and_Final

Biofilm_df_37<-Biofilm_Ancestor_and_Final[which(Biofilm_Ancestor_and_Final$Temp==37),]

BF_Evo_37C_Plot_Fig_2B<-ggplot(data=Biofilm_df_37, aes(x=Treatment, y=FinalOD550, col=Treatment))+
  #geom_jitter(color="#E3E4E8")+
  geom_jitter(alpha=0.3)+
  scale_color_manual(values=c("#2A6DD2","#6AA0EF","#a7c9df","#99282E","#c8686d", "#eda4a7"))+
  stat_summary(fun=mean, geom="point", pch=20,size=4)+
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar",size=1)+
  ylab("Evolved Biofilm (Abs: 550 nm)")+
  xlab("Treatment")+
  scale_y_continuous(limits=c(0,0.4))+
  #geom_hline(yintercept=0.132, linetype="dashed", color="#6AA0EF")+ 
  #geom_hline(yintercept=0.612, linetype="dashed", color="#c8686d")+
  theme_bw()+
  theme(text = element_text(size=10),axis.text.y = element_text(size=10),panel.grid = element_blank(),axis.text.x = element_text(angle=45,hjust=1,size=10),axis.title = element_text(size=8,face="bold"), legend.position = "none", legend.title = element_text(size=10,face="bold"), legend.background = element_rect(color="black"))

BF_Evo_37C_Plot_Fig_2B

###Anova for 37C, linear model OD550~Treatment
BF_37.aov <-aov(FinalOD550~Treatment,data=Biofilm_df_37)
summary(BF_37.aov) 
TukeyHSD(BF_37.aov)

###Figure S1

Biofilm_df_25<-Biofilm_Ancestor_and_Final[which(Biofilm_Ancestor_and_Final$Temp==25),]

BF_Evo_25C_Plot_Fig_S1<-ggplot(data=Biofilm_df_25, aes(x=Treatment, y=FinalOD550, col=Treatment))+
  #geom_jitter(color="#E3E4E8")+
  geom_jitter(alpha=0.3)+
  scale_color_manual(values=c("#2A6DD2","#6AA0EF","#a7c9df","#99282E","#c8686d", "#eda4a7"))+
  stat_summary(fun=mean, geom="point", pch=20,size=4)+
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar",size=1)+
  ylab("Evolved Biofilm (Abs: 550 nm)")+
  xlab("Treatment")+
  scale_y_continuous(limits=c(0,1.5))+
  #geom_hline(yintercept=0.132, linetype="dashed", color="#6AA0EF")+ 
  #geom_hline(yintercept=0.612, linetype="dashed", color="#c8686d")+
  theme_bw()+
  facet_wrap(.~Temp)+
  theme(text = element_text(size=10),axis.text.y = element_text(size=10),panel.grid = element_blank(),axis.text.x = element_text(angle=45,hjust=1,size=10),axis.title = element_text(size=8,face="bold"), legend.position = "right", legend.title = element_text(size=10,face="bold"), legend.background = element_rect(color="black"))

BF_Evo_25C_Plot_Fig_S1


##Figure 2C

Biofilm_WT_Tube<-Biofilm_df[which(Biofilm_df$Treatment=="Ancestor, WT"|Biofilm_df$Treatment=="Tube, WT"),]
Biofilm_WT_Tube_37<-Biofilm_WT_Tube[which(Biofilm_WT_Tube$Temp==37),]
Biofilm_WT_Tube_37$Time<-as.character(Biofilm_WT_Tube_37$Time)
  
BF_WT_Tube_37_Plot_Fig_2C<-ggplot(data=Biofilm_WT_Tube_37, aes(x=Time, y=FinalOD550))+
  #geom_jitter(color="#E3E4E8")+
  geom_jitter(alpha=0.3)+
  #scale_color_manual(values=c("#2A6DD2","#6AA0EF","#a7c9df","#99282E","#c8686d", "#eda4a7"))+
  stat_summary(fun=mean, geom="point", pch=20,size=4, color="#99282E")+
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar",size=1, color="#99282E")+
  ylab("Evolved Biofilm (Abs: 550 nm)")+
  xlab("Time")+
  scale_y_continuous(limits=c(0,0.5))+
  #scale_x_continuous(breaks=c(0,7,90))+
  theme_bw()+
  theme(text = element_text(size=10),axis.text.y = element_text(size=10),panel.grid = element_blank(),axis.text.x = element_text(angle=45,hjust=1,size=10),axis.title = element_text(size=8,face="bold"), legend.position = "none", legend.title = element_text(size=10,face="bold"), legend.background = element_rect(color="black"))

BF_WT_Tube_37_Plot_Fig_2C

##Read in Culture Density Data for Figure 2D and S2
CultureOD<-read.table("Figure_2/Culture_Density.txt", sep="\t", header=TRUE)

CultureOD.melt<-reshape2::melt(CultureOD, id=c("Time","Temp","Sample"))
names(CultureOD.melt)<-c("Time","Temperature","Sample","Treatment","OD600")
CultureOD.melt$Treatment <- factor(CultureOD.melt$Treatment, levels = c("FCT","FFL","WCT","WFL"))
CultureOD.melt$Temperature <- factor(CultureOD.melt$Temperature, levels = c("25","37"))
CultureOD.melt$Generations<-CultureOD.melt$Time*3.3

#Figure 2D

CultureOD_37<-CultureOD.melt[which(CultureOD.melt$Temperature==37),]

Culture_Density_37_Plot_Fig_2D<-ggplot(data=CultureOD_37, aes(x=Generations,y=OD600,color=Treatment))+
  geom_point()+
  stat_summary(geom="line",fun.y=median)+
  scale_y_continuous(limits=c(0.5,1.6))+
  ylab("Culture Density (Abs: 600nm)")+
  xlab("Generations")+
  scale_color_manual(values=c("#6AA0EF","#a7c9df","#c8686d","#db9a9d"),labels=c( "Tube, FimA-","Flask, FimA-","Tube, WT","Flask, WT"))+
  theme_bw()+
  theme(legend.position="none",legend.background = element_rect(fill = "white", color = "black"),axis.text.x = element_text(size=10),text = element_text(size=10),axis.text.y = element_text(size=10), axis.title=element_text(size=8,face="bold"))

Culture_Density_37_Plot_Fig_2D

#Plot Figure 2
plot_grid(BF_Ancestor_Plot_Fig_2A,BF_Evo_37C_Plot_Fig_2B,BF_WT_Tube_37_Plot_Fig_2C,Culture_Density_37_Plot_Fig_2D,nrow=1,rel_widths = c(0.75,1,0.5,1))


###Figure S2
CultureOD_25<-CultureOD.melt[which(CultureOD.melt$Temperature==25),]

Culture_Density_25_Plot_Fig_S2<-ggplot(data=CultureOD_25, aes(x=Generations,y=OD600,color=Treatment))+
  geom_point()+
  stat_summary(geom="line",fun.y=median)+
  scale_y_continuous(limits=c(0.5,1.6))+
  ylab("Culture Density (Abs: 600nm)")+
  xlab("Generations")+
  scale_color_manual(values=c("#6AA0EF","#a7c9df","#c8686d","#db9a9d"),labels=c( "Tube, FimA-","Flask, FimA-","Tube, WT","Flask, WT"))+
  theme_bw()+
  facet_wrap(~Temperature)+
  theme(legend.position="right",legend.background = element_rect(fill = "white", color = "black"),axis.text.x = element_text(size=10),text = element_text(size=10),axis.text.y = element_text(size=10), axis.title=element_text(size=8,face="bold"))

Culture_Density_25_Plot_Fig_S2
