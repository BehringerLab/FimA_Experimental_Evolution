#Load Required Libraries not provided by Growth Curve R script.
library(ggplot2)
library(ggpubr)
library(reshape2)
library(cowplot)
library(dplyr)
library(chron)
library(ggrepel)
library(data.table)
library(rstatix)

# Setup Work Environment
rm(list=ls())
#setwd("<Path to Github Folder>")

#create Fitness data tables
Comp_EvolvedCondition <- read.table("Figure_3/Evolved_Condition_Comp_Data.txt", header=T)
Comp_ReciprocalCondition <- read.table("Figure_3/Reciprocal_Condition_Comp_Data_New.txt", header=T)
Comp_ReciprocalCondition
Comp_EvolvedCondition$Day<-as.character(Comp_EvolvedCondition$Day)
Comp_EvolvedCondition$Temp<-as.character(Comp_EvolvedCondition$Temp)

Comp_ReciprocalCondition$Day<-as.character(Comp_ReciprocalCondition$Day)
Comp_ReciprocalCondition$Temp<-as.character(Comp_ReciprocalCondition$Temp)

#Create Mean Fitness Tables 
MeanFit_EvolvedCondition<-as.data.table(Comp_EvolvedCondition%>% group_by(Treatment,Temp) %>% summarise(Evolved=mean(Fitness)))
MeanFit_EvolvedCondition
MeanFit_ReciprocalCondition<-as.data.table(Comp_ReciprocalCondition%>% group_by(Treatment,Temp) %>% summarise(Evolved=mean(Fitness)))

#Combine tables
MeanFit_EvolvedCondition$Environment<-"Evolved"
MeanFit_ReciprocalCondition$Environment<-"Reciprocal"
rbind(MeanFit_EvolvedCondition,MeanFit_ReciprocalCondition)

#Create Evolved Condition Summary Table
Comp_EvolvedCondition %>% group_by(Treatment,Temp) %>% summarise(count.data=n())
EvolvedFitness<-Comp_EvolvedCondition %>% group_by(Treatment,Temp,Evolved) %>% summarise(mean(Fitness))
EvolvedFitness<-as.data.frame(EvolvedFitness)
colnames(EvolvedFitness) <- c('Treatment','Temp','Evolved','Fitness')

#Create Reciprocal Condition Summary Table
ReciprocalFitness<-Comp_ReciprocalCondition%>% group_by(Treatment,Temp,Evolved) %>% summarise(mean(Fitness))
ReciprocalFitness<-as.data.frame(ReciprocalFitness)

Comp_ReciprocalCondition %>% group_by(Temp,Treatment) %>% summarise(median.fit=median(Fitness))

#ANOVA for evolved condition
evolved.aov <-aov(Fitness~Temp*Treatment,data=EvolvedFitness)
summary(evolved.aov) 
TukeyHSD(evolved.aov)


##Create Data table for Ancestors
Comp_Ancestors <- read.table("Figure_3/Ancestor_Competition.txt", header=T)
Comp_Ancestors.Mean<-Comp_Ancestors  %>% summarise(Ancestor.mean = mean(W_fimA),Ancestor.sd = sd(W_fimA),Ancestor.n = n(),Ancestor.se = Ancestor.sd / sqrt(Ancestor.n))

supp.labs <- c("25C", "37C")
names(supp.labs) <- c("25", "37")


###plot Figure 3A
Comp_EvolvedCondition_37C<-Comp_EvolvedCondition[which(Comp_EvolvedCondition$Temp==37),]

Figure_3A_Evolved_Fitness_37<-ggplot(data=Comp_EvolvedCondition_37C, aes(x=Treatment,y=Fitness,col=Treatment))+
  geom_jitter(width=0.15, alpha=0.5)+
  stat_summary(fun.data = "mean_se", geom = "errorbar")+
  geom_hline(data=Comp_Ancestors.Mean,aes(yintercept=Ancestor.mean+Ancestor.sd),linetype="dashed", color="darkgrey")+
  geom_hline(data=Comp_Ancestors.Mean,aes(yintercept=Ancestor.mean-Ancestor.sd),linetype="dashed", color="darkgrey")+
  scale_color_manual(values=c("#6AA0EF","#a7c9df","#c8686d","#db9a9d"),name="Evolution Treatment",labels=c("Tube, FimA-", "Flask, FimA-", "Tube, WT", "Flask, WT"))+
  scale_x_discrete(labels=c("Tube, FimA-", "Flask, FimA-", "Tube, WT", "Flask, WT"))+
  scale_y_continuous(limits=c(0.85,1.55))+
  xlab("Evolution Treatment")+
  theme_bw()+
  theme(legend.position="none",legend.background = element_rect(fill = "white", color = "black"), axis.text.x = element_text(angle=22,hjust=1),axis.title = element_text(face="bold",size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_grid(.~Temp,labeller=labeller(Temp=supp.labs))

Figure_3A_Evolved_Fitness_37

###plot Figure S3A
Comp_EvolvedCondition_25C<-Comp_EvolvedCondition[which(Comp_EvolvedCondition$Temp==25),]

Figure_S3A_Evolved_Fitness_25<-ggplot(data=Comp_EvolvedCondition_25C, aes(x=Treatment,y=Fitness,col=Treatment))+
  geom_jitter(width=0.15, alpha=0.5)+
  stat_summary(fun.data = "mean_se", geom = "errorbar")+
  geom_hline(data=Comp_Ancestors.Mean,aes(yintercept=Ancestor.mean+Ancestor.sd),linetype="dashed", color="darkgrey")+
  geom_hline(data=Comp_Ancestors.Mean,aes(yintercept=Ancestor.mean-Ancestor.sd),linetype="dashed", color="darkgrey")+
  scale_color_manual(values=c("#6AA0EF","#a7c9df","#c8686d","#db9a9d"),name="Evolution Treatment",labels=c("Tube, FimA-", "Flask, FimA-", "Tube, WT", "Flask, WT"))+
  scale_x_discrete(labels=c("Tube, FimA-", "Flask, FimA-", "Tube, WT", "Flask, WT"))+
  scale_y_continuous(limits=c(0.85,1.55))+
  xlab("Evolution Treatment")+
  theme_bw()+
  theme(legend.position="none",legend.background = element_rect(fill = "white", color = "black"), axis.text.x = element_text(angle=22,hjust=1),axis.title = element_text(face="bold",size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_grid(.~Temp,labeller=labeller(Temp=supp.labs))

Figure_S3A_Evolved_Fitness_25


###Import Data for Competitions vs FimA Ancestor
Comp_WT_Evo_Fim_Ans<-read.table("Figure_3/WT_Evo_v_Fim_Ans.txt",sep='\t',header=TRUE)
Comp_WT_Evo_Fim_Ans$Condition <- factor(Comp_WT_Evo_Fim_Ans$Condition , levels = c("Tube","Flask"))

Comp_WT_Evo_Fim_Ans %>% group_by(Condition,Population) %>% summarise(mean.fit=mean(Fitness))
Comp_WT_Evo_Fim_Ans

Comp_EvolvedCondition[which(Comp_EvolvedCondition$Treatment=="FCT" & Comp_EvolvedCondition$Temp==37),]


###Plot Figure 3B

Figure_3B_Comp_vs_FimA_Ans<-ggplot(data=Comp_WT_Evo_Fim_Ans, aes(x=Condition,y=Fitness,col=Condition))+
  geom_jitter(width=0.15, alpha=0.5)+
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar")+
  geom_hline(yintercept=0.825,linetype="dashed", color="#c8686d")+
  geom_hline(yintercept=0.224,linetype="dashed", color="#db9a9d")+
  scale_color_manual(values=c("#c8686d","#db9a9d"),name="Evolution Treatment",labels=c("Tube, WT", "Flask, WT"))+
  scale_x_discrete(labels=c("Tube, WT", "Flask, WT"))+
  #scale_y_continuous(limits=c(0.85,1.55))+
  xlab("Evolution Treatment")+
  ylab("Fitness vs. fimA Ancestor")+
  theme_bw()+
  theme(legend.position="none",legend.background = element_rect(fill = "white", color = "black"), axis.text.x = element_text(angle=22,hjust=1),axis.title = element_text(face="bold",size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_grid(.~Population)

Figure_3B_Comp_vs_FimA_Ans


###Plot top row of Figure 3

Figure_3_Row1<-plot_grid(Figure_3A_Evolved_Fitness_37,Figure_3B_Comp_vs_FimA_Ans, labels = c("A","B"))
Figure_3_Row1

#ANOVA for Figure 3B
WTvsFimA.aov <-aov(Fitness~Population*Condition,data=Comp_WT_Evo_Fim_Ans)
summary(WTvsFimA.aov) 
TukeyHSD(WTvsFimA.aov)


#Reciprocal Conditions
Comp_ReciprocalCondition %>% group_by(Treatment,Temp) %>% summarise(count.data=n())

ggplot(data=Comp_ReciprocalCondition, aes(x=Treatment,y=Fitness,col=Treatment))+
  geom_jitter(width=0.15, alpha=0.5)+
  stat_summary(fun.data = "mean_se", geom = "errorbar")+
  geom_hline(data=Comp_Ancestors.Mean,aes(yintercept=Ancestor.mean+Ancestor.sd),linetype="dashed", color="darkgrey")+
  geom_hline(data=Comp_Ancestors.Mean,aes(yintercept=Ancestor.mean-Ancestor.sd),linetype="dashed", color="darkgrey")+
  scale_color_manual(values=c("#6AA0EF","#a7c9df","#c8686d","#db9a9d"),name="Evolution Treatment",labels=c("Tube, FimA-", "Flask, FimA-", "Tube, WT", "Flask, WT"))+
  #scale_x_discrete(labels=c("Flask, FimA-(Tube)", "Tube, FimA-(Flask)", "Flask, WT(Tube)", "Tube, WT(Flask)"))+
  #scale_y_continuous(limits=c(0.85,1.55))+
  theme_bw()+
  theme(legend.position=c(0.88, 0.75),,legend.background = element_rect(fill = "white", color = "black"),axis.title.x = element_blank(), axis.text.x = element_text(angle=30,hjust=1))+
  facet_wrap(.~Temp,labeller=labeller(Temp=supp.labs),scales="free_y")

Comp_EvolvedCondition
Comp_ReciprocalCondition

myvars <- names(Comp_EvolvedCondition) %in% c("Day", "Evolved","Ans", "Treatment","Temp","EC_RC","Fitness")
Comp_EvolvedCondition.sub <- Comp_EvolvedCondition[myvars]

Fim1<-Comp_EvolvedCondition.sub$Treatment %in% c("FCT","FFL")
Comp_EvolvedCondition.sub$Mutation[Fim1] <- "FimA-"

WT1<-Comp_EvolvedCondition.sub$Treatment %in% c("WFL","WCT")
Comp_EvolvedCondition.sub$Mutation[WT1] <- "WT"

Fim2<-Comp_ReciprocalCondition$Treatment %in% c("FCT","FFL")
Comp_ReciprocalCondition$Mutation[Fim2] <- "FimA-"

WT2<-Comp_ReciprocalCondition$Treatment %in% c("WFL","WCT")
Comp_ReciprocalCondition$Mutation[WT2] <- "WT"

Tube1<-Comp_EvolvedCondition.sub$Treatment %in% c("FCT","WCT")
Comp_EvolvedCondition.sub$Environment[Tube1] <- "Tube"

Flask1<-Comp_EvolvedCondition.sub$Treatment %in% c("FFL","WFL")
Comp_EvolvedCondition.sub$Environment[Flask1] <- "Flask"

Flask2<-Comp_ReciprocalCondition$Treatment %in% c("FCT","WCT")
Comp_ReciprocalCondition$Environment[Flask2] <- "Flask"

Tube2<-Comp_ReciprocalCondition$Treatment %in% c("FFL","WFL")
Comp_ReciprocalCondition$Environment[Tube2] <- "Tube"

Comp_ReciprocalCondition.sub<-select(Comp_ReciprocalCondition,Day,Evolved,Ans,Treatment,Temp,EC_RC,Fitness,Mutation,Environment)



#combine Evolved and Reciprocal data frames
FimA_Fitness<-rbind(Comp_EvolvedCondition.sub,Comp_ReciprocalCondition.sub)

#ANOVA for reciprocal condition
reciprocal.aov <-aov(Fitness~Temp*Treatment,data=Comp_ReciprocalCondition)
summary(reciprocal.aov) 
TukeyHSD(reciprocal.aov)

#ANOVA for evolved and reciprocal conditions
ECRC.aov <-aov(Fitness~Temp*Treatment*EC_RC,data=FimA_Fitness)
summary(ECRC.aov) 
Tukey_ECRC <- TukeyHSD(ECRC.aov)

#Combined plot for evolved and reciprocal conditions
LocalAdaptPlot<-ggplot(data=FimA_Fitness, aes(x=Environment,y=Fitness,col=Treatment,group=Treatment))+
  geom_jitter(width=0.15, alpha=0.5)+
  stat_summary(fun.data = "mean_se", geom = "errorbar", width=0.5)+
  stat_summary(fun.y = "mean", geom = "line")+
  #geom_hline(data=GP_MY_Comp_Ancestors.Mean,aes(yintercept=Ancestor.mean+Ancestor.sd),linetype="dashed", color="darkgrey")+
  #geom_hline(data=GP_MY_Comp_Ancestors.Mean,aes(yintercept=Ancestor.mean-Ancestor.sd),linetype="dashed", color="darkgrey")+
  scale_color_manual(values=c("#6AA0EF","#a7c9df","#c8686d","#eda4a7"),name="Evolution Treatment",labels=c("Tube, FimA-", "Flask, FimA-", "Tube, WT", "Flask, WT"))+
  #scale_x_discrete(labels=c("Flask, FimA-(Tube)", "Tube, FimA-(Flask)", "Flask, WT(Tube)", "Tube, WT(Flask)"))+
  #scale_y_continuous(limits=c(0.85,1.55))+
  theme_bw()+
  theme(legend.position="bottom",legend.title = element_text(face="bold",size=10),legend.background = element_rect(fill = "white", color = "black"),axis.title = element_text(face="bold",size=10))+
  facet_grid(Temp~Mutation,labeller=labeller(Temp=supp.labs))
#,scales="free_y")

## Plot Figure 3C
Local_Adapt_Fitness_37<-FimA_Fitness[which(FimA_Fitness$Temp==37),]
Local_Adapt_Fitness_37

Figure_3C_LocalAdaptPlot_37<-ggplot(data=Local_Adapt_Fitness_37, aes(x=Environment,y=Fitness,col=Treatment,group=Treatment))+
  geom_jitter(width=0.15, alpha=0.5)+
  stat_summary(fun.data = "mean_se", geom = "errorbar", width=0.5)+
  stat_summary(fun.y = "mean", geom = "line")+
  #geom_hline(data=GP_MY_Comp_Ancestors.Mean,aes(yintercept=Ancestor.mean+Ancestor.sd),linetype="dashed", color="darkgrey")+
  #geom_hline(data=GP_MY_Comp_Ancestors.Mean,aes(yintercept=Ancestor.mean-Ancestor.sd),linetype="dashed", color="darkgrey")+
  scale_color_manual(values=c("#6AA0EF","#a7c9df","#c8686d","#eda4a7"),name="Evolution Treatment",labels=c("Tube, FimA-", "Flask, FimA-", "Tube, WT", "Flask, WT"))+
  #scale_x_discrete(labels=c("Flask, FimA-(Tube)", "Tube, FimA-(Flask)", "Flask, WT(Tube)", "Tube, WT(Flask)"))+
  scale_y_continuous(limits=c(0.75,1.8))+
  theme_bw()+
  theme(legend.position="bottom",legend.title = element_text(face="bold",size=10),legend.background = element_rect(fill = "white", color = "black"),axis.title = element_text(face="bold",size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_grid(Temp~Mutation,labeller=labeller(Temp=supp.labs))

Figure_3C_LocalAdaptPlot_37

Figure_3_Row1


plot_grid(Figure_3_Row1,Figure_3C_LocalAdaptPlot_37,nrow=2,rel_heights = c(1,1),labels = c("","C"))

###Figure S3B Local Adaptation 25C

Local_Adapt_Fitness_25<-FimA_Fitness[which(FimA_Fitness$Temp==25),]
Local_Adapt_Fitness_25

Figure_S3B_LocalAdaptPlot_25<-ggplot(data=Local_Adapt_Fitness_25, aes(x=Environment,y=Fitness,col=Treatment,group=Treatment))+
  geom_jitter(width=0.15, alpha=0.5)+
  stat_summary(fun.data = "mean_se", geom = "errorbar", width=0.5)+
  stat_summary(fun.y = "mean", geom = "line")+
  #geom_hline(data=GP_MY_Comp_Ancestors.Mean,aes(yintercept=Ancestor.mean+Ancestor.sd),linetype="dashed", color="darkgrey")+
  #geom_hline(data=GP_MY_Comp_Ancestors.Mean,aes(yintercept=Ancestor.mean-Ancestor.sd),linetype="dashed", color="darkgrey")+
  scale_color_manual(values=c("#6AA0EF","#a7c9df","#c8686d","#eda4a7"),name="Evolution Treatment",labels=c("Tube, FimA-", "Flask, FimA-", "Tube, WT", "Flask, WT"))+
  #scale_x_discrete(labels=c("Flask, FimA-(Tube)", "Tube, FimA-(Flask)", "Flask, WT(Tube)", "Tube, WT(Flask)"))+
  scale_y_continuous(limits=c(0.75,1.8))+
  theme_bw()+
  theme(legend.position="bottom",legend.title = element_text(face="bold",size=10),legend.background = element_rect(fill = "white", color = "black"),axis.title = element_text(face="bold",size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_grid(Temp~Mutation,labeller=labeller(Temp=supp.labs))

Figure_S3B_LocalAdaptPlot_25

plot_grid(Figure_S3A_Evolved_Fitness_25,Figure_S3B_LocalAdaptPlot_25, rel_widths = c(1,2), labels=c("A","B") )


###Interaction Terms

FimA_Interaction<-Local_Adapt_Fitness_37[which(Local_Adapt_Fitness_37$Mutation=="FimA-"),]
WT_Interaction<-Local_Adapt_Fitness_37[which(Local_Adapt_Fitness_37$Mutation=="WT"),]


###P-value FimA- Interaction
Local_Effect_FimA <- lm(Fitness~Treatment*Environment, data = FimA_Interaction)

Local_Effect_FimA
summary(Local_Effect_FimA,confint = TRUE, digits = 3)
partial_eta_squared(car::Anova(Local_Effect_FimA, type = "3"))



###P-value WT Interaction

Local_Effect_WT<- lm(Fitness~Treatment*Environment, data = WT_Interaction)
summary(Local_Effect_WT,confint = TRUE, digits = 3)
partial_eta_squared(car::Anova(Local_Effect_WT, type = 3))
