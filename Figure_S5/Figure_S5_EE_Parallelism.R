##Load Required Libraries

library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)


setwd("<Path to GitHub Directory")

Good_et_al_LTEE_2000Gen<- read.table("Figure_S5/All_LTEE_Mut_2000_gens",sep="\t",header=TRUE)

Good_et_al_LTEE_2000Gen

Good_et_al_LTEE_2000Gen_melt<-reshape2::melt(Good_et_al_LTEE_2000Gen, id=c("Line", "sample","position","locus","MutType","PositionHit_No_Ans","Parallel"))

Not_parallel<-Good_et_al_LTEE_2000Gen_melt[which(Good_et_al_LTEE_2000Gen_melt$Parallel=="no"),]
Parallel<-Good_et_al_LTEE_2000Gen_melt[which(Good_et_al_LTEE_2000Gen_melt$Parallel=="yes"),]

Figure_S5_Good_et_al_LTEE_2000Gen_Plot<-ggplot()+
  geom_point(data=Not_parallel, aes(x=variable, y=value), color="#888888")+
  geom_line(data=Not_parallel, aes(x=variable, y=value, group=position),color="#888888")+
  geom_point(data=Parallel, aes(x=variable, y=value), color="#aa0000")+
  geom_line(data=Parallel, aes(x=variable, y=value, group=position),color="#aa0000")+
  xlab("Number of Generations")+
  ylab("Relative Frequency")+
  theme_bw()+
  theme(axis.title.y = element_text(face="bold"),axis.title.x = element_text(face="bold"), panel.grid.major.x = element_blank(),panel.grid.minor=element_blank())+
  facet_wrap(~sample, nrow=2)



Patton_et_al<- read.table("Figure_S5/All_Patton_Mut",sep="\t",header=TRUE)

Patton_et_al$Parallelism<-factor(Patton_et_al$Parallelism, levels=c("Among_Treatments","Within_Treatment","None"))


Figure_S5_Patton_et_al_Plot<-ggplot(Patton_et_al, aes(x=Population,fill=Parallelism))+
  geom_bar(stat="count", color="black")+
  ylab("Mutation Count")+
  theme_bw()+
  scale_fill_manual(values=c("#0000bb","#aa0000","#aaaaaa"))+
  theme(legend.background = element_rect(color="black"),legend.position = "inside",legend.position.inside = c(0.1,0.68),axis.title.y = element_text(face="bold"),axis.title.x = element_text(face="bold"),axis.text.x = element_text(angle=45,hjust=1), panel.grid.major.x = element_blank(),panel.grid.minor=element_blank())
  
plot_grid(Figure_S5_Good_et_al_LTEE_2000Gen_Plot,Figure_S5_Patton_et_al_Plot, nrow=2,labels = c("A","B"))



