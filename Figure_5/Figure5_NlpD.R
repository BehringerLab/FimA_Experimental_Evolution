#Load Required Libraries 
library(ggplot2)
library(ggpubr)
library(reshape2)
library(cowplot)
library(dplyr)
library(chron)
library(ggrepel)
library(data.table)
library(rstatix)

#setwd("<Path to GitHub Dir>")


###Plot Competition Data for nlpD mutants
nlpD_Comps<-read.table("Figure_5/NlpD_Comps.txt",sep="\t",header=TRUE)
nlpD_Comps

nlpD_Comps$Competition <- factor(nlpD_Comps$Competition , levels = c("fimA_fimAnlpD","WT_fimE","WT_nlpD","WT_fimEnlpD"))


Figure_5B_NlpD_Comp_Plot<-ggplot(nlpD_Comps,aes(x=Competition,y=W_Strain2,fill=Competition))+
  stat_summary(fun = "mean", geom = "bar", color="black")+
  stat_summary(fun.data = "mean_se", geom = "errorbar", width=0.2, color="black")+
  geom_hline(yintercept = 1, color="black")+
  scale_y_log10(limits=c(0.4,1.6),breaks=c(0.5,0.75,1,1.25,1.5))+
  ylab("Relative Fitness")+
  scale_fill_manual(values=c("#6AA0EF","#c8686d","#c8686d","#c8686d"))+
  theme_bw()+
  theme(legend.position = "none", axis.title.x = element_blank(),axis.title.y = element_text(face="bold"),panel.grid.major.x = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(angle=33, hjust = 1))+
  facet_grid(~Growth.Vessel)

Figure_5B_NlpD_Comp_Plot

nlpD_Comps %>% group_by(Growth.Vessel,Competition) %>% summarise(mean.fit=mean(W_Strain2))


###Plot Biofilm Data for nlpD mutants
nlpD_Biofilm<-read.table("Figure_5/NlpD_Biofilm.txt",sep="\t",header=TRUE)
nlpD_Biofilm

nlpD_Biofilm$Strain <- factor(nlpD_Biofilm$Strain, levels = c("FimA","FimA/NlpD","WT","FimE","NlpD","FimE/NlpD"))


Figure_5C_NlpD_Biofilm_Plot<-ggplot(nlpD_Biofilm,aes(x=Strain,y=OD550,fill=Strain))+
  stat_summary(fun = "mean", geom = "bar", color="black")+
  stat_summary(fun.data = "mean_se", geom = "errorbar", width=0.2, color="black")+
  ylab("Biofilm: Abs (550nm)")+
  scale_fill_manual(values=c("#6AA0EF","#6AA0EF","#c8686d","#c8686d","#c8686d","#c8686d"))+
  theme_bw()+
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_text(face="bold"),panel.grid.major.x = element_blank(),panel.grid.minor = element_blank(), axis.text.x = element_text(angle=33, hjust = 1))

Figure_5C_NlpD_Biofilm_Plot
plot_grid(Figure_5B_NlpD_Comp_Plot,Figure_5C_NlpD_Biofilm_Plot,nrow=1,align = "vh", labels=c("B","C"))

