library(vegan)
library(ggfortify)
library(ggdendro)
library(ggvegan)
library(ggplot2)
library(cowplot)
library(FactoMineR)
library(ggrepel)
library(ggpubr)

setwd("<Path to GitHub Dir>")

###Import Data For Mutation Spectrum
Spectrum_Percent<-read.table("Figure_4/MutationSpectrum_Percent.txt",sep="\t",header=TRUE)
Spectrum_Percent
Spectrum_Percent.melt<-reshape2::melt(Spectrum_Percent, id=c("Type"))

names(Spectrum_Percent.melt)<-c("Type","Treatment","Percent")

Spectrum_Percent.melt$Type<-factor(Spectrum_Percent.melt$Type,levels=c("mobile_element_insertion","small_indel","intergenic","noncoding","nonsense","nonsynonymous","synonymous"))


Figure_4A_Mutation_Spectrum<-ggplot(data=Spectrum_Percent.melt, aes(x=Treatment,y=Percent,fill=Type))+
  geom_bar(stat="identity", color="black")+
  coord_cartesian(expand=FALSE)+
  ylab("Mutation Fraction")+
  scale_fill_manual(values=c("#5e4fa2","#3288bd","#66c2a5","#c5ff96","#Fee08b","#f46d43","#9e0142"),name="Mutation Type",labels=c("IS-element insertion","Small indel","Intergenic SNP","Noncoding SNP","Nonsense SNP","Nonsynonymous SNP","Synonymous SNP"))+
  theme_bw()+
  theme(legend.background = element_rect(color="black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank())

Figure_4A_Mutation_Spectrum


##Import Matrix for Bray Curtis analysis of mutations
BrayCurtis<-read.table("Figure_4/BrayCurtisCalc_Matrix.txt",sep="\t",header=TRUE)
BrayCurtis

BrayCurtis.meta<-read.table("Figure_4/BrayCurtisMeta_Matrix.txt",sep="\t",header=TRUE)

BrayCurtis.meta

rownames(BrayCurtis)<-BrayCurtis[,1]
BrayCurtis<-BrayCurtis[,-1]
BrayCurtis.t<-t(BrayCurtis)


##Figure_4B

Bray.back.md <- with(BrayCurtis.meta, meandist(vegdist(BrayCurtis.t,method="bray"), Background))
Bray.back.md
plot(Bray.back.md)

zz <- hclust(as.dist(Bray.back.md), method = "average")
dhc <- as.dendrogram(zz)
# Rectangular lines
ddata <- dendro_data(dhc, type = "rectangle")


Figure_4B_Bray_Distance_Dendogram <- ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), size=1)+
  geom_text(data = ddata$labels, aes(x = x, y = y-0.05, label = label,fontface="bold"), size = 3) +
  scale_y_continuous(limits=c(-0.1,1.1), breaks=c(0,0.2,0.4,0.6,0.8,1))+
  scale_x_continuous(limits=c(0.5,4.5))+
  ylab("Average Distance")+
  theme_bw()+
  theme(axis.title.y = element_text(face="bold"),axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(),panel.grid.minor=element_blank())
Figure_4B_Bray_Distance_Dendogram 

Figure_4_Top_Row<-plot_grid(Figure_4A_Mutation_Spectrum,Figure_4B_Bray_Distance_Dendogram, rel_widths = c(1,0.5), labels=c("A","B") )
Figure_4_Top_Row
##Figure_4C
BrayCurtis.Geno.ano <- with(BrayCurtis.meta, anosim(BrayCurtis.dist, Genotype))
summary(BrayCurtis.Geno.ano)
plot(BrayCurtis.Geno.ano)
Figure_4C_Geno.plot<-autoplot(BrayCurtis.Geno.ano,notch=FALSE)+
  aes(fill=Class)+
  ylab("Rank Distance (Bray-Curtis)")+
  scale_fill_manual(values=c("#656563","#95595A","#6C81AF"))+
  theme_bw()+
  theme(axis.title.y = element_text(face="bold"), plot.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none")
Figure_4C_Geno.plot


##Figure_4D
BrayCurtis.Treat.ano <- with(BrayCurtis.meta, anosim(BrayCurtis.dist, Treatment))
summary(BrayCurtis.Treat.ano)
BrayCurtis.Treat.ano$signif
Figure_4D_Treat.plot<-autoplot(BrayCurtis.Treat.ano,notch=FALSE)+
  aes(fill=Class)+
  ylab("Rank Distance (Bray-Curtis)")+
  scale_fill_manual(values=c("#656563","#bbbbbb","#888888"))+
  theme_bw()+
  theme(axis.title.y = element_text(face="bold"), plot.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none")

Figure_4D_Treat.plot


##Figure_4E
BrayCurtis.Planned.ano <- with(BrayCurtis.meta, anosim(BrayCurtis.dist, PlannedContrast))
summary(BrayCurtis.Planned.ano)
plot(BrayCurtis.Planned.ano)
Figure_4E_Plan.plot<-autoplot(BrayCurtis.Planned.ano,notch=FALSE)+
  aes(fill=Class)+
  ylab("Rank Distance (Bray-Curtis)")+
  scale_fill_manual(values=c("#656563","#92A8CA","#bbbbbb","#C08586"))+
  theme_bw()+
  theme(axis.title.y = element_text(face="bold"), plot.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none")

Figure_4E_Plan.plot

Figure_4_Middle_Row<-plot_grid(Figure_4C_Geno.plot,Figure_4D_Treat.plot,Figure_4E_Plan.plot,nrow=1,labels=c("C","D","E"))

plot_grid(Figure_4_Top_Row,Figure_4_Middle_Row, nrow=2)


###PCA of Populations

BrayCurtis
BrayCurtis$Total<-rowSums(BrayCurtis[,c("FCT_1", "FCT_2", "FCT_3" ,"FCT_4" ,"FCT_5" ,"FCT_6","WCT_1", "WCT_2", "WCT_3", "WCT_4", "WCT_5" ,"WCT_6","FFL_1", "FFL_2" ,"FFL_3", "FFL_4", "FFL_5" ,"FFL_6","WFL_1", "WFL_2" ,"WFL_3", "WFL_4", "WFL_5" ,"WFL_6")])

BrayCurtis.Gene<-BrayCurtis[which(BrayCurtis$Total>0),]
BrayCurtis.Gene<-BrayCurtis.Gene[,-25]
#BrayCurtis.Gene<-BrayCurtis
BrayCurtis.Gene

pca1=PCA(BrayCurtis.Gene, graph=TRUE)
pca1$eig
pca1$var

pca.dim<-pca1$var$coord
pca.dim
pca.dim<-as.data.frame(pca.dim)
pca.dim$Background<-c("FCT","FCT","FCT","FCT","FCT","FCT","FFL","FFL","FFL","FFL","FFL","FFL","WCT","WCT","WCT","WCT","WCT","WCT","WFL","WFL","WFL","WFL","WFL","WFL")
names(pca.dim)<-c("Dim.1","Dim.2","Dim.3","Dim.4","Dim.5","Background")
pca.dim

pca.dim$Dim.1<-as.numeric(pca.dim$Dim.1)
pca.dim$Dim.2<-as.numeric(pca.dim$Dim.2)


Figure_4F_PCAplot_Populations <- ggplot(data=pca.dim, aes(x=Dim.1,y=Dim.2,fill=Background,label=rownames(pca.dim)))+
  geom_point(size=2, pch=21,color="black")+
  scale_fill_manual(values=c("#a7c9df","#6AA0EF","#eda4a7","#c8686d"))+
  geom_text_repel(max.overlaps = Inf,size=2)+
  scale_y_continuous(limits=c(-0.25,0.9))+
  scale_x_continuous(limits=c(-0.23,0.9))+
  theme_bw()+
  theme(legend.position="none")
Figure_4F_PCAplot_Populations

gene.pca <- pca1$ind$coord
gene.pca<-as.data.frame(gene.pca)
#names(gene.mca)<-c("Dim.1","Dim.2","Dim.3","Dim.4","Dim.5","Genes")
gene.pca$Dim.1<-as.numeric(gene.pca$Dim.1)
gene.pca$Dim.2<-as.numeric(gene.pca$Dim.2)
gene.pca$Dim.3<-as.numeric(gene.pca$Dim.3)

Figure_4G_PCAplot_Genes<-ggplot(data=gene.pca, aes(x=Dim.1,y=Dim.2,label=rownames(gene.pca)))+
  geom_point(alpha=0.4, pch=23, fill="black", size=2)+
  geom_text_repel(size=2)+
  scale_y_continuous(limits=c(-3,15))+
  scale_x_continuous(limits=c(-3,23))+
  theme_bw()+
  theme(legend.position="none")



Figure_4_Bottom_Row<-plot_grid(Figure_4F_PCAplot_Populations,Figure_4G_PCAplot_Genes, labels=c("F","G"))
Figure_4_Bottom_Row

plot_grid(Figure_4_Top_Row,Figure_4_Middle_Row, Figure_4_Bottom_Row, nrow=3)
