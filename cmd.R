setwd("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/")
rm(list = ls())
library(DESeq2)
count=read.delim("ART-VP/Input.txt",row.names = 1, check.names = FALSE)
Design=read.delim("ART-VP/Info.txt")
ds=DESeqDataSetFromMatrix(countData = count,colData = Design,design = ~sex+batch+group)
ds=DESeq(ds)
res=results(ds,contrast = c("group","ART","VP"),independentFiltering = FALSE)
write.table(res,file="ART-VP/Result.txt",sep="\t", quote=FALSE,col.names = NA)
?results
###################
count=read.delim("EC-VP/Input.txt",row.names = 1, check.names = FALSE)
Design=read.delim("EC-VP/Info.txt")
ds=DESeqDataSetFromMatrix(countData = count,colData = Design,design = ~sex+batch+group)
ds=DESeq(ds)
res=results(ds,contrast = c("group","VP","EC"),independentFiltering = FALSE)
write.table(res,file="EC-VP/Result.txt",sep="\t", quote=FALSE,col.names = NA)

########################
count=read.delim("HC-ART/Input.txt",row.names = 1, check.names = FALSE)
Design=read.delim("HC-ART/Info.txt")
ds=DESeqDataSetFromMatrix(countData = count,colData = Design,design = ~sex+batch+group)
ds=DESeq(ds)
res=results(ds,contrast = c("group","ART","HC"),independentFiltering = FALSE)
write.table(res,file="HC-ART/Result.txt",sep="\t", quote=FALSE,col.names = NA)

############################
count=read.delim("HC-EC/Input.txt",row.names = 1, check.names = FALSE)
Design=read.delim("HC-EC/Info.txt")
ds=DESeqDataSetFromMatrix(countData = count,colData = Design,design = ~sex+batch+group)
ds=DESeq(ds)
res=results(ds,contrast = c("group","HC","EC"),independentFiltering = FALSE)
write.table(res,file="HC-EC/Result.txt",sep="\t", quote=FALSE,col.names = NA)

###########################
count=read.delim("HC-VP/Input.txt",row.names = 1, check.names = FALSE)
Design=read.delim("HC-VP/Info.txt")
ds=DESeqDataSetFromMatrix(countData = count,colData = Design,design = ~sex+batch+group)
ds=DESeq(ds)
res=results(ds,contrast = c("group","VP","HC"),independentFiltering = FALSE)
write.table(res,file="HC-VP/Result.txt",sep="\t", quote=FALSE,col.names = NA)

####################
count=read.delim("EC-ART/Input.txt",row.names = 1, check.names = FALSE)
Design=read.delim("EC-ART/Info.txt")
head(Design)
ds=DESeqDataSetFromMatrix(countData = count,colData = Design,design = ~ group + sex+group:nestedBatch+group:sex)
ds=DESeq(ds,modelMatrixType="standard")
res=results(ds,contrast = c("group","ART","EC"),independentFiltering = FALSE)
write.table(res,file="EC-ART/Result.txt",sep="\t", quote=FALSE,col.names = NA)
#https://support.bioconductor.org/p/64480/

install.packages("UpSetR")
library(library(PCAtools))
library(PCAtools)
library(Biobase)
library(edgeR)
count=read.delim("ART-Specific/PCA/ART_Coding_Count.txt",row.names = 1,check.names = FALSE)
meta=read.delim("ART-Specific/PCA/MetaData.txt",row.names = 1)
p <- pca(count, metadata = meta, removeVar = 0.1)
biplot(p,lab = FALSE,axisLabSize = 10,legendLabSize = 8,legendIconSize = 3,
       gridlines.major = FALSE, gridlines.minor = FALSE,legendPosition = "right",colby = "group")

head(p$rotated)
write.table(p$rotated,file="ART-Specific/PCA/PCA.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("ART-Specific/PCA/PCA.txt",row.names = 1,check.names = FALSE)
pdf("ART-Specific/PCA/PCA.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=group)) + geom_point(size=3.7,shape=21,aes(fill=group))+
  scale_color_manual(values=c(HC="#595959",VP="#99792f",EC="#315b7d",ART="#99665b"))+
  scale_fill_manual(values=c(HC="#808080",VP="#e5b647",EC="#4682b4",ART="#FFAA99"))+
  labs(x="PC1, 19.72% variation",y="PC2, 14.8% variation")+theme(axis.title = element_text(size=10))
dev.off()

library(gplots)

data=read.table("ART-Specific/PCA/ART_Coding_Count.txt",sep="\t",header=T)
rnames <- data[,1] 
mat_data <- data.matrix(data[,2:ncol(data)])
rownames(mat_data) <- rnames
my_palette <- colorRampPalette(c("#003d00","#19a419","white","#ff0000","#b20000"))(n=20)
sampleinfo <- read.delim("ART-Specific/PCA/MetaData.txt")
col.cell <- c(HC="#808080",VP="#e5b647",EC="#4682b4",ART="#FFAA99")[sampleinfo$group]
col.cell <- c("#FFAA99","#4682b4","#808080","#e5b647")[sampleinfo$group]
pdf("ART-Specific/PCA/HeatMap.pdf",width=32,height=35)
heatmap.2(mat_data,tracecol=NA,margins =c(11,16),col=my_palette,cexCol=1.5,
          ColSideColors=col.cell,labRow=FALSE,keysize = 1,Rowv = TRUE,scale="row",lhei=c(1,11))
dev.off()

############################ UpSet
setwd("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/")
library(UpSetR)
ART_VP=read.delim("ART-VP/Significant_Result.txt",row.names = 1,header = T)
EC_ART=read.delim("EC-ART/Significant_Result.txt",row.names = 1,header = T)
HC_ART=read.delim("HC-ART/Significant_Result.txt",row.names = 1,header = T)
EC_VP=read.delim("EC-VP/Significant_Result.txt",row.names = 1,header = T)
HC_VP=read.delim("HC-VP/Significant_Result.txt",row.names = 1,header = T)

list=list(VP_vs_ART=as.vector(ART_VP$GeneName),EC_vs_ART=as.vector(EC_ART$GeneName),HC_vs_ART=as.vector(HC_ART$GeneName),EC_vs_VP=as.vector(EC_VP$GeneName),HC_vs_VP=as.vector(HC_VP$GeneName))
pdf("Figure_1_a.pdf",width=10,height = 5)
upset(fromList(list), order.by = "freq",set_size.show = TRUE,set_size.scale_max=5000,main.bar.color = "#009999",sets.bar.color = "#009999")
dev.off()

####### #############################Log2CPM
setwd("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Batch_Effect_Exercise")
rm(list = ls())
data=read.table("Log2cpm.txt",sep="\t",header=T,row.names = 1)
dataX=t(data)
rnames <- data[,1] 
mat_data <- data.frame(dataX[,2:ncol(dataX)])
Art.umap = umap(mat_data)
head(art.umap,3)
write.table(Art.umap$layout,file="l2cpm_Umap.txt",sep="\t",col.names = NA,quote = FALSE)
library(ggplot2)
dat=read.delim("l2cpm_Umap.txt",row.names = 1,check.names = FALSE)
pdf("Suppl_Fig_1_a.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Batch)) + geom_point(size=4,shape=21,aes(fill=Batch))+
  theme(axis.title = element_text(size=9,hjust = 0.5),legend.position = c(0.08, 0.13),plot.margin = margin(0.7,0.5,0.7,0.5, "cm"))+labs(x = "UMAP1",y="UMAP2")
dev.off()

prot=read.delim("Log2cpm.txt",row.names = 1)
mat <- as.matrix(prot)
des=read.table("MetaData.txt",sep="\t",header=T)
designCombat = model.matrix(~ des$group)
rnaseqCombat = ComBat(mat, batch = des$batch, mod = designCombat, par.prior = TRUE, prior.plots = TRUE)
write.table(rnaseqCombat,file="lcpm/Combat_Result.txt", sep="\t",quote=FALSE,col.names = NA)

data=read.table("Combat_Result.txt",sep="\t",header=T,row.names = 1)
dataX=t(data)
rnames <- data[,1] 
mat_data <- data.frame(dataX[,2:ncol(dataX)])
Art.umap = umap(mat_data)
head(art.umap,3)
write.table(Art.umap$layout,file="Combat_Umap.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("Combat_Umap.txt",row.names = 1,check.names = FALSE)
pdf("Suppl_Fig_1_b.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Batch)) + geom_point(size=4,shape=21,aes(fill=Batch))+
  theme(axis.title = element_text(size=9,hjust = 0.5),legend.position = c(0.08, 0.87),plot.margin = margin(0.7,0.5,0.7,0.5, "cm"))+labs(x = "UMAP1",y="UMAP2")
dev.off()

pdf("Suppl_Fig_1_c.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=4,shape=21,aes(fill=Group))+
  scale_color_manual(values=c(VP="#b20000",HC="#0d3a1b",EC="#315b7d",ART="#b27300"))+
  scale_fill_manual(values=c(VP="#ff4c4c",HC="#96a94e",EC="#4682b4",ART="#FFA500"))+
  theme(axis.title = element_text(size=9,hjust = 0.5),legend.position = c(0.07, 0.88),plot.margin = margin(0.7,0.5,0.7,0.5, "cm"))+labs(x = "UMAP1",y="UMAP2")
dev.off()

pdf("ComBat_Sex.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Sex)) + geom_point(size=3.7,shape=21,aes(fill=Sex))+
  scale_color_manual(values=c(M="#315b7d", F="#997379"))+
  scale_fill_manual(values=c(M="#4682B4",F="#FFC0CB"))+
  theme(plot.title = element_text(size=15,hjust = 0.5))+labs(title = "Log2 cpm + ComBat: Protein Coding Genes")
dev.off()

count=read.delim("Combat_Result.txt",row.names = 1,check.names = FALSE)
countX=as.matrix(count)
calc_coef_var <- function(x) sd(x) / mean(x)
coef_var <- apply(countX, 1, calc_coef_var)
HVG_10 <- countX[rank(coef_var) / length(coef_var) > 0.9, ]
write.table(HVG_10,file="lcpm/Combat_Top.txt",sep="\t",col.names = NA,quote = FALSE)

data=read.table("Combat_Top.txt",sep="\t",header=T,row.names = 1)
dataX=t(data)
rnames <- data[,1] 
mat_data <- data.frame(dataX[,2:ncol(dataX)])
Art.umap = umap(mat_data)
head(art.umap,3)
write.table(Art.umap$layout,file="lcpm/Combat_Top_Umap.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Batch_Effect_Exercise/Combat_Top_Umap.txt",row.names = 1,check.names = FALSE)
pdf("Combat_Top_Batch.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Batch)) + geom_point(size=3.7,shape=21,aes(fill=Batch))+
  theme(plot.title = element_text(size=12,hjust = 0.5))+labs(title = "Log2 cpm + ComBat: Top 10% HVG")
dev.off()
library(ggplot2)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Batch_Effect_Exercise/Suppl_Fig_2_a.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=4,shape=21,aes(fill=Group))+
  scale_color_manual(values=c(VP="#b20000",HC="#0d3a1b",EC="#315b7d",ART="#b27300"))+
  scale_fill_manual(values=c(VP="#ff4c4c",HC="#96a94e",EC="#4682b4",ART="#FFA500"))+
  theme(axis.title = element_text(size=9,hjust = 0.5),legend.position = c(0.93, 0.12),plot.margin = margin(0.7,0.5,0.7,0.5, "cm"))+labs(x = "UMAP1",y="UMAP2")
dev.off()

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Batch_Effect_Exercise/Suppl_Fig_2_a.pdf")
ggplot(dat, aes(x=V1, y=V2)) + geom_point(size=5,aes(shape=Sex,color=Group,fill=Group))+ scale_shape_manual(values=c(21,22))+
  scale_color_manual(values=c(VP="#ff4c4c",HC="#96a94e",EC="#4682b4",ART="#FFA500"))+
  scale_fill_manual(values=c(VP="#ff4c4c",HC="#96a94e",EC="#4682b4",ART="#FFA500"))+
  theme(axis.title = element_text(size=15,hjust = 0.5),legend.position = c(0.92, 0.18),
        legend.title =element_text(size=15),legend.text = element_text(size = 15),plot.margin = margin(0.7,0.5,0.7,0.5, "cm"))+labs(x = "UMAP1",y="UMAP2")
dev.off()




pdf("Suppl_Fig_2_b.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Sex)) + geom_point(size=4,shape=21,aes(fill=Sex))+
  scale_color_manual(values=c(Male="#315b7d", Female="#997379"))+
  scale_fill_manual(values=c(Male="#4682B4",Female="#FFC0CB"))+
  theme(axis.title = element_text(size=9,hjust = 0.5),legend.position = c(0.91, 0.08),plot.margin = margin(0.7,0.5,0.7,0.5, "cm"))+labs(x = "UMAP1",y="UMAP2")
dev.off()

data=read.table("Combat_Top.txt",sep="\t",header=T,row.names = 1)
meta=read.table("MetaData.txt",sep="\t",header=T,row.names = 1)
P <- pca(data, metadata = meta, removeVar = 0.1)
biplot(P,colby = 'group',lab = FALSE,gridlines.major = FALSE, gridlines.minor = FALSE,
       legendPosition = 'right',axisLabSize = 10,legendLabSize = 8,legendIconSize = 3)+
  theme(plot.subtitle = element_text(hjust = 0.5))

###################### without ART
count=read.delim("withoutART/Combat_Result.txt",row.names = 1,check.names = FALSE)
countX=as.matrix(count)
calc_coef_var <- function(x) sd(x) / mean(x)
coef_var <- apply(countX, 1, calc_coef_var)
HVG_10 <- countX[rank(coef_var) / length(coef_var) > 0.9, ]
write.table(HVG_10,file="withoutART/Combat_Top.txt",sep="\t",col.names = NA,quote = FALSE)

data=read.table("withoutART/Combat_Top.txt",sep="\t",header=T,row.names = 1)
dataX=t(data)
rnames <- data[,1] 
mat_data <- data.frame(dataX[,2:ncol(dataX)])
Art.umap = umap(mat_data)
head(art.umap,3)
write.table(Art.umap$layout,file="withoutART/Combat_Top_Umap.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("withoutART/Combat_Top_Umap.txt",row.names = 1,check.names = FALSE)
pdf("Suppl_Fig_3_a.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=4,shape=21,aes(fill=Group))+
  scale_color_manual(values=c(VP="#b20000",HC="#0d3a1b",EC="#315b7d",ART="#b27300"))+
  scale_fill_manual(values=c(VP="#ff4c4c",HC="#96a94e",EC="#4682b4",ART="#FFA500"))+
  theme(axis.title = element_text(size=9,hjust = 0.5),legend.position = c(0.06, 0.099),plot.margin = margin(0.7,0.5,0.7,0.5, "cm"))+labs(x = "UMAP1",y="UMAP2")
dev.off()

pdf("Suppl_Fig_3_b.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Sex)) + geom_point(size=4,shape=21,aes(fill=Sex))+
  scale_color_manual(values=c(Male="#315b7d", Female="#997379"))+
  scale_fill_manual(values=c(Male="#4682B4",Female="#FFC0CB"))+
  theme(axis.title = element_text(size=9,hjust = 0.5),legend.position = c(0.09, 0.08),plot.margin = margin(0.7,0.5,0.7,0.5, "cm"))+labs(x = "UMAP1",y="UMAP2")
dev.off()
##########################  ART Specific
setwd("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Umap")
data=read.table("ART_Coding_Count.txt",sep="\t",header=T,row.names = 1)
dataX=t(data)
rnames <- data[,1] 
mat_data <- data.frame(dataX[,2:ncol(dataX)])
Art.umap = umap(mat_data)
head(art.umap,3)
write.table(Art.umap$layout,file="ART_Umap.txt",sep="\t",col.names = NA,quote = FALSE)
n

dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Umap/ART_Umap.txt",row.names = 1,check.names = FALSE)
pdf("ART_Batch.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Batch)) + geom_point(size=3.7,shape=21,aes(fill=Batch))+
  theme(plot.title = element_text(size=12,hjust = 0.5))+labs(title = "ART Specific genes")
dev.off()
library(ggrepel)
library(ggplot2)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Umap/Figure_2_a.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=7,shape=21,aes(fill=Group))+
  scale_color_manual(values=c(VP="#b20000",HC="#0d3a1b",EC="#315b7d",ART="#b27300"))+
  scale_fill_manual(values=c(VP="#ff4c4c",HC="#96a94e",EC="#4682b4",ART="#FFA500"))+
  theme(axis.title = element_text(size=20,hjust = 0.5),legend.title = element_blank(),legend.text = element_text(size=20),
        plot.margin = margin(0.7,0.5,0.7,0.5, "cm"),
        legend.position = c(0.91, 0.12))+labs(x = "UMAP1",y="UMAP2")#+
  #geom_text_repel(data=dat,aes(x=V1,y=V2,label = ID),size=1.9,box.padding=0.3,show.legend=FALSE,colour="black")
dev.off()

############### heatMap
setwd("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Batch_Effect_Exercise")
data=read.table("Log2cpm.txt",sep="\t",header=T)
rnames <- data[,1] 
mat_data <- data.matrix(data[,2:ncol(data)])
rownames(mat_data) <- rnames
my_palette <- colorRampPalette(c("#003d00","#19a419","white","#ff0000","#b20000"))(n=20)
sampleinfo <- read.delim("MetaData.txt")
col.cell <- c("#FFAA99","#4682b4","#808080","#e5b647")[sampleinfo$group]
library(gplots)
pdf("HeatMap.pdf",width=32,height=35)
heatmap.2(mat_data,tracecol=NA,margins =c(11,16),col=my_palette,cexCol=1.5,
          ColSideColors=col.cell,labRow=FALSE,keysize = 1,Rowv = TRUE,scale="row",lhei=c(1,11))
dev.off()



pdf("ART_Sex.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Sex)) + geom_point(size=3.7,shape=21,aes(fill=Sex))+
  scale_color_manual(values=c(M="#315b7d", F="#997379"))+
  scale_fill_manual(values=c(M="#4682B4",F="#FFC0CB"))+
  theme(plot.title = element_text(size=12,hjust = 0.5))+labs(title = "ART Specific genes")
dev.off()

setwd("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/True_ART/Umap")
data=read.table("64_Coding_Count.txt",sep="\t",header=T,row.names = 1)
dataX=t(data)
rnames <- data[,1] 
mat_data <- data.frame(dataX[,2:ncol(dataX)])
Art.umap = umap(mat_data)
head(art.umap,3)
write.table(Art.umap$layout,file="64_Umap.txt",sep="\t",col.names = NA,quote = FALSE)

data=read.table("64_50_Coding_Count.txt",sep="\t",header=T,row.names = 1)
dataX=t(data)
rnames <- data[,1] 
mat_data <- data.frame(dataX[,2:ncol(dataX)])
Art.umap = umap(mat_data)
head(art.umap,3)
write.table(Art.umap$layout,file="64_50_Umap.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("64_Umap.txt",row.names = 1,check.names = FALSE)
pdf("64_Group.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=3.7,shape=21,aes(fill=Group))+
  scale_color_manual(values=c(HC="#595959",VP="#99792f",EC="#315b7d",ART="#99665b"))+
  scale_fill_manual(values=c(HC="#808080",VP="#e5b647",EC="#4682b4",ART="#FFAA99"))+
  theme(plot.title = element_text(size=12,hjust = 0.5))+labs(title = "ART Specific genes")
dev.off()

dat=read.delim("64_50_Umap.txt",row.names = 1,check.names = FALSE)
pdf("64_50_Group.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=3.7,shape=21,aes(fill=Group))+
  scale_color_manual(values=c(HC="#595959",VP="#99792f",EC="#315b7d",ART="#99665b"))+
  scale_fill_manual(values=c(HC="#808080",VP="#e5b647",EC="#4682b4",ART="#FFAA99"))+
  theme(plot.title = element_text(size=12,hjust = 0.5))+labs(title = "ART Specific genes")
dev.off()

##### Pathway bargraphs
library(ggplot2)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/KEGG/BarPlot_input.txt",header = T)
#head(data)
data$Term <- factor(data$Term, levels = data$Term)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/KEGG/Figure_1_b.pdf",width=8,height = 10)
ggplot(data = data,aes(x=Term,y=overlap))+geom_bar(stat = "identity",fill="#687132")+coord_flip()+
  geom_text(aes(label=Pval), vjust=0.5, hjust=-0.1, color="black", size=3)+
  theme(axis.title.y = element_blank(),axis.text.y = element_text(size=12))+labs(y="Percentage of Gene Overlap")+ylim(0,25)
dev.off()

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/HC-ART/KEGG/BarPlot_Input.txt",header = T)
#head(data)
data$Term <- factor(data$Term, levels = data$Term)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/HC-ART/KEGG/Suppl_Fig_4_a.pdf",width=8,height = 10)
ggplot(data = data,aes(x=Term,y=overlap))+geom_bar(stat = "identity",fill="#86a3ac")+coord_flip()+
  geom_text(aes(label=Pval), vjust=0.5, hjust=-0.1, color="black", size=3)+
  theme(axis.title.y = element_blank(),axis.text.y = element_text(size=12))+labs(y="Percentage of Gene Overlap")+ylim(0,28)
dev.off()

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/KEGG/BarPlot_Input.txt",header = T)
#head(data)
data$Term <- factor(data$Term, levels = data$Term)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/KEGG/Figure_2_c.pdf",width=8,height = 10)
ggplot(data = data,aes(x=Term,y=overlap))+geom_bar(stat = "identity",fill="#7c977c")+coord_flip()+
  geom_text(aes(label=Pval), vjust=0.5, hjust=-0.1, color="black", size=3)+
  theme(axis.title.y = element_blank(),axis.text.y = element_text(size=12))+labs(y="Percentage of Gene Overlap")+ylim(0,28)
dev.off()



library(UpSetR)
setwd("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/")
Dat=read.delim("Tmp.txt",header = T,check.names = FALSE)
head(Dat)
pdf("Upset.pdf",width=10,height = 5)
upset(Dat, nsets = 6,sets.x.label="Number of regulated genes",order.by = "freq",set_size.show = TRUE,query.legend = "bottom",set_size.scale_max=5000,main.bar.color = "#4682b4",sets.bar.color="#4682b4",
      queries = list(list(query = intersects, params = list("EC_vs_ART"), color = "#ffb910",active = T,query.name="Unique Genes : EC-vs-ART"), 
                     list(query = intersects, params = list("HC_vs_ART"), color = "#800000", active = T,query.name="Unique Genes : HC-vs-ART"),
                     list(query = intersects, params = list("ART","TRUE"), color = "red", active = F,query.name="ART Specific Genes")))
dev.off()

setwd("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/True_ART/KEGG_EC_ART")

data=read.delim("Input.txt",header = T)
head(data)
library(ggplot2)
data$Term <- factor(data$Term, levels = data$Term)
pdf("Pathway.pdf")
ggplot(data, aes(x = Term, y = Overlap,size=pval))+geom_point(data=data,aes(color=pval))+coord_flip()+
  scale_color_gradient(low="#a6a6a6", high="#446655")+scale_size(range = c(2, 5))+labs(y="Number of Genes Regulated",color="-log10(Adj.Pvalue)",size="-log10(Adj.Pvalue)")+
  theme(axis.title.y = element_blank(),legend.title = element_text(size=7),legend.key.size = unit(0.6,"line"),legend.text = element_text(size=7),
        axis.title.x = element_text(size=7.5),axis.text = element_text(size=7),plot.margin = margin(5,1.2,5,1.2, "cm"))+
  guides(size=guide_legend(override.aes=list(colour="grey")))
dev.off()

setwd("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/KEGG")
data=read.delim("Input.txt",header = T)
head(data)
library(ggplot2)
data$Term <- factor(data$Term, levels = data$Term)
pdf("Pathway.pdf")
ggplot(data, aes(x = Term, y = Overlap,size=pval))+geom_point(data=data,aes(color=pval))+coord_flip()+
  scale_color_gradient(low="#a6a6a6", high="#6CA1E4")+scale_size(range = c(2, 5))+labs(y="Number of Genes Regulated",color="-log10(Adj.Pvalue)",size="-log10(Adj.Pvalue)")+
  theme(axis.title.y = element_blank(),legend.title = element_text(size=8),legend.key.size = unit(0.6,"line"),legend.text = element_text(size=7),
        axis.title.x = element_text(size=7.5),axis.text = element_text(size=7),plot.margin = margin(5,1.2,5,1.2, "cm"))+
  guides(size=guide_legend(override.aes=list(colour="grey")))
dev.off()

setwd("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/True_ART/KEGG_HC_ART")
data=read.delim("Input.txt",header = T)
data$Term <- factor(data$Term, levels = data$Term)
pdf("Pathway.pdf")
ggplot(data, aes(x = Term, y = Overlap,size=pval))+geom_point(data=data,aes(color=pval))+coord_flip()+
  scale_color_gradient(low="#a6a6a6", high="#00CC99")+scale_size(range = c(2, 6))+labs(y="Number of Genes Regulated",color="-log10(Pvalue)",size="-log10(Pvalue)")+
  theme(axis.title.y = element_blank(),legend.title = element_text(size=7),legend.key.size = unit(0.6,"line"),legend.text = element_text(size=7),
        axis.title.x = element_text(size=7.5),axis.text = element_text(size=7),plot.margin = margin(5,1.2,5,1.2, "cm"))+
  guides(size=guide_legend(override.aes=list(colour="grey")))
dev.off()

setwd("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Umap")
library(gplots)
rm(list = ls())
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Batch_Effect_Exercise/Input.txt",sep="\t",header=T)
data=t(dat)
write.table(data,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Batch_Effect_Exercise/IP_HeatMap.txt",sep="\t",col.names = NA,quote = FALSE)
head(dat)
rnames <- data[,1] 
mat_data <- data.matrix(data[,2:ncol(data)])
rownames(mat_data) <- rnames
my_palette <- colorRampPalette(c("#00FFFF","#32ffff","black","#cccc00","#7f7f00"))(n=20)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Batch_Effect_Exercise/Info_HeatMap.txt",row.names = 1)

Group <- c("ART"="#FFAA99","EC"="#4682b4","HC"="#808080","VP"="#e5b647")[sampleinfo$group]
Sex <- c("Female"="#ffc0cb","Male"="#296DC3")[sampleinfo$sex]
Ethnicity<- c("Black"="#3e2b13","Caucasian"="#eac086","Latin"="#fffbae")[sampleinfo$Ethnicity]
clab=cbind(Group,Sex)
class(clab)
nrow(mat_)
colnames(clab)=c("Group","Sex")
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}
head(clab)
nrow(clab)
my_palette <- colorRampPalette(c("#00FFFF","#7f7f7f","#ffff00"))(n=25)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Batch_Effect_Exercise/HeatMap.pdf",width=37,height=35)
heatmap.3(mat_data,tracecol=NA,margins  =c(8,35),lhei=c(1,11),lwid =c(2,15),col=my_palette,cexCol=1.5,
          labRow=TRUE,keysize = 1.2,Rowv = TRUE,hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="row",
          dendrogram="row",Colv=FALSE, ColSideColors=clab,ColSideColorsSize=4,labCol = FALSE)
legend("topright",legend=c("Black","Caucasian","Latin","","Female","Male","","ART","EC","HC","VP"),
       fill=c("#3e2b13","#eac086","#fffbae","white","#ffc0cb","#296DC3","white","#FFAA99","#4682b4","#808080","#e5b647"), border=FALSE, bty="n", x.intersp=0.2,y.intersp = 0.7, cex=3)
dev.off()
?heatmap.3
library("devtools")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")



data(mtcars)
x  <- as.matrix(mtcars)
rc <- rainbow(nrow(x), start=0, end=.3)
cc <- rainbow(ncol(x), start=0, end=.3)
pdf("temp.pdf",width=37,height=35)
heatmap.3(x,margins  =c(8,35),lhei=c(1,11),lwid =c(2,15))
dev.off()
getwd()

library(countToFPKM)
file.readcounts <- system.file("extdata", "RNA-seq.read.counts.csv", package="countToFPKM")
file.annotations <- system.file("extdata", "Biomart.annotations.hg38.txt", package="countToFPKM")
file.sample.metrics <- system.file("extdata", "RNA-seq.samples.metrics.txt", package="countToFPKM")
gene.annotations <- read.table(file.annotations, sep="\t", header=TRUE)
featureLength <- gene.annotations$length

head(gene.annotations)
head(featureLength)
samples.metrics <- read.table(file.sample.metrics, sep="\t", header=TRUE)

head(samples.metrics)
library(NOISeq)
getwd()
setwd("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Correlations")
data=read.delim("ART_Count.txt",header = T,row.names = 1)
dataX=as.matrix(data)
norm=tmm(dataX)
head(norm)

write.table(norm,file="TMM_Normed.txt",sep="\t",col.names = NA,quote = FALSE)
library(psych)
clini=read.csv("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Correlations/New/Input_Clinical.txt",check.names=FALSE)
head(clini)
gene=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Correlations/New/Input_Genes.txt",check.names=FALSE)
cliniX=as.matrix(clini)
geneX=as.matrix(gene)
Res=corr.test(cliniX,geneX,use = "pairwise",method="spearman",adjust="bonferroni")
?corr.test
write.table(Res$r,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Correlations/New/FPKM/Correlation_Pair.txt",sep="\t",quote=FALSE)
write.table(Res$p,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Correlations/New/FPKM/Probabilty__Pair.txt",sep="\t",quote=FALSE)
tmp=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Correlations/New/FPKM/Correlation_Pair.txt",check.names = FALSE)
head(tmp)
library(reshape2)
m=melt(tmp,id.vars = "ID")
head(m)
write.table(m,file ="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Correlations/New/FPKM/Correlation.txt",col.names = NA,quote = FALSE,sep="\t")

tmp=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Correlations/New/FPKM/Probabilty__Pair.txt",check.names = FALSE)
head(tmp)
library(reshape2)
m=melt(tmp,id.vars = "ID")
head(m)
write.table(m,file ="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Correlations/New/FPKM/Pvalues.txt",col.names = NA,quote = FALSE,sep="\t")


#################### Chord diagram
library(circlize)
data=read.csv("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/KEGG/Chord_Input.txt")
head(data)
col_fun = colorRamp2(c(-1,0,1), c("#004c00", "white","#ae0001"), transparency = 0)
grid.col = c(ATP6V1G1="#A67C00",NDUFA1="#A67C00",ATP5MC2="#A67C00",ATP5PF="#A67C00",ATP5PD="#A67C00",NDUFA6="#A67C00",COX7C="#A67C00",NDUFA4="#A67C00",COX11="#A67C00",COX5A="#A67C00",COX10="#A67C00",COX7B="#A67C00",ATP5MC1="#A67C00",NDUFA13="#A67C00",UQCRFS1="#A67C00",ATP5F1E="#A67C00",COX7A2="#A67C00",NDUFB5="#A67C00",COX6A1="#A67C00",UQCRB="#A67C00",NDUFB4="#A67C00",NDUFS3="#A67C00",COX6C="#A67C00",NDUFS5="#A67C00",ATP5MC3="#A67C00",NDUFS4="#A67C00",SDHA="#A67C00",UQCRH="#A67C00",NDUFB1="#A67C00",COX6B1="#A67C00",NDUFA1="#A67C00",ATP5PF="#A67C00",ATP5MC2="#A67C00",ATP5PD="#A67C00",SOS1="#A67C00",NDUFA6="#A67C00",COX7C="#A67C00",ARID1A="#A67C00",NDUFA4="#A67C00",RPS6KA3="#A67C00",COX11="#A67C00",COX14="#A67C00",GNAS="#A67C00",COX5A="#A67C00",COX10="#A67C00",COX7B="#A67C00",ATP5MC1="#A67C00",NDUFA13="#A67C00",ADCY7="#A67C00",ADCY4="#A67C00",UQCRFS1="#A67C00",ATP5F1E="#A67C00",COX7A2="#A67C00",NDUFB5="#A67C00",COA1="#A67C00",COX6A1="#A67C00",UQCRB="#A67C00",PPARGC1A="#A67C00",PRKG2="#A67C00",NDUFB4="#A67C00",NDUFS3="#A67C00",COX6C="#A67C00",CPT1B="#A67C00",NDUFS5="#A67C00",ATP5MC3="#A67C00",NDUFS4="#A67C00",SDHA="#A67C00",MAPK12="#A67C00",UQCRH="#A67C00",NDUFB1="#A67C00",COX6B1="#A67C00",Thermogenesis="#FFA560",Oxidative_phosphorylation="#99CCFF")
grid.col= c(ATP6V1G1="#A67C00",NDUFA1="#A67C00",ATP5MC2="#A67C00",ATP5PF="#A67C00",ATP5PD="#A67C00",NDUFA6="#A67C00",COX7C="#A67C00",NDUFA4="#A67C00",COX11="#bfbfbf",COX5A="#A67C00",COX10="#bfbfbf",COX7B="#A67C00",ATP5MC1="#A67C00",NDUFA13="#bfbfbf",UQCRFS1="#A67C00",ATP5F1E="#A67C00",COX7A2="#A67C00",NDUFB5="#bfbfbf",COX6A1="#A67C00",UQCRB="#A67C00",NDUFB4="#A67C00",NDUFS3="#A67C00",COX6C="#A67C00",NDUFS5="#A67C00",ATP5MC3="#A67C00",NDUFS4="#A67C00",SDHA="#bfbfbf",UQCRH="#A67C00",NDUFB1="#A67C00",COX6B1="#A67C00",NDUFA1="#A67C00",ATP5PF="#A67C00",ATP5MC2="#A67C00",ATP5PD="#A67C00",SOS1="#A67C00",NDUFA6="#A67C00",COX7C="#A67C00",ARID1A="#A67C00",NDUFA4="#A67C00",RPS6KA3="#bfbfbf",COX11="#A67C00",COX14="#A67C00",GNAS="#A67C00",COX5A="#A67C00",COX10="#A67C00",COX7B="#A67C00",ATP5MC1="#A67C00",NDUFA13="#A67C00",ADCY7="#A67C00",ADCY4="#bfbfbf",UQCRFS1="#A67C00",ATP5F1E="#A67C00",COX7A2="#A67C00",NDUFB5="#A67C00",COA1="#A67C00",COX6A1="#A67C00",UQCRB="#A67C00",PPARGC1A="#A67C00",PRKG2="#bfbfbf",NDUFB4="#A67C00",NDUFS3="#A67C00",COX6C="#A67C00",CPT1B="#bfbfbf",NDUFS5="#A67C00",ATP5MC3="#A67C00",NDUFS4="#A67C00",SDHA="#A67C00",MAPK12="#bfbfbf",UQCRH="#A67C00",NDUFB1="#A67C00",COX6B1="#A67C00",Thermogenesis="#FFA560",Oxidative_phosphorylation="#99CCFF")
chordDiagram(data)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/KEGG/test.pdf")
chordDiagram(data, annotationTrack = "grid",col=col_fun,grid.col=grid.col,
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(data))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex=0.5)
}, bg.border = NA)
dev.off()

############ Correlation Heatmap
setwd("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Correlations/")
library(gplots)
library(RColorBrewer)
input=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Correlations/New/FPKM/Correlation_Pair.txt",header = T,check.names = FALSE)
head(input)
rnames <- input[,1]
mat_data <- data.matrix(input[,2:ncol(input)])
rownames(mat_data) <- rnames
head(mat_data)
not=read.csv("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Correlations/New/pval.txt",header = T,check.names = FALSE,sep = ",",stringsAsFactors = FALSE,na.strings = FALSE,row.names = 1)


pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Correlations/New/FPKM/HeatMap.pdf",width =13,height=5)
heatmap.2(mat_data,trace="none",cellnote = not,notecol="#4c4c4c",notecex=2,Rowv = FALSE,cexRow=1.1,cexCol=1.1,col=brewer.pal(n = 9, name = "YlGnBu"),lwid  = c(1,8),dendrogram = "none",margins = c(7, 15),key.title = NA,key.xlab = NA,keysize = 1,key.ylab = NA)
dev.off()

head(rnames)

############ cluster HeatMap
library(gplots)
rm(list = ls())
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Batch_Effect_Exercise/Input.txt",sep="\t",header=T,row.names = 1)
rnames <- data[,1]
head(mat_data)
mat_data <- data.matrix(data[,2:ncol(data)])
rownames(mat_data) <- rnames
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Batch_Effect_Exercise/Info_HeatMap.txt",row.names = 1)

Group <- c("ART"="#FFAA99","EC"="#4682b4","HC"="#808080","VP"="#e5b647")[sampleinfo$group]
Sex <- c("Female"="#ffc0cb","Male"="#296DC3")[sampleinfo$sex]

clab=cbind(Group,Sex)

colnames(clab)=c("Group","Sex")
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}
nrow(clab)
my_palette <- colorRampPalette(c("#00FFFF","#7f7f7f","#ffff00"))(n=25)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Batch_Effect_Exercise/HeatMap.pdf",width=37,height=35)
heatmap.3(mat_data,tracecol=NA,margins  =c(8,35),lhei=c(1,11),lwid =c(2,15),col=my_palette,cexCol=1.5,
          labRow=TRUE,keysize = 1.2,Rowv = TRUE,hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="row",
          dendrogram="row",Colv=FALSE, ColSideColors=clab,ColSideColorsSize=4,labCol = FALSE)
legend("topright",legend=c("Black","Caucasian","Latin","","Female","Male","","ART","EC","HC","VP"),
       fill=c("#3e2b13","#eac086","#fffbae","white","#ffc0cb","#296DC3","white","#FFAA99","#4682b4","#808080","#e5b647"), border=FALSE, bty="n", x.intersp=0.2,y.intersp = 0.7, cex=3)
dev.off()


#######################  Group wise UMAP

data=read.table("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Batch_Effect_Exercise/GroupWiseUMAP/Count_lcpm_combat.txt",sep="\t",header=T,row.names = 1)
dataX=t(data)
rnames <- data[,1] 
mat_data <- data.frame(dataX[,2:ncol(dataX)])
Art.umap = umap(mat_data)
head(art.umap,3)
write.table(Art.umap$layout,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Batch_Effect_Exercise/GroupWiseUMAP/ART_Umap.txt",sep="\t",col.names = NA,quote = FALSE)


dat=read.csv("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Batch_Effect_Exercise/GroupWiseUMAP/ART_Umap.txt",row.names = 1,check.names = FALSE)
library(ggplot2)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Batch_Effect_Exercise/GroupWiseUMAP/ART_Umap.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Gender)) + geom_point(size=4,shape=21,aes(fill=Gender))+
  scale_color_manual(values=c(Male="#315b7d", Female="#997379"))+
  scale_fill_manual(values=c(Male="#4682B4",Female="#FFC0CB"))+
  theme(axis.title = element_text(size=9,hjust = 0.5),legend.position = c(0.09, 0.08),plot.margin = margin(0.7,0.5,0.7,0.5, "cm"))+labs(x = "UMAP1",y="UMAP2")
dev.off()

data=read.table("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Batch_Effect_Exercise/GroupWiseUMAP/EC.txt",sep="\t",header=T,row.names = 1)
dataX=t(data)
rnames <- data[,1] 
mat_data <- data.frame(dataX[,2:ncol(dataX)])
Art.umap = umap(mat_data)
head(art.umap,3)
write.table(Art.umap$layout,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Batch_Effect_Exercise/GroupWiseUMAP/EC_Umap.txt",sep="\t",col.names = NA,quote = FALSE)


dat=read.csv("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Batch_Effect_Exercise/GroupWiseUMAP/EC_Umap.txt",row.names = 1,check.names = FALSE)
library(ggplot2)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Batch_Effect_Exercise/GroupWiseUMAP/EC_Umap.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Gender)) + geom_point(size=4,shape=21,aes(fill=Gender))+
  scale_color_manual(values=c(Male="#315b7d", Female="#997379"))+
  scale_fill_manual(values=c(Male="#4682B4",Female="#FFC0CB"))+
  theme(axis.title = element_text(size=9,hjust = 0.5),legend.position = c(0.91, 0.08),plot.margin = margin(0.7,0.5,0.7,0.5, "cm"))+labs(x = "UMAP1",y="UMAP2")
dev.off()


library(countToFPKM)
getwd()
file.readcounts <- system.file("extdata", "RNA-seq.read.counts.csv", package="countToFPKM")
file.annotations <- system.file("extdata", "Biomart.annotations.hg38.txt", package="countToFPKM")
file.sample.metrics <- system.file("extdata", "RNA-seq.samples.metrics.txt", package="countToFPKM")

counts <- as.matrix(read.csv(file.readcounts))
head(counts)
gene.annotations <- read.table(file.annotations, sep="\t", header=TRUE)
featureLength <- gene.annotations$length
head(gene.annotations)
samples.metrics <- read.table(file.sample.metrics, sep="\t", header=TRUE)




cnt=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/NetWorkAnalysis/ART_NetWork/Coding_Count_Hg38.txt",row.names = 1)
head(cnt,2)
anno=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/NetWorkAnalysis/ART_NetWork/GeneLength.txt")
head(anno)
sam=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/NetWorkAnalysis/ART_NetWork/FragLen.txt")
head(sam)
meanFragmentLength <- sam$meanFragmentLength
featureLength <- anno$length
fpkm_matrix <- fpkm (as.matrix(cnt), featureLength, meanFragmentLength)
head(fpkm_matrix)
write.table(fpkm_matrix,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/NetWorkAnalysis/ART_NetWork/FPKM.txt",sep="\t",quote = FALSE,col.names = NA)
?fpkm


###################################
library(psych)
?corr.test

gene=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/ART_Correlation/ART_TMM.txt",check.names=FALSE,row.names = 1)
geneX=as.matrix(t(gene))
Res=corr.test(geneX,use = "pairwise",method="spearman",adjust="BH")
write.table(Res$r,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Correlations/New/Correlation_Pair.txt",sep="\t",quote=FALSE)
write.table(Res$p,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Correlations/New/Probabilty__Pair.txt",sep="\t",quote=FALSE)
tmp=read.csv("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Correlations/New/Correlation_Pair.txt",check.names = FALSE)
head(tmp)
library(reshape2)
m=melt(tmp,id.vars = "ID")
head(m)
write.table(m,file ="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Correlations/New/Correlation.txt",col.names = NA,quote = FALSE,sep="\t")


############################# Non Coding Analysis #########################
setwd("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/")
rm(list = ls())
library(DESeq2)
count=read.delim("VP_ART//Input.txt",row.names = 1, check.names = FALSE)
Design=read.delim("VP_ART/Info.txt")
ds=DESeqDataSetFromMatrix(countData = count,colData = Design,design = ~sex+batch+group)
ds=DESeq(ds)
res=results(ds,contrast = c("group","ART","VP"),independentFiltering = FALSE)
write.table(res,file="VP_ART/Result.txt",sep="\t", quote=FALSE,col.names = NA)
?results
###################
count=read.delim("EC_VP/Input.txt",row.names = 1, check.names = FALSE)
Design=read.delim("EC_VP/Info.txt")
ds=DESeqDataSetFromMatrix(countData = count,colData = Design,design = ~sex+batch+group)
ds=DESeq(ds)
res=results(ds,contrast = c("group","VP","EC"),independentFiltering = FALSE)
write.table(res,file="EC_VP/Result.txt",sep="\t", quote=FALSE,col.names = NA)

########################
count=read.delim("HC_ART/Input.txt",row.names = 1, check.names = FALSE)
Design=read.delim("HC_ART/Info.txt")
ds=DESeqDataSetFromMatrix(countData = count,colData = Design,design = ~sex+batch+group)
ds=DESeq(ds)
res=results(ds,contrast = c("group","ART","HC"),independentFiltering = FALSE)
write.table(res,file="HC_ART/Result.txt",sep="\t", quote=FALSE,col.names = NA)

############################
count=read.delim("HC_EC/Input.txt",row.names = 1, check.names = FALSE)
Design=read.delim("HC_EC/Info.txt")
ds=DESeqDataSetFromMatrix(countData = count,colData = Design,design = ~sex+batch+group)
ds=DESeq(ds)
res=results(ds,contrast = c("group","EC","HC"),independentFiltering = FALSE)
write.table(res,file="HC_EC/Result.txt",sep="\t", quote=FALSE,col.names = NA)

###########################
count=read.delim("HC_VP/Input.txt",row.names = 1, check.names = FALSE)
Design=read.delim("HC_VP/Info.txt")
ds=DESeqDataSetFromMatrix(countData = count,colData = Design,design = ~sex+batch+group)
ds=DESeq(ds)
res=results(ds,contrast = c("group","VP","HC"),independentFiltering = FALSE)
write.table(res,file="HC_VP/Result.txt",sep="\t", quote=FALSE,col.names = NA)

####################
count=read.delim("EC_ART/Input.txt",row.names = 1, check.names = FALSE)
Design=read.delim("EC_ART/Info.txt")
head(Design)
ds=DESeqDataSetFromMatrix(countData = count,colData = Design,design = ~ group + sex+group:nestedBatch+group:sex)
ds=DESeq(ds,modelMatrixType="standard")
res=results(ds,contrast = c("group","ART","EC"),independentFiltering = FALSE)
write.table(res,file="EC_ART/Result.txt",sep="\t", quote=FALSE,col.names = NA)

library(LncPath)
NetLncPath <- getNet()
print(head(NetLncPath), row.names = FALSE)
SigLncs <- getExampleData("SigLncs")
print(head(SigLncs), row.names = FALSE)

linc=read.delim("EC_ART/tmp.txt",header = T)
head(linc)
avector <- as.character(linc[['ID']])
class(avector)
Result <- lncPath(avector, NetLncPath, Weighted = TRUE, PathwayDataSet = "KEGG", nperm = 2000,minPathSize = 0, maxPathSize = 500)
names(Result$KEGG_GLYCOLYSIS_GLUCONEOGENESIS)
head(Result$KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY$Obs.ES)
PathwaySummaryTable <- lncPath2Table(Result)
head(PathwaySummaryTable)
write.table(PathwaySummaryTable,file="EC_ART/Pwy_linRNA.txt",sep="\t",col.names = NA,quote = FALSE)
?lncPath
class(Result)
head(Result$KEGG_OXIDATIVE_PHOSPHORYLATION)
.libPaths()

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/KEGG_Selected/pwy.txt",header = T,check.names = FALSE)
head(data)
library(ggplot2)
data$Term <- factor(data$Term, levels = data$Term)
ggplot(data=data, aes(x=Term, y=gene)) +geom_bar(stat="identity", fill="steelblue")+geom_text(aes(label=pval), vjust=-0.3, size=3.5)+
  theme_minimal()+coord_flip()+theme(axis.title.y = element_blank())+ylab("# Regulated genes")


### Network analysis

library(ggplot2)
ip=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/NetWorkAnalysis/From_rui/network analysis/results/art/association_analysis/New/CentralComm/PWY_All/Input.txt",header = T)
ip$Term <- factor(ip$Term, levels = ip$Term)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/NetWorkAnalysis/From_rui/network analysis/results/art/association_analysis/New/CentralComm/PWY_All/Pathway.pdf")
ggplot(ip, aes(y=Term)) + 
  geom_point(data=ip,aes(x=1,y=Term,size=Gene,color=pval))+scale_x_discrete(limits=c("1"))+scale_color_gradient(high="#bfdbd8",low="#1e5f58",breaks=c(0.0026,0.0087,0.038,0.14,0.16),limits=c(0.0026,0.2))+
  scale_y_discrete(position = "right")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=9,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(4,3.7,5.2,5, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "# Genes"),color=guide_legend(title = "Adj.Pvalue"))
dev.off()

ip=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/NetWorkAnalysis/From_rui/network analysis/results/global/association_analysis/MostCentral/PWY/Input.txt",header = T)
ip$Term <- factor(ip$Term, levels = ip$Term)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/NetWorkAnalysis/From_rui/network analysis/results/global/association_analysis/MostCentral/PWY/Pathway.pdf")
ggplot(ip, aes(y=Term)) + 
  geom_point(data=ip,aes(x=1,y=Term,size=Gene,color=pval))+scale_x_discrete(limits=c("1"))+scale_color_gradient(high="#99e1ff",low="#006c99",breaks=c(0.00031,0.0004,0.005,0.007,0.01),limits=c(0.00031,0.01))+
  scale_y_discrete(position = "right")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=9,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(3,5,4,5.5, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),legend.position = c(4.5, -0.2),legend.box="vertical",
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "# Genes",nrow = 1),color=guide_legend(title = "Adj.Pvalue",nrow=2))
dev.off()


library(gplots)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/KEGG_Selected/HeatMap/Oxphos.txt",sep="\t",header=T,row.names = 1)
mat_data <- data.matrix(data[,1:ncol(data)])
head(mat_data)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/KEGG_Selected/HeatMap/Info.txt",row.names = 1)

Group <- c("#FFA500","#4682b4","#96a94e","#ff4c4c")[sampleinfo$group]

sampleinfo$group
ncol(mat_data)
my_palette <- colorRampPalette(c("#00FFFF","#7f7f7f","#ffff00"))(n=25)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/KEGG_Selected/HeatMap/HeatMap.pdf",width=37,height=35)
heatmap.2(mat_data,tracecol=NA,margins  =c(95,15),lhei=c(0.7,5),lwid =c(2,15),col=my_palette,cexRow=3,cexCol = 1.5,keysize = 1,Rowv = TRUE, na.rm = TRUE,
          ColSideColors=Group,scale="row",key.title=NA,labCol = FALSE)
dev.off()


##### Regression ####

library(tidyverse)
library(ggpubr)
data("marketing", package = "datarium")
install.packages("datarium")
head(marketing)
setwd("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Correlations/New/Regression")

data=read.delim("Data.txt",header = T, row.names = 1,check.names = FALSE)
head(data)
ggplot(data, aes(data$Age, data$NDUFS5)) + geom_point() + stat_smooth(method = lm)

model <- lm(NDUFS5 ~ Age, data = data)
S=summary(model)
S$adj.r.squared

library(FSA)
Data=read.csv("Result.txt",header = T)
head(Data)
Data$Bonferroni =p.adjust(Data$pvalue,method = "bonferroni")
Data$BH =p.adjust(Data$pvalue,method = "BH")
write.table(Data,file="Regression.txt",sep="\t",col.names = NA,quote = FALSE)

clini=read.delim("Input_Clinical.txt",check.names=FALSE)
head(geneX)
gene=read.delim("Input_Genes.txt",check.names=FALSE)
cliniX=as.matrix(clini)
geneX=as.matrix(gene)
Res=corr.test(cliniX,geneX,use = "pairwise",method="spearman",adjust="bonferroni")
?corr.test
write.table(Res$r,file="Correlation_Pair.txt",sep="\t",quote=FALSE)
write.table(Res$p,file="Probabilty__Pair.txt",sep="\t",quote=FALSE)
tmp=read.delim("Correlation_Pair.txt",check.names = FALSE)
head(tmp)
library(reshape2)
m=melt(tmp,id.vars = "ID")
head(m)
write.table(m,file ="Correlation.txt",col.names = NA,quote = FALSE,sep="\t")

tmp=read.delim("Probabilty__Pair.txt",check.names = FALSE)
head(tmp)
library(reshape2)
m=melt(tmp,id.vars = "ID")
head(m)
write.table(m,file ="Pvalues.txt",col.names = NA,quote = FALSE,sep="\t")
?lm()

library(countToFPKM)
cnt=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/NonCoding_Count_Hg38.txt",row.names = 1)
head(cnt,2)
anno=read.csv("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/ART_Specific/GeneLength.txt")
head(anno)
sam=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/NetWorkAnalysis/ART_NetWork/FragLen.txt")
head(sam)
meanFragmentLength <- sam$meanFragmentLength
featureLength <- anno$length
fpkm_matrix <- fpkm (as.matrix(cnt), featureLength, meanFragmentLength)
?fpkm
write.table(fpkm_matrix,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/FPKM.txt",sep="\t",quote = FALSE,col.names = NA)


library(ggplot2)
ip=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/ART_Specific/For_figure.txt",header = T)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/ART_Specific/Pathway.pdf")
ggplot(ip, aes(y=Term)) + 
  geom_point(data=ip,aes(x=1,y=Term,size=Gene,color=Correlation))+scale_x_discrete(limits=c("1"))+scale_color_manual(values=c("#007300","#FF704C"))+   #+scale_color_gradient(high="#99e1ff",low="#006c99",breaks=c(0.00031,0.0004,0.005,0.007,0.01),limits=c(0.00031,0.01))+
  scale_y_discrete(position = "right")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=9,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(3,5,4,5.5, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),legend.position = c(4.5, -0.2),legend.box="vertical",
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "# Genes",nrow = 1),color=guide_legend(title = "Correlation",nrow=1))
dev.off()




library(gplots)
data=read.delim("/home/anoop/Desktop/covid/Coagulation_cascade_pwy/heatmap/Pwy_RNASeq.txt",sep="\t",header=T,row.names = 1)
mat_data <- data.matrix(data[,1:ncol(data)])
head(mat_data)
sampleinfo <- read.delim("/home/anoop/Desktop/covid/Coagulation_cascade_pwy/heatmap/design.txt",row.names = 1)

Group <- c("#FFA500","#4682b4","#96a94e","#ff4c4c")[sampleinfo$group]

sampleinfo$group
ncol(mat_data)
my_palette <- colorRampPalette(c("#7e3f12","#d67834","#fef2c1","#3f75a2","#23415a"))(n=30)
pdf("/home/anoop/Desktop/covid/Coagulation_cascade_pwy/heatmap/HeatMap.pdf",width=37,height=35)
heatmap.2(mat_data,tracecol=NA,margins  =c(95,15),lhei=c(0.7,5),lwid =c(2,15),col=my_palette,cexRow=3,cexCol = 1.5,keysize = 1,Rowv = TRUE, na.rm = TRUE,
          ColSideColors=Group,scale="row",key.title=NA,labCol = FALSE)
dev.off()

library(NOISeq)
setwd("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Correlations")
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Coding_Count_Hg38.txt",header = T,row.names = 1)
head(data)
dataX=as.matrix(data)
norm=tmm(dataX)
head(norm)
write.table(norm,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/TMM.txt",sep="\t",col.names = NA,quote = FALSE)

YY=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NFE2L2.txt",header = T)
head(YY)
library(ggpubr)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NFE2L2.pdf")
ggboxplot(YY, x = "Group", y = "NFE2L2",fill="Group",order=c("ART","EC","HC","VP"),outlier.size = 0.8,add=c("mean"),
          add.params = list(color = "white",size=0.2),ggtheme = theme_gray())+
  scale_fill_manual(values=c("ART"="#FFAA99","EC"="#4682b4","HC"="#808080","VP"="#e5b647"))+labs(y="TMM normalized read count")+
  theme(axis.title.x = element_blank(),plot.margin = margin(3,0.5,3,0.5, "cm"),
        legend.position = "none",axis.text = element_text(colour = "black"),axis.title.y  = element_text(colour = "black"))
dev.off()



ip=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/ART_Specific/Network/AL161785.1/PWY/Input.txt",header = T)
head(ip)
ip$Term <- factor(ip$Term, levels = ip$Term)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/ART_Specific/Network/AL161785.1/PWY/Pathway.pdf")
ggplot(ip, aes(y=Term)) + 
  geom_point(data=ip,aes(x=1,y=Term,size=Gene,color=pval))+scale_x_discrete(limits=c("1"))+scale_color_gradient(high="#99e1ff",low="#006c99",breaks=c(0.000059093,0.004,0.01,0.05),limits=c(0.000059093,0.05))+
  scale_y_discrete(position = "right")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=9,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(3,5,4,5.5, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),legend.position = c(4.5, -0.2),legend.box="vertical",
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "# Genes",nrow = 1),color=guide_legend(title = "Adj.Pvalue",nrow=2))
dev.off()


ip=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/ART_Specific/Network/EIF2AK3-DT/Input.txt",header = T)
head(ip)
ip$Term <- factor(ip$Term, levels = ip$Term)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/ART_Specific/Network/EIF2AK3-DT/Pathway.pdf")
ggplot(ip, aes(y=Term)) + 
  geom_point(data=ip,aes(x=1,y=Term,size=Gene,color=pval))+scale_x_discrete(limits=c("1"))+scale_color_gradient(high="#99e1ff",low="#006c99",limits=c(0.1,0.6))+
  scale_y_discrete(position = "right")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=9,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(3,4,4,5.5, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),legend.position = c(4.5, -0.2),legend.box="vertical",
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "# Genes",nrow = 1),color=guide_legend(title = "Adj.Pvalue",nrow=1))
dev.off()

ip=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/ART_Specific/Network/LINC00893/Input.txt",header = T)
head(ip)
ip$Term <- factor(ip$Term, levels = ip$Term)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/ART_Specific/Network/LINC00893/Pathway.pdf")
ggplot(ip, aes(y=Term)) + 
  geom_point(data=ip,aes(x=1,y=Term,size=Gene,color=pval))+scale_x_discrete(limits=c("1"))+scale_color_gradient(high="#99e1ff",low="#006c99",breaks=c(0.03,0.1,0.4),limits=c(0.03,0.35))+
  scale_y_discrete(position = "right")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=9,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(3,4,4,7.7, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),legend.position = c(4.5, -0.2),legend.box="vertical",
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "# Genes",nrow = 1),color=guide_legend(title = "Adj.Pvalue",nrow=1))
dev.off()



dat=read.delim("/home/anoop/Desktop/tmp/Bai/Temp.txt",check.names = FALSE,header = T)
library(ggplot2)
head(dat)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Batch_Effect_Exercise/GroupWiseUMAP/EC_Umap.pdf")
ggplot(dat, aes(x=PC1, y=PC2)) + geom_point(size=4,shape=21,fill="darkgreen")+stat_ellipse(level = 0.5)
dev.off()

cls <- kmeans(x = dat[-1], centers = 3)

dat$cluster <- as.character(cls$cluster)
head(dat)
ggplot() +geom_point(data = dat, mapping = aes(x = PC1, y = PC2, colour = cluster))

write.table(dat,file="/home/anoop/Desktop/tmp/Bai/Cluster.txt",sep="\t",col.names = NA,quote = FALSE)

library(ggplot2)

ip=read.delim("/home/anoop/Desktop/tmp/BE.txt",header = T)
head(ip)
ip$Term <- factor(ip$Name, levels = ip$Name)
pdf("/home/anoop/Desktop/tmp/BE.pdf")
ggplot(ip, aes(y=Name)) + 
  geom_point(data=ip,aes(x=1,y=Term,size=-BE,color=abs(BE)))+scale_x_discrete(limits=c("1"))+scale_color_gradient(low="#99e1ff",high="#006c99")+
  scale_y_discrete(position = "right")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=9,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(1,7.3,6.5,6.5, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),legend.position = c(3.5, -0.1),legend.box="vertical",
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Binding Energy",nrow = 1,reverse=TRUE),color=FALSE)+
  scale_size("New legend",breaks=c(5,6,7,10,12),labels=c(-1,-6,-8,-10,-12))
dev.off()

XX=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/Uniq/PWY/Input.txt",header = T)
library(ggplot2)
library(ggalluvial)
head(XX)
col=c(rep("black",24),rep("#92bacf",28))
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/Uniq/PWY/Pathway.pdf",height = 10,width = 5.5)
ggplot(XX,
       aes(axis1 = Pathway, axis2 = Genes)) +
  geom_alluvium(aes(fill=Pathway))+scale_fill_manual(values = c("#92bacf","#cbaeb8"))+
  geom_stratum(fill = c("#975d72","#2675a0","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2"),width = 0.2) +
  geom_label(stat = "stratum", infer.label = TRUE,size=3) +
  scale_x_discrete(limits = c("Pathway", "Genes"), expand = c(.05, .05)) +
  theme(axis.text.y = element_blank(),axis.ticks.y=element_blank())+theme(legend.position = "none",axis.text.x = element_text(size=12))
dev.off()

ip=read.delim("/home/anoop/Desktop/tmp/Bai/Figure/Pvalue.txt",header = T)
head(ip)
pdf("/home/anoop/Desktop/tmp/Bai/Figure/Pvalue.pdf",height = 12,width = 5)
ggplot(ip, aes(y=Gene)) + 
  geom_point(data=ip,aes(x=1,y=Gene,size=SIZE,color=pval))+scale_x_discrete(limits=c("1"))+scale_color_gradient(high="#99e1ff",low="#006c99")+
  scale_y_discrete(position = "right")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=6,color="black"),
        axis.ticks = element_blank(),plot.margin = margin(0.1,2,0.1,1, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),legend.box="vertical",
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=FALSE,color=guide_legend(title = "Adjusted Pvalue",nrow=2))+scale_size(range = c(1, 3))
dev.off()
library(countToFPKM)
cnt=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/rest.txt",row.names = 1)
head(cnt,2)
anno=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/Length.txt")
head(anno)
sam=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/FragLen.txt")
head(sam)
meanFragmentLength <- sam$meanFragmentLength
featureLength <- anno$length
fpkm_matrix <- fpkm (as.matrix(cnt), featureLength, meanFragmentLength)
?fpkm
write.table(fpkm_matrix,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/ART_Specific/FPKM.txt",sep="\t",quote = FALSE,col.names = NA)


library(ggplot2)
library(reshape2)
library(plyr)
library(scales)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/log2cpm/LCPM.txt",header = T)
data1=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/log2cpm/OxphoseTMM.txt",header = T)
head(data)
dat=log2(data)
m=melt(data)
m$log=log2(m$value)
?log()
head(m)
m <- ddply(m, .(variable), transform,rescale = scale(value))
ggplot(m, aes(variable, Name)) + geom_tile(aes(fill = rescale),color="white") + scale_fill_gradient(low = "green",high = "red")+theme(axis.text.x = element_text(angle = 45))

sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/sampleInfo.txt",row.names = 1,header = T)
geneinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/log2cpm/GeneInfo.txt",row.names = 1,header = T)
Group <- c("#FFA500","#4682b4","#96a94e","#ff4c4c")[sampleinfo$Group]
tail(geneinfo)
Gene<-c("#a3dd9e","#dcdd9e","#ddbb9e","#dd9ec0","white")[geneinfo$Complex]
ncol(dat)
library(gplots)
my_palette <- colorRampPalette(c("#004000","#008000","#fde6ed","#ea1313","#8c0b0b"))(n=30)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/HM.pdf",width = 12,height = 17)
heatmap.2(as.matrix(dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 0.5,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
          key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(0.7,8.5),margins  =c(5,10),ColSideColors=Group,RowSideColors = Gene)
dev.off()

library(gplots)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/log2cpm/LCPM.txt",header = T,row.names = 1,check.names = FALSE)
geneinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/log2cpm/Gene.txt",row.names = 1,header = T)
Gene<-c("#16605f","#891b69","#e5ba33","#332a6c")[geneinfo$complex]

my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/log2cpm/HeatMap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(data),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="col",
          key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),ColSideColors=Gene,RowSideColors = Group,labRow = FALSE)
dev.off()

head(X$carpet)

input=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/Diff/LFC_Trans.txt",header = T,check.names = FALSE,row.names = 1)

not=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/Diff/pval_trans.txt",header = T,check.names = FALSE,row.names = 1)
geneinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/Diff/Gene.txt",row.names = 1,header = T)
Gene<-c("#16605f","#891b69","#e5ba33","#332a6c")[geneinfo$complex]
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/Diff/HeatMap.pdf",width =20,height=5)
heatmap.2(as.matrix(input),trace="none",cellnote = not,notecol="#4c4c4c",notecex=2,
          col=my_palette,cexRow=1.5,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="none",
          key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,4),margins  =c(9,10),lwid=c(1,11),ColSideColors=Gene)
dev.off()


pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/Diff/Test.pdf",width = 20,height=15)
ggarrange((heatmap.2(as.matrix(data),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="col",
                     key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),ColSideColors=Gene,RowSideColors = Group,labRow = FALSE)),
          (heatmap.2(as.matrix(input),trace="none",cellnote = not,notecol="#4c4c4c",notecex=2,
                     col=my_palette,cexRow=1.5,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="none",
                     key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,4),margins  =c(9,10),lwid=c(1,11),ColSideColors=Gene)),nrow = 2,ncol = 1)
dev.off()
library(ComplexHeatmap)
m1=Heatmap(as.matrix(input))
m2=Heatmap(as.matrix(data))
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/Diff/Test.pdf",width = 20,height=15)
m2 %v% m1
dev.off()
packageVersion("complexheatmap")
## ART Specific non-coding
library(edgeR)
count=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/NonCoding_Count_Hg38.txt",row.names = 1, check.names = FALSE)
countX=as.matrix(count)
meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Batch_Effect_Exercise/MetaData.txt",row.names = 1)
expSet=ExpressionSet(countX, phenoData=AnnotatedDataFrame(data=meta))
expSet <- expSet[rowSums(exprs(expSet)) != 0, ] # to remove genes having zeros in all samples
log2cpm <- cpm(exprs(expSet), log = TRUE)

write.table(log2cpm,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/log2CPM.txt",sep="\t",quote = FALSE,col.names = NA)

library(sva)
BiocManager::install("sva")
n
library(genefilter)
prot=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/log2CPM.txt",row.names = 1)
mat <- as.matrix(prot)
des=read.table("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Batch_Effect_Exercise/MetaData.txt",sep="\t",header=T)
designCombat = model.matrix(~ des$group)
rnaseqCombat = ComBat(mat, batch = des$batch, mod = designCombat, par.prior = TRUE, prior.plots = TRUE)
write.table(rnaseqCombat,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/Combat_Result.txt", sep="\t",quote=FALSE,col.names = NA)

library(umap)
data=read.table("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/Combat_ARTs.txt",sep="\t",header=T,row.names = 1)
dataX=t(data)
rnames <- data[,1] 
mat_data <- data.frame(dataX[,2:ncol(dataX)])
Art.umap = umap(mat_data)
head(Art.umap$layout,3)
write.table(Art.umap$layout,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/ART_Umap.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/ART_Umap.txt",header = T)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/ART_Umap.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=4,shape=21,aes(fill=Group))+
  scale_color_manual(values=c(VP="#b20000",HC="#0d3a1b",EC="#315b7d",ART="#b27300"))+
  scale_fill_manual(values=c(VP="#ff4c4c",HC="#96a94e",EC="#4682b4",ART="#FFA500"))+
  theme(axis.title = element_text(size=9,hjust = 0.5),plot.margin = margin(0.7,0.5,0.7,0.5, "cm"),
        legend.position = c(0.06, 0.12))+labs(x = "UMAP1",y="UMAP2")#+
#geom_text_repel(data=dat,aes(x=V1,y=V2,label = ID),size=1.9,box.padding=0.3,show.legend=FALSE,colour="black")
dev.off()

DAT=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Correlations/New/Regression/Signi.txt",header = T,check.names = FALSE,row.names = 1)
head(DAT)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Correlations/New/Regression/Test.pdf")
p1=ggplot(DAT, aes(x=DAT$`Duration of Suppressive ART`, y=DAT$NDUFB5)) + geom_point(color="#3f4b3b",size=2.5) +
  stat_smooth(method = lm)+labs(x="Duration of Suppressive ART",y="NDUFB5")+theme_bw()+theme(axis.title = element_text(color = "black"))

p2=ggplot(DAT, aes(x=DAT$`Duration of Suppressive ART`, y=DAT$NDUFA13)) + geom_point(color="#3f4b3b",size=2.5) +
  stat_smooth(method = lm)+labs(x="Duration of Suppressive ART",y="NDUFA13")+theme_bw()+theme(axis.title = element_text(color = "black"))

p3=ggplot(DAT, aes(x=DAT$`Duration of ART`, y=DAT$NDUFB1)) + geom_point(color="#3f4b3b",size=2.5) +
  stat_smooth(method = lm)+labs(x="Duration of ART",y="NDUFB1")+theme_bw()+theme(axis.title = element_text(color = "black"))
p3
library(ggpubr)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Correlations/New/Regression/Test.pdf",width = 15,height=5)
ggarrange(p1,p2,p3,nrow = 1,ncol = 3)
dev.off()


############################ Complex HeatMap

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/log2cpm/LCPM.txt",header = T,row.names = 1,check.names = FALSE)
head(data)

my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/log2cpm/HeatMap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(data),tracecol=NA,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="col",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

head(X$carpet)

sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/sampleInfo.txt",row.names = 1,header = T)
geneinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/Diff/GeneInfo.txt",header = T,row.names = 1)

input=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/Diff/LFC_Trans.txt",header = T,check.names = FALSE,row.names = 1)
write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/log2cpm/ZScore.txt",col.names = NA,quote = FALSE,sep = "\t")

DAT=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/log2cpm/ZScore.txt",row.names = 1,header = T,check.names = FALSE)
head(DAT)
library(ComplexHeatmap)
row.names(DAT) <- factor(row.names(DAT), levels = row.names(DAT))
colours <- list("Cohort"=c("HC"="#96a94e","VP"="#ff4c4c","EC"="#4682b4","ART"="#FFA500"))
comp<- list("Complex"=c("Complex1"="#16605f","Complex2"="#891b69","Complex3"="#e5ba33","Complex4"="#332a6c","Complex5"="#d17294"))
H1=Heatmap(as.matrix(DAT),col=col_fun1,cluster_rows=FALSE,cluster_columns = FALSE,left_annotation = rowAnnotation(Cohort = sampleinfo$Group,col=colours,show_legend=FALSE,show_annotation_name=FALSE),
           top_annotation =columnAnnotation(Complex = geneinfo$Complex,col=comp,show_legend=FALSE,show_annotation_name=FALSE),
           name = "Z-Score",show_row_names = FALSE,row_names_gp=gpar(fontsize = 5),row_order = order(rownames(DAT)),
           row_split=c(rep("HC",19),rep("VP",19),rep("EC",19),rep("ART",19)),column_split = geneinfo$Complex,height=10)
H2=Heatmap(as.matrix((input)),col=col_fun,cluster_rows=FALSE,cluster_columns = FALSE,name="Log2FoldChange",row_names_gp=gpar(fontsize = 12),height=3,column_names_gp =gpar(fontsize = 8)) 

tt=H1 %v% H2
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/Diff/Test.pdf",width = 20,height=7)
draw(tt, merge_legend = TRUE)
dev.off()

col_fun1 = colorRamp2(c(5, 2, 0, -2, -5), c("#0000cc","#4c4cff" ,"white","#e5e500", "#b2b200"))
col_fun = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))
library(circlize)



XX=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/NetWorkAnalysis/From_rui/network analysis/results/global/association_analysis/Comm_0001/results/global/C5/PWY/Input.txt",header = T)
library(ggplot2)
library(ggalluvial)
head(XX)
col=c(rep("black",24),rep("#92bacf",28))
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/NetWorkAnalysis/From_rui/network analysis/results/global/association_analysis/Comm_0001/results/global/C5/PWY/Pathway.pdf",height = 10,width = 5.5)
ggplot(XX,
       aes(axis1 = Pathway, axis2 = Genes)) +
  geom_alluvium(aes(fill=Pathway))+scale_fill_manual(values = c("#92bacf","#cbaeb8"))+
  geom_stratum(fill = c("#975d72","#2675a0","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2"),width = 0.2) +
  geom_label(stat = "stratum", infer.label = TRUE,size=3) +
  scale_x_discrete(limits = c("Pathway", "Genes"), expand = c(.05, .05)) +
  theme(axis.text.y = element_blank(),axis.ticks.y=element_blank())+theme(legend.position = "none",axis.text.x = element_text(size=12))
dev.off()
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/NetWorkAnalysis/From_rui/network analysis/results/global/association_analysis/Comm_0001/results/global/C5/PWY/Pathway.pdf",height = 10,width = 5.5)

ggplot(XX,
       aes(axis1 = Pathway, axis2 = Genes)) +
  geom_alluvium(aes(fill=Pathway))+scale_fill_manual(values = c("#92bacf","#cbaeb8"))+
  geom_stratum(fill = c("#975d72","#2675a0","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2","#b2b2b2"),width = 0.2) +
  geom_label(stat = "stratum", infer.label = TRUE,size=3) +
  scale_x_discrete(limits = c("Pathway", "Genes"), expand = c(.05, .05)) +
  theme(axis.text.y = element_blank(),axis.ticks.y=element_blank())+theme(legend.position = "none",axis.text.x = element_text(size=12))
dev.off()




YY=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/For_boxplot.txt",header = T)
head(YY)
library(ggpubr)
library(ggplot2)
library(reshape2)
M=melt(YY)
head(M)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/Boxplot.pdf")
ggboxplot(M, x = "ID", y = "value", facet.by = "variable",fill="ID",order=c("ART","EC","HC","VP"),outlier.shape =NA,ggtheme = theme_gray())+
  scale_fill_manual(values=c(VP="#ff4c4c",HC="#96a94e",EC="#4682b4",ART="#FFA500"))+
  geom_jitter(color="black",shape=16,size=0.7, position=position_jitter(0.08))+labs(y="log2CPM")+
  theme(axis.title.x = element_blank(),plot.margin = margin(0.7,0.5,10,0.5, "cm"),
        legend.position = "none",axis.text = element_text(colour = "black"),axis.title.y  = element_text(colour = "black"))
dev.off()

X1=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/AL161785.1.txt",header = T)
head(X1)
p1=ggboxplot(X1, x = "ID", y = "AL161785.1",fill="ID",order=c("ART","EC","HC","VP"),outlier.shape =NA,ggtheme = theme_gray())+
  scale_fill_manual(values=c(VP="#ff4c4c",HC="#96a94e",EC="#4682b4",ART="#FFA500"))+
  geom_jitter(color="black",shape=16,size=0.7, position=position_jitter(0.08))+labs(y="FPKM",title = "AL161785.1")+
  theme(axis.title.x = element_blank(),plot.margin = margin(0.7,0.5,10,0.5, "cm"),plot.title = element_text(hjust = 0.5,size = 11),
        legend.position = "none",axis.text = element_text(colour = "black"),axis.title.y  = element_text(colour = "black",size = 10))

X2=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/EIF2AK3-DT.txt",header = T)
head(X2)
p2=ggboxplot(X2, x = "ID", y = "EIF2AK3.DT",fill="ID",order=c("ART","EC","HC","VP"),outlier.shape =NA,ggtheme = theme_gray())+
  scale_fill_manual(values=c(VP="#ff4c4c",HC="#96a94e",EC="#4682b4",ART="#FFA500"))+
  geom_jitter(color="black",shape=16,size=0.7, position=position_jitter(0.08))+labs(y="FPKM",title = "EIF2AK3-DT")+
  theme(axis.title.x = element_blank(),plot.margin = margin(0.7,0.5,10,0.5, "cm"),plot.title = element_text(hjust = 0.5,size = 11),
        legend.position = "none",axis.text = element_text(colour = "black"),axis.title.y  = element_blank())

X3=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/LINC00893.txt",header = T)
p3=ggboxplot(X3, x = "ID", y = "LINC00893",fill="ID",order=c("ART","EC","HC","VP"),outlier.shape =NA,ggtheme = theme_gray())+
  scale_fill_manual(values=c(VP="#ff4c4c",HC="#96a94e",EC="#4682b4",ART="#FFA500"))+
  geom_jitter(color="black",shape=16,size=0.7, position=position_jitter(0.08))+labs(y="FPKM",title = "LINC00893")+
  theme(axis.title.x = element_blank(),plot.margin = margin(0.7,0.5,10,0.5, "cm"),plot.title = element_text(hjust = 0.5,size = 11),
        legend.position = "none",axis.text = element_text(colour = "black"),axis.title.y  = element_blank())

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/Boxplot.pdf")
ggarrange(p1,p2,p3,nrow = 1,ncol = 3)
dev.off()


cnt=read.csv("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/VP.txt",row.names = 1)
head(cnt,2)
anno=read.csv("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/GeneLength.txt")
head(anno)
sam=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/VP_len.txt")
head(sam)
meanFragmentLength <- sam$meanFragmentLength
featureLength <- anno$length
fpkm_matrix <- fpkm (as.matrix(cnt), featureLength, meanFragmentLength)
head(fpkm_matrix)
write.table(fpkm_matrix,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/VP_FPKM.txt",sep="\t",quote = FALSE,col.names = NA)
?fpkm



data=read.csv("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/Volcano_Result_1.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/Volcano1.pdf")
ggplot(data, aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(data=subset(data, padj<.05 & log2FoldChange <= -1.5),aes(x=log2FoldChange,y=-log10(padj)),color="#004000",size=1.2)+
  geom_point(data=subset(data, padj<.05 & log2FoldChange >= 1),aes(x=log2FoldChange,y=-log10(padj)),color="#b20000",size=1.2)+
  geom_point(data=subset(data, padj<.05 & log2FoldChange > 0 & log2FoldChange < 1.5),aes(x=log2FoldChange,y=-log10(padj)),color="#b2b200",size=1.2)+
  geom_point(data=subset(data, padj<.05 & log2FoldChange < 0 & log2FoldChange > -1.5),aes(x=log2FoldChange,y=-log10(padj)),color="#b2b200",size=1.2)+
  geom_point(data=subset(data, padj>=.05),aes(x=log2FoldChange,y=-log10(padj)),color="#90b4d2",size=1.2)+
  geom_vline(xintercept=1.5, linetype="dashed",size=0.35)+
  geom_vline(xintercept=-1.5, linetype="dashed",size=0.35)+scale_x_continuous(limits = c(-8, 8), breaks = seq(-8, 8, by = 2))+
  geom_hline(yintercept=1.3010299957, linetype="dashed",size=0.35)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
                   plot.title = element_text(hjust = 0.5,size =9),panel.grid.minor = element_blank(),
                   axis.title=element_text(size=18),axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),plot.margin = margin(0.9,1,0.9,1, "cm"))+
  labs(x="Log2 Fold Change",y="-log10 (Adj.pvalue)")
dev.off()

data=read.csv("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/HC-ART/Volcano_Result.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/HC-ART/Volcano1.pdf")
ggplot(data, aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(data=subset(data, padj<.05 & log2FoldChange <= -1.5),aes(x=log2FoldChange,y=-log10(padj)),color="#004000",size=1.2)+
  geom_point(data=subset(data, padj<.05 & log2FoldChange >= 1),aes(x=log2FoldChange,y=-log10(padj)),color="#b20000",size=1.2)+
  geom_point(data=subset(data, padj<.05 & log2FoldChange > 0 & log2FoldChange < 1.5),aes(x=log2FoldChange,y=-log10(padj)),color="#b2b200",size=1.2)+
  geom_point(data=subset(data, padj<.05 & log2FoldChange < 0 & log2FoldChange > -1.5),aes(x=log2FoldChange,y=-log10(padj)),color="#b2b200",size=1.2)+
  geom_point(data=subset(data, padj>=.05),aes(x=log2FoldChange,y=-log10(padj)),color="#90b4d2",size=1.2)+
  geom_vline(xintercept=1.5, linetype="dashed",size=0.35)+
  geom_vline(xintercept=-1.5, linetype="dashed",size=0.35)+scale_x_continuous(limits = c(-8, 8), breaks = seq(-8, 8, by = 2))+
  geom_hline(yintercept=1.3010299957, linetype="dashed",size=0.35)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),panel.grid.minor = element_blank(),
        axis.title=element_text(size=18),axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),plot.margin = margin(0.9,1,0.9,1, "cm"))+
  labs(x="Log2 Fold Change",y="-log10 (Adj.pvalue)")
dev.off()


Sam=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/sam.txt",header = T)
library(ComplexHeatmap)
col_fun1 = colorRamp2(c(5, 2, 0, -2, -5), c("#0000cc","#4c4cff" ,"white","#e5e500", "#b2b200"))

colours <- list("Cohort"=c("EC"="#4682b4","ART"="#FFA500"))
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/HeatMap.pdf",height = 10,width= 20)
Heatmap((X$carpet),col=col_fun1,cluster_rows=FALSE,cluster_columns = TRUE,show_column_names = FALSE,column_dend_height  = unit(4, "cm"),height  = unit(12, "cm"),
           left_annotation  =rowAnnotation(Cohort = Sam$group,col=colours,show_legend=FALSE,show_annotation_name=FALSE),column_split = 3,width = unit(40, "cm"),
           name = "Z-Score",show_row_names = FALSE,row_names_gp=gpar(fontsize = 5),
        row_split=c(rep("ART",19),rep("EC",19)))

dev.off()
PWy=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/PWY/G5_1.txt",header = T, check.names = FALSE,row.names = 1)
head(PWy)
col = colorRamp2(c(0, 1), c("grey","#b03e69"))
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/PWY/G1.pdf")
H2=Heatmap(t(PWy),col=col,cluster_rows=FALSE,cluster_columns = FALSE,show_heatmap_legend = FALSE,height = unit(3, "cm"),width = unit(50, "cm"),show_column_names = FALSE)
dev.off()
ncol(X$carpet)
ncol(t(PWy))
tt=H1 %v% H2
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/HeatMap.pdf",height = 8,width= 25)
draw(tt, merge_legend = TRUE)
dev.off()

dend = hclust(dist(t(X$carpet)))
C=cutree(dend, k = 3)
D=as.data.frame(C)
write.table(D,file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/Clustr.txt",sep = "\t",quote = FALSE,col.names = NA)
library(dendextend)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/Names_HeatMap.txt",header = T,row.names = 1,check.names = FALSE)
head(data)
library(gplots)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/HeatMap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(data),tracecol=NA,cexRow=0.8,col=my_palette,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()




Sam=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/HeatMap/grp.txt",header = T)
library(ComplexHeatmap)
col_fun1 = colorRamp2(c(5, 2, 0, -2, -5), c("#0000cc","#4c4cff" ,"white","#e5e500", "#b2b200"))

colours <- list("Cohort"=c("HC"="#96a94e","VP"="#ff4c4c","EC"="#4682b4","ART"="#FFA500"))
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/HeatMap/HeatMap.pdf",height = 10,width= 20)
Heatmap((X$carpet),col=col_fun1,cluster_rows=FALSE,cluster_columns = TRUE,show_column_names = TRUE,height  = unit(10, "cm"),column_split = 2,
        left_annotation  =rowAnnotation(Cohort = Sam$Group,col=colours,show_legend=FALSE,show_annotation_name=FALSE),width = unit(40, "cm"),
        name = "Z-Score",show_row_names = FALSE,row_names_gp=gpar(fontsize = 5),
        row_split=c(rep("ART",19),rep("EC",19),rep("HC",19),rep("VP",19)))

dev.off()
data=read.csv("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/HeatMap/ART_lcpm.txt",row.names = 1,header = T)

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/HeatMap/HeatMap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(data),tracecol=NA,cexRow=0.8,col=my_palette,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()



Sam=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/HeatMap/grp.txt",header = T)
library(ComplexHeatmap)
col_fun1 = colorRamp2(c(5, 2, 0, -2, -5), c("#8c0b0b","#ea1313" ,"white","#008000","#004000"))

colours <- list("Cohort"=c("HC"="#96a94e","VP"="#ff4c4c","EC"="#4682b4","ART"="#FFA500"))
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/HeatMap/HeatMapCmb.pdf",height = 10,width= 20)
Heatmap((X$carpet),col=col_fun1,cluster_rows=FALSE,cluster_columns = TRUE,show_column_names = TRUE,height  = unit(10, "cm"),column_split = 2,column_dend_height = unit(2, "cm"),
        left_annotation  =rowAnnotation(Cohort = Sam$Group,col=colours,show_legend=FALSE,show_annotation_name=FALSE),width = unit(40, "cm"),
        name = "Z-Score",show_row_names = FALSE,row_names_gp=gpar(fontsize = 5),
        row_split=c(rep("ART",19),rep("EC",19),rep("HC",19),rep("VP",19)))

dev.off()
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/HeatMap/ART_Combat.txt",row.names = 1,header = T)

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NonCoding/HeatMap/HeatMapcmb.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(data),tracecol=NA,cexRow=0.8,col=my_palette,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()





#################################

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/BatchCorrected/Oxp_Combat.txt",header = T,row.names = 1,check.names = FALSE)
head(data)

my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/BatchCorrected/HeatMap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(data),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/BatchCorrected/Zscore.txt",sep="\t",col.names = NA,quote = FALSE)
###################
library(circlize)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/sampleInfo.txt",row.names = 1,header = T)
geneinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/BatchCorrected/GeneInfo_Sorted.txt",header = T,row.names = 1)

input=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/BatchCorrected/LFC.txt",header = T,check.names = FALSE,row.names = 1)

DAT=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/BatchCorrected/Combat_Sorted.txt",check.names = FALSE,header = T,row.names = 1)
col_fun = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))
library(ComplexHeatmap)
nrow(t(X$carpet))
colours <- list("Cohort"=c("ART"="#FFA500","EC"="#4682b4","HC"="#96a94e","VP"="#ff4c4c"))
comp<- list("Complex"=c("Complex1"="#16605f","Complex2"="#891b69","Complex3"="#e5ba33","Complex4"="#332a6c","Complex5"="#d17294"))
col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7F7F00","#B2B200" ,"#E5E500","white","#BF7FBF","#993299","#590059"))
H1=Heatmap(as.matrix(t(DAT)),col=col_fun1,cluster_rows=FALSE,cluster_columns = FALSE,left_annotation = rowAnnotation(Cohort = sampleinfo$Group,col=colours,show_legend=FALSE,show_annotation_name=FALSE),
           top_annotation =columnAnnotation(Complex = geneinfo$Complex,col=comp,show_legend=FALSE,show_annotation_name=FALSE),
           name = "Z-Score",show_row_names = FALSE,row_names_gp=gpar(fontsize = 5),row_order = order(rownames(t(DAT))),
           row_split=c(rep("ART",19),rep("EC",19),rep("HC",19),rep("VP",19)),column_split = geneinfo$Complex,height=10)

H2=Heatmap(as.matrix((input)),col=col_fun,cluster_rows=FALSE,cluster_columns = FALSE,name="Log2FoldChange",row_names_gp=gpar(fontsize = 12),height=3,column_names_gp =gpar(fontsize = 8)) 

colnames(t(DAT))
colnames(input)
tt=H1 %v% H2
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/BatchCorrected/HeatMap.pdf",width = 20,height=7)
draw(tt, merge_legend = TRUE)
dev.off()

col_fun1 = colorRamp2(c(5, 2, 0, -2, -5), c("#0000cc","#4c4cff" ,"white","#e5e500", "#b2b200"))



################## EC-ART

library(ComplexHeatmap)
library(gplots)

dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/HeatMap/Names_HeatMap.txt",check.names = FALSE,row.names = 1,header = T)

my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/HeatMap/HeatMap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()
library(circlize)
col_fun1 = colorRamp2(c(5, 2, 0, -2, -5), c("#23415a","#3f75a2" ,"white","#d67834", "#7e3f12"))
sample=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/HeatMap/sample.txt",header = T)
colours <- list("Cohort"=c("ART"="#FFA500","EC"="#4682b4"))
H1=Heatmap(as.matrix(t(X$carpet)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,top_annotation = columnAnnotation(Cohort = sample$Group,col=colours,show_legend=FALSE,show_annotation_name=FALSE),
           name = "Z-Score",show_row_names = FALSE,row_names_gp=gpar(fontsize = 5),row_dend_width = unit(4, "cm"),row_gap = unit(2, "mm"),
           column_split=c(rep("ART",19),rep("EC",19)), width = unit(15, "cm"),show_column_names = FALSE,height  = unit(35, "cm"),
           right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 2:4,fill=c("#c4e5d6","#fff1d6","#ecd1fc")),
                                                             labels = c("Oxidative phosphorylation\nThermogenesis", "No significant enrichment", "Oxidative phosphorylation\nThermogenesis"), 
                                                             labels_gp = gpar(col = "blue", fontsize = 18))),row_split  = 3)

?anno_block
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/HeatMap/HeatMap.pdf",width = 15,height = 20)
H1
dev.off()


library(circlize)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/HeatMap/Circos/Chord_Input.txt")
head(data)
col_fun = colorRamp2(c(-2.5,0,1.5), c("#004c00", "white","#ae0001"), transparency = 0)

grid=c(HSPA1B = "#A67C00",H2AFX = "#A67C00",RPL17 = "#A67C00",TUBB4B = "#A67C00",ATP13A2 = "#A67C00",RAB5IF = "#A67C00",NEK7 = "#A67C00",
           EP400 = "#A67C00",BICD2 = "#A67C00",RBM14 = "#A67C00",TBL1XR1 = "#A67C00",SENP3 = "#A67C00",XPO1 = "#A67C00",EDEM3 = "#A67C00",
           SLC30A7 = "#A67C00",FMR1 = "#A67C00",GAK = "#A67C00",PAK2 = "#A67C00",CPSF2 = "#A67C00",MYH9 = "#A67C00",HERC2 = "#A67C00",
           DYNC1H1 = "#A67C00",PDCD6IP = "#A67C00",UBR4 = "#A67C00",ARID1A = "#A67C00",AP1M1 = "#A67C00",MLLT1 = "#A67C00",
           TRRAP = "#A67C00",BRD4 = "#A67C00",SMARCA4 = "#A67C00",CAPN1 = "#A67C00",BRD3 = "#A67C00",XRCC6 = "#A67C00",POLR2G = "#A67C00",
           KARS = "#A67C00",DDX19B = "#A67C00",HNRNPC = "#A67C00",DAP3 = "#A67C00",EIF3I = "#A67C00",YBX1 = "#A67C00",H2AFV = "#A67C00",
           GNL3 = "#A67C00",SSBP1 = "#A67C00",B2M = "#A67C00",POLR2D = "#A67C00",PSMA7 = "#A67C00",POLE3 = "#A67C00",MCM7 = "#A67C00",
           POLR2C = "#A67C00",EIF3H = "#A67C00",SAP18 = "#A67C00",HNRNPA1 = "#A67C00",PSMB6 = "#A67C00",ANP32B = "#A67C00",DHX9 = "#A67C00",
           EIF3F = "#A67C00",CCDC32 = "#A67C00",EIF3E = "#A67C00",PSMA3 = "#A67C00",BTF3 = "#A67C00",PSMB7 = "#A67C00",AIMP1 = "#A67C00",
           EIF3K = "#A67C00",GTF2F2 = "#A67C00",PSMB4 = "#A67C00",PPIB = "#A67C00",FAM133B = "#A67C00",NPM1 = "#A67C00",RNF7 = "#A67C00",
           FRG1 = "#A67C00",VDAC2 = "#A67C00",PRDX4 = "#A67C00",EIF3D = "#A67C00",RAN = "#A67C00",ATRAID = "#A67C00",PSPC1 = "#A67C00",
           NOP58 = "#A67C00",EMG1 = "#A67C00",GTPBP4 = "#A67C00",NUP88 = "#A67C00",C1QBP = "#A67C00",CD63 = "#A67C00",FUT8 = "#A67C00",
           COX5A = "#A67C00",PSMC1 = "#A67C00",CDK7 = "#A67C00",BRIX1 = "#A67C00",RBX1 = "#A67C00",CD59 = "#A67C00",EEF1A1 = "#A67C00",
           ILF2 = "#A67C00",RTRAF = "#A67C00",RANBP6 = "#A67C00",PPIA = "#A67C00",DCAF16 = "#A67C00",CCNH = "#A67C00",RPL6 = "#A67C00",
           SETMAR = "#A67C00",RPS3A = "#A67C00",RPS4X = "#A67C00",PSMB1 = "#A67C00",COX4I1 = "#A67C00",SNRPF = "#A67C00",ERH = "#A67C00",
           RPL27A = "#A67C00",SNRPE = "#A67C00",COX6B1 = "#A67C00",RPS3 = "#A67C00",RPS6 = "#A67C00",RPL11 = "#A67C00",DDX24 = "#A67C00",
           RPSA = "#A67C00",RPS16 = "#A67C00",RPL35 = "#A67C00",RPS23 = "#A67C00",RPL32 = "#A67C00",
       SNRPG ="#A67C00",RPS25 = "#A67C00",RPS27A = "#A67C00",RPL34 = "#A67C00",TOE1 = "#A67C00",RPL23 = "#A67C00",RPL27 = "#A67C00",GTF2B = "#A67C00",
       RPL37A = "#A67C00",RPL35A = "#A67C00",Tat = "#2f80ed",Rev = "#2f80ed",Nef = "#2f80ed",Env = "#2f80ed",Gag_Pol = "#2f80ed",
       Protease = "#2f80ed",Vif ="#2f80ed" , Integrase = "#2f80ed",Gag = "#2f80ed",p6 = "#2f80ed",Vpu = "#2f80ed",p24 = "#2f80ed")

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/HeatMap/Circos/test.pdf")
chordDiagram(data, annotationTrack = "grid",col=col_fun,grid.col=grid,
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(data))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex=0.4)
}, bg.border = NA)
dev.off()

##################################


library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/GeneSets.txt",header = T,row.names = 1)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/HeatMap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/Zscore.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/Zscore.txt",header = T, row.names = 1)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/samples.txt",row.names = 1,header = T)

input=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/Out3.txt",header = T,check.names = FALSE,row.names = 1)



col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#990000","#e50000" ,"#ff9999","white","#b5cde1","#315b7d","#23415a"))

col_fun= c("white","#004000")
library(ComplexHeatmap)
library(circlize)

colours <- list("Cohort"=c("ART"="#FFA500","EC"="#4682b4","HC"="#96a94e","VP"="#ff4c4c"))

H1=Heatmap(as.matrix(Zscore),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_dend_width = unit(3, "cm"),
           top_annotation =columnAnnotation(Cohort = sampleinfo$group,col=colours,show_legend=FALSE,show_annotation_name=FALSE),
           name = "Z-Score",show_row_names = FALSE,row_names_gp=gpar(fontsize = 5),height  = unit(35, "cm"),width  = unit(20, "cm"),row_split = 6,
           column_split =c(rep("ART",19),rep("EC",19),rep("HC",19),rep("VP",19)))
H2=Heatmap(as.matrix((input)),col=col_fun,cluster_rows=FALSE,cluster_columns = FALSE,name="NULL",show_row_names = FALSE,
           column_names_gp=gpar(fontsize = 15),height  = unit(35, "cm"),width  = unit(3, "cm")) 

tt=H1 + H2
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/HeatMap1.pdf",height = 20,width =15)
draw(tt, merge_legend = TRUE)
dev.off()
library(umap)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/Data.txt",header=T,row.names = 1)
head(dataX)
dataX=t(data)
rnames <- data[,1] 
mat_data <- data.frame(dataX[,2:ncol(dataX)])
head(data)
Art.umap = umap(data)
head(art.umap,3)
write.table(Art.umap$layout,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/Umap.txt",sep="\t",col.names = NA,quote = FALSE)
dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/Umap.txt",row.names = 1)
library(ggplot2)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=4,shape=21,aes(fill=Group))+
  scale_color_manual(values=c(VP="#b20000",HC="#0d3a1b",EC="#315b7d",ART="#b27300"))+
  scale_fill_manual(values=c(VP="#ff4c4c",HC="#96a94e",EC="#4682b4",ART="#FFA500"))+
  theme(axis.title = element_text(size=9,hjust = 0.5),legend.position = c(0.07, 0.12),plot.margin = margin(0.7,0.5,0.7,0.5, "cm"))+labs(x = "UMAP1",y="UMAP2")
dev.off()

###################### Maike's project

library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/GeneSets.txt",header = T,row.names = 1)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/HeatMap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/Zscore.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/Zscore.txt",header = T, row.names = 1)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/samples.txt",row.names = 1,header = T)

input=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/Out2.txt",header = T,check.names = FALSE,row.names = 1)



col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#990000","#e50000" ,"#ff9999","white","#b5cde1","#315b7d","#23415a"))

col_fun= c("#e0e0e0","orange")
library(ComplexHeatmap)
library(circlize)

colours <- list("Cohort"=c("EC"="#4682b4","HC"="#96a94e","VP"="#ff4c4c"))
LFC=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/Exp3.txt",header = T,check.names = FALSE,row.names = 1)

H1=Heatmap(as.matrix(Zscore),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_dend_width = unit(3, "cm"),
           top_annotation =columnAnnotation(Cohort = sampleinfo$Group,col=colours,show_legend=FALSE,show_annotation_name=FALSE),
           name = "Z-Score",show_row_names = FALSE,row_names_gp=gpar(fontsize = 5),height  = unit(35, "cm"),width  = unit(20, "cm"),row_split = 6,
           column_split =c(rep("EC",14),rep("HC",12),rep("VP",16)))
H2=Heatmap(as.matrix((input)),col=col_fun,cluster_rows=FALSE,cluster_columns = FALSE,name="NULL",show_row_names = FALSE,
           column_names_gp=gpar(fontsize = 15),height  = unit(35, "cm"),width  = unit(3, "cm")) 

col_fun_lfc = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))
H3=Heatmap(as.matrix((LFC)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="Log2FoldChange",width  = unit(3, "cm"),show_row_names = FALSE,
           row_names_gp=gpar(fontsize = 12),height=3,column_names_gp =gpar(fontsize = 15),na_col = "#e0e0e0") 

tt=H1 + H2 + H3
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/HeatMap1.pdf",height = 20,width =15)
draw(tt, merge_legend = TRUE)
dev.off()
BiocManager::install(c("DESeq2"))
BiocManager::valid()
s

BiocManager::install(c("DESeq2"
), update = TRUE, ask = FALSE)

unloadNamespace("DESeq2")
library(DESeq2)
count=read.csv("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/EC-VP/Input.txt",row.names = 1, check.names = FALSE)
Design=read.csv("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/EC-VP/Info.txt")
ds=DESeqDataSetFromMatrix(countData = count,colData = Design,design = ~sex+batch+group)
ds1=DESeq2::DESeq(ds)
?results
res=results(ds1,contrast = c("group","VP","EC"),independentFiltering = FALSE)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/EC-VP/Result.txt",sep="\t", quote=FALSE,col.names = NA)


count=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/HC-EC/Input.txt",row.names = 1, check.names = FALSE)
Design=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/HC-EC/Info.txt")
ds=DESeqDataSetFromMatrix(countData = count,colData = Design,design = ~sex+batch+group)
ds=DESeq(ds)
res=results(ds,contrast = c("group","HC","EC"),independentFiltering = FALSE)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/HC-EC/Result.txt",sep="\t", quote=FALSE,col.names = NA)


count=read.csv("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/HC-VP/Input.txt",row.names = 1, check.names = FALSE)
Design=read.csv("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/HC-VP/Info.txt")
ds=DESeqDataSetFromMatrix(countData = count,colData = Design,design = ~sex+batch+group)
ds=DESeq(ds)
res=results(ds,contrast = c("group","VP","HC"),independentFiltering = FALSE)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/HC-VP/Result.txt",sep="\t", quote=FALSE,col.names = NA)



data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/Data.txt",header=T,row.names = 1)

Art.umap = umap(data)
head(art.umap,3)
write.table(Art.umap$layout,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/Umap.txt",sep="\t",col.names = NA,quote = FALSE)
dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/Umap.txt",row.names = 1)
library(ggplot2)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=4,shape=21,aes(fill=Group))+
  scale_color_manual(values=c(VP="#b20000",HC="#0d3a1b",EC="#315b7d"))+
  scale_fill_manual(values=c(VP="#ff4c4c",HC="#96a94e",EC="#4682b4"))+
  theme(axis.title = element_text(size=9,hjust = 0.5),legend.position = c(0.91, 0.12),plot.margin = margin(0.7,0.5,0.7,0.5, "cm"))+labs(x = "UMAP1",y="UMAP2")
dev.off()



um=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/UMAP/Combat_Result.txt",header = T, row.names = 1)
library(umap)

data=read.table("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/UMAP/Combat_Result.txt",sep=",",header=T,row.names = 1)
Art.umap = umap(t(data))
head(Art.umap$layout)

write.table(Art.umap$layout,file ="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/UMAP/UMAP.txt",sep = "\t",col.names = NA,quote = FALSE )


data=read.table("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/UMAP/Combat_Filt.txt",sep=",",header=T,row.names = 1)
Art.umap = umap(t(data))
head(Art.umap$layout)

write.table(Art.umap$layout,file ="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/UMAP/UMAP_Filt.txt",sep = "\t",col.names = NA,quote = FALSE )



dat=read.csv("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/UMAP/UMAP_Filt.txt",row.names = 1)
library(ggplot2)
library(ggrepel)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/UMAP/UMAP_Filt_NewVP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=4,shape=21,aes(fill=Group))+
  scale_color_manual(values=c(VP="#b20000",HC="#0d3a1b",EC="#315b7d"))+
  scale_fill_manual(values=c(VP="#ff4c4c",HC="#96a94e",EC="#4682b4"))+geom_text_repel(data=dat,aes(x=V1,y=V2,label = row.names(dat)),size=1.9,box.padding=0.3,show.legend=FALSE,colour="black")+
  theme(axis.title = element_text(size=9,hjust = 0.5),legend.position = c(0.91, 0.9),plot.margin = margin(0.7,0.5,0.7,0.5, "cm"))+labs(x = "UMAP1",y="UMAP2")
dev.off()



data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/GeneSets.txt",header=T,row.names = 1)

Art.umap = umap(t(data))
head(art.umap,3)
write.table(Art.umap$layout,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/Umap.txt",sep="\t",col.names = NA,quote = FALSE)


dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/Umap.txt",row.names = 1)
library(ggplot2)
library(ggrepel)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=4,shape=21,aes(fill=Group))+
  scale_color_manual(values=c(VP="#b20000",HC="#0d3a1b",EC="#315b7d"))+
  scale_fill_manual(values=c(VP="#ff4c4c",HC="#96a94e",EC="#4682b4"))+
  theme(axis.title = element_text(size=9,hjust = 0.5),legend.position = c(0.93, 0.1),plot.margin = margin(0.7,0.5,0.7,0.5, "cm"))+labs(x = "UMAP1",y="UMAP2")
dev.off()


###################################################

library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/NIRF2/GeneSets.txt",header = T,row.names = 1)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/NIRF2/HeatMap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/NIRF2/Zscore.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/NIRF2/Zscore.txt",header = T, row.names = 1)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/samples.txt",row.names = 1,header = T)

col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#990000","#e50000" ,"#ff9999","white","#b5cde1","#315b7d","#23415a"))
library(ComplexHeatmap)
library(circlize)
col_fun1 = colorRamp2(c(5, 2, 0, -2, -5), c("#0000cc","#4c4cff" ,"white","#e5e500", "#b2b200"))
colours <- list("Cohort"=c("EC"="#4682b4","HC"="#96a94e","VP"="#ff4c4c"))
LFC=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/NIRF2/Exp3.txt",header = T,check.names = FALSE,row.names = 1)

H1=Heatmap(as.matrix(Zscore),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_dend_width = unit(3, "cm"),
           top_annotation =columnAnnotation(Cohort = sampleinfo$Group,col=colours,show_legend=FALSE,show_annotation_name=FALSE),
           name = "Z-Score",show_row_names = FALSE,row_names_gp=gpar(fontsize = 5),height  = unit(45, "cm"),width  = unit(20, "cm"),row_split = 5,
           column_split =c(rep("EC",14),rep("HC",12),rep("VP",16)))

col_fun_lfc = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))
H3=Heatmap(as.matrix((LFC)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="Log2FoldChange",width  = unit(3, "cm"),show_row_names = TRUE,
           row_names_gp=gpar(fontsize = 10),column_names_gp =gpar(fontsize = 15),na_col = "#e0e0e0",height  = unit(45, "cm")) 

tt=H1 + H3
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/NIRF2/HeatMap1.pdf",height = 20,width =15)
draw(tt, merge_legend = TRUE)
dev.off()


YY=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/NIRF2/BoxInput.txt",header = T)
head(YY)
library(ggpubr)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/NIRF2/NFE2L2.pdf")
ggboxplot(YY, x = "sample", y = "NFE2L2",fill="sample",order=c("HC","EC","VP"),outlier.size = 0.8,add=c("mean"),
          add.params = list(color = "white",size=0.2),ggtheme = theme_gray())+
  scale_fill_manual(values=c("EC"="#4682b4","HC"="#808080","VP"="#e5b647"))+labs(y="Log2 CPM",title = "NFE2L2")+
  theme(axis.title.x = element_blank(),plot.margin = margin(3,3,3,3, "cm"),plot.title = element_text(hjust = 0.5),
        legend.position = "none",axis.text = element_text(colour = "black"),axis.title.y  = element_text(colour = "black"))
dev.off()

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/NIRF2/KEAP1.pdf")
ggboxplot(YY, x = "sample", y = "KEAP1",fill="sample",order=c("HC","EC","VP"),outlier.size = 0.8,add=c("mean"),
          add.params = list(color = "white",size=0.2),ggtheme = theme_gray())+
  scale_fill_manual(values=c("EC"="#4682b4","HC"="#808080","VP"="#e5b647"))+labs(y="Log2 CPM",title = "KEAP1")+
  theme(axis.title.x = element_blank(),plot.margin = margin(3,3,3,3, "cm"),plot.title = element_text(hjust = 0.5),
        legend.position = "none",axis.text = element_text(colour = "black"),axis.title.y  = element_text(colour = "black"))
dev.off()
#####################################################################################


library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/BatchCorrected/New_PWY/Pwy_Combat.txt",header = T,row.names = 1)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/BatchCorrected/New_PWY/HeatMap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/BatchCorrected/New_PWY/Zscore.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/BatchCorrected/New_PWY/Zscore.txt",header = T, row.names = 1)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/BatchCorrected/New_PWY/samples.txt",row.names = 1,header = T)

col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#990000","#e50000" ,"#ff9999","white","#b5cde1","#315b7d","#23415a"))
library(ComplexHeatmap)
library(circlize)
col_fun1 = colorRamp2(c(5, 2, 0, -2, -5), c("#0000cc","#4c4cff" ,"white","#e5e500", "#b2b200"))
colours <- list("Cohort"=c("ART"="#FFA500","EC"="#4682b4","HC"="#96a94e","VP"="#ff4c4c"))
LFC=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/BatchCorrected/New_PWY/Exp5.txt",header = T,check.names = FALSE,row.names = 1)

H1=Heatmap(as.matrix(Zscore),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_dend_width = unit(3, "cm"),
           top_annotation =columnAnnotation(Cohort = sampleinfo$Group,col=colours,show_legend=FALSE,show_annotation_name=FALSE),
           name = "Z-Score",show_row_names = FALSE,row_names_gp=gpar(fontsize = 5),height  = unit(35, "cm"),width  = unit(20, "cm"),row_split = 4,
           column_split =sampleinfo$Group)

col_fun_lfc = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))
H3=Heatmap(as.matrix((LFC)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="Log2FoldChange",width  = unit(3, "cm"),show_row_names = TRUE,
           row_names_gp=gpar(fontsize = 10),column_names_gp =gpar(fontsize = 12),na_col = "#e0e0e0",height  = unit(35, "cm")) 

col_fun= c("#e0e0e0","orange")
input=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/BatchCorrected/New_PWY/Out10.txt",header = T,check.names = FALSE,row.names = 1)

H2=Heatmap(as.matrix((input)),col=col_fun,cluster_rows=FALSE,cluster_columns = FALSE,name="NULL",show_row_names = FALSE,
           column_names_gp=gpar(fontsize = 12),height  = unit(35, "cm"),width  = unit(5, "cm")) 


tt=H1 + H2+ H3
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/BatchCorrected/New_PWY/HeatMap.pdf",height = 20,width =15)
draw(tt, merge_legend = TRUE)
dev.off()

#################################################################################################

library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/GeneSets.txt",header = T,row.names = 1)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/HeatMap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/Zscore.txt",sep="\t",col.names = NA, quote = FALSE)


library(ComplexHeatmap)
library(circlize)
library(NormalyzerDE)
package.version("limma")

col_runif = colorRamp2(c(-2,-1,0,1, 2), c("#000099","#0000ff","black","#ffff4c" ,"yellow"))

colours <- list("Cohort"=c("ART"="#FFA500","EC"="#4682b4","HC"="#96a94e","VP"="#ff4c4c"))

sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/Samples.txt",row.names = 1,header = T)
colours <- list("Cohort"=c("ART"="#FFA500","EC"="#4682b4","HC"="#96a94e","VP"="#ff4c4c"))

H1=Heatmap(as.matrix(t(X$carpet)),col=col_runif,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_dend_width = unit(3, "cm"),
           top_annotation =columnAnnotation(Cohort = sampleinfo$Group,col=colours,show_legend=FALSE,show_annotation_name=FALSE),row_gap = unit(2, "mm"),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 5),height  = unit(40, "cm"),width  = unit(20, "cm"),row_split = 4,
           column_split =c(rep("ART",19),rep("EC",19),rep("HC",19),rep("VP",19)))


#col_fun_lfc = colorRamp2(c(-1,0, 1), c("#007300","#f2f2f2","#e50000"))

col_fun_lfc = colorRamp2(c(-1,-0.5, 0,0.5, 1), c("#007300","#198119","#f2f2f2","#e71919","#e50000"))

corr1=read.table("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/Correlation/HeatMap/Merged.txt",row.names = 1,header = T)
head(corr1)
H3=Heatmap(as.matrix(t(corr1)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="Correlation",width  = unit(3, "cm"),show_row_names = TRUE,
           row_names_gp=gpar(fontsize = 5),column_names_gp =gpar(fontsize = 12),na_col = "#e0e0e0",height  = unit(35, "cm")) 

tt=H1+H3
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/HeatMapNew.pdf",height = 20,width =15)
draw(tt, merge_legend = TRUE)
dev.off()

art=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/ART.txt",header = T,row.names = 1)
head(art)
library(psych)
library(reshape2)
library(reshape)
Res=corr.test(t(art),use = "pairwise",method="spearman",adjust="BH")
cor1=melt(Res$r)
pval=melt(Res$p)
cor$Adj_Pval=pval$value
write.table(cor,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/Correlation.txt",sep="\t",col.names = NA, quote = FALSE)

##########################################################

C1=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/cls1.txt")
C=melt(C1)
head(C)
write.table(C,file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/MW/cls1.txt",sep="\t",col.names = NA,quote = FALSE)
library(ggplot2)
library(ggridges)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/Cluster1.pdf")
ggplot(C, aes(x=value,y=GeneName,fill=GeneName,alpha=0.1,color=GeneName))+scale_fill_manual(values = c("ART"="#FFA500","EC"="#4682b4","HC"="#96a94e","VP"="#ff4c4c"))+
  scale_color_manual(values = c("ART"="#FFA500","EC"="#4682b4","HC"="#96a94e","VP"="#ff4c4c"))+
  stat_density_ridges(aes(fill = GeneName),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+xlim(-4,4)+
  theme(plot.margin = margin(3.5,3.5,3.5,2.5, "cm"),legend.position = "none",axis.title.y = element_blank())+labs(x="Zscore")
dev.off()


C1=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/cls2.txt")
C=melt(C1)
write.table(C,file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/MW/cls2.txt",sep="\t",col.names = NA,quote = FALSE)

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/Cluster2.pdf")
ggplot(C, aes(x=value,y=GeneName,fill=GeneName,alpha=0.1,color=GeneName))+scale_fill_manual(values = c("ART"="#FFA500","EC"="#4682b4","HC"="#96a94e","VP"="#ff4c4c"))+
  scale_color_manual(values = c("ART"="#FFA500","EC"="#4682b4","HC"="#96a94e","VP"="#ff4c4c"))+
  stat_density_ridges(aes(fill = GeneName),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+xlim(-4,4)+
  theme(plot.margin = margin(3.5,3.5,3.5,2.5, "cm"),legend.position = "none",axis.title.y = element_blank())+labs(x="Zscore")
dev.off()


C1=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/cls3.txt")
C=melt(C1)
write.table(C,file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/MW/cls3.txt",sep="\t",col.names = NA,quote = FALSE)

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/Cluster3.pdf")
ggplot(C, aes(x=value,y=GeneName,fill=GeneName,alpha=0.1,color=GeneName))+scale_fill_manual(values = c("ART"="#FFA500","EC"="#4682b4","HC"="#96a94e","VP"="#ff4c4c"))+
  scale_color_manual(values = c("ART"="#FFA500","EC"="#4682b4","HC"="#96a94e","VP"="#ff4c4c"))+
  stat_density_ridges(aes(fill = GeneName),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+xlim(-4,4)+
  theme(plot.margin = margin(3.5,3.5,3.5,2.5, "cm"),legend.position = "none",axis.title.y = element_blank())+labs(x="Zscore")
dev.off()

C1=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/cls4.txt")
C=melt(C1)
write.table(C,file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/MW/cls4.txt",sep="\t",col.names = NA,quote = FALSE)

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/Cluster4.pdf")
ggplot(C, aes(x=value,y=GeneName,fill=GeneName,alpha=0.1,color=GeneName))+scale_fill_manual(values = c("ART"="#FFA500","EC"="#4682b4","HC"="#96a94e","VP"="#ff4c4c"))+
  scale_color_manual(values = c("ART"="#FFA500","EC"="#4682b4","HC"="#96a94e","VP"="#ff4c4c"))+
  stat_density_ridges(aes(fill = GeneName),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+xlim(-4,4)+
  theme(plot.margin = margin(3.5,3.5,3.5,2.5, "cm"),legend.position = "none",axis.title.y = element_blank())+labs(x="Zscore")
dev.off()


blah=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/blah.txt")
head(blah)

wilcox.test(value ~ Group, data=blah)
chisq.test(blah)
library("dplyr")
mu <- C %>% 
  group_by(GeneName) %>%
  summarise(grp.mean = mean(value))
?wilcox.test
#############################################

X3=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/For_Box.txt",header = T)
head(X3)
library(ggpubr)
my_comparisons <- list( c("ART", "EC"), c("ART", "VP"), c("ART", "HC"),c("EC", "HC"),c("EC", "VP"),c("HC", "VP") )
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/SOD1.pdf")
ggboxplot(X3, x = "Group", y = "SOD1",fill="Group",order=c("ART","EC","HC","VP"),outlier.shape =NA,ggtheme = theme_gray())+
  scale_fill_manual(values=c(VP="#ff4c4c",HC="#96a94e",EC="#4682b4",ART="#FFA500"))+
  geom_jitter(color="black",shape=16,size=0.7, position=position_jitter(0.08))+labs(y="Log2CPM",title = "SOD1")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "wilcox.test", paired = FALSE)+
  theme(axis.title.x = element_blank(),plot.margin = margin(3,3,3,3, "cm"),plot.title = element_text(hjust = 0.5,size = 11),
        legend.position = "none",axis.text = element_text(colour = "black"))
dev.off()


pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/BNIP3.pdf")
ggboxplot(X3, x = "Group", y = "BNIP3",fill="Group",order=c("ART","EC","HC","VP"),outlier.shape =NA,ggtheme = theme_gray())+
  scale_fill_manual(values=c(VP="#ff4c4c",HC="#96a94e",EC="#4682b4",ART="#FFA500"))+
  geom_jitter(color="black",shape=16,size=0.7, position=position_jitter(0.08))+labs(y="Log2CPM",title = "BNIP3")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "wilcox.test", paired = FALSE)+
  theme(axis.title.x = element_blank(),plot.margin = margin(3,3,3,3, "cm"),plot.title = element_text(hjust = 0.5,size = 11),
        legend.position = "none",axis.text = element_text(colour = "black"))
dev.off()

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/APOE.pdf")
ggboxplot(X3, x = "Group", y = "APOE",fill="Group",order=c("ART","EC","HC","VP"),outlier.shape =NA,ggtheme = theme_gray())+
  scale_fill_manual(values=c(VP="#ff4c4c",HC="#96a94e",EC="#4682b4",ART="#FFA500"))+
  geom_jitter(color="black",shape=16,size=0.7, position=position_jitter(0.08))+labs(y="Log2CPM",title = "APOE")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "wilcox.test", paired = FALSE)+
  theme(axis.title.x = element_blank(),plot.margin = margin(3,3,3,3, "cm"),plot.title = element_text(hjust = 0.5,size = 11),
        legend.position = "none",axis.text = element_text(colour = "black"))
dev.off()

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/MAPK7.pdf")
ggboxplot(X3, x = "Group", y = "MAPK7",fill="Group",order=c("ART","EC","HC","VP"),outlier.shape =NA,ggtheme = theme_gray())+
  scale_fill_manual(values=c(VP="#ff4c4c",HC="#96a94e",EC="#4682b4",ART="#FFA500"))+
  geom_jitter(color="black",shape=16,size=0.7, position=position_jitter(0.08))+labs(y="Log2CPM",title = "MAPK7")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "wilcox.test", paired = FALSE)+
  theme(axis.title.x = element_blank(),plot.margin = margin(3,3,3,3, "cm"),plot.title = element_text(hjust = 0.5,size = 11),
        legend.position = "none",axis.text = element_text(colour = "black"))
dev.off()

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/HMOX1.pdf")
ggboxplot(X3, x = "Group", y = "HMOX1",fill="Group",order=c("ART","EC","HC","VP"),outlier.shape =NA,ggtheme = theme_gray())+
  scale_fill_manual(values=c(VP="#ff4c4c",HC="#96a94e",EC="#4682b4",ART="#FFA500"))+
  geom_jitter(color="black",shape=16,size=0.7, position=position_jitter(0.08))+labs(y="Log2CPM",title = "HMOX1")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "wilcox.test", paired = FALSE)+
  theme(axis.title.x = element_blank(),plot.margin = margin(3,3,3,3, "cm"),plot.title = element_text(hjust = 0.5,size = 11),
        legend.position = "none",axis.text = element_text(colour = "black"))
dev.off()

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/NQO1.pdf")
ggboxplot(X3, x = "Group", y = "NQO1",fill="Group",order=c("ART","EC","HC","VP"),outlier.shape =NA,ggtheme = theme_gray())+
  scale_fill_manual(values=c(VP="#ff4c4c",HC="#96a94e",EC="#4682b4",ART="#FFA500"))+
  geom_jitter(color="black",shape=16,size=0.7, position=position_jitter(0.08))+labs(y="Log2CPM",title = "NQO1")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "wilcox.test", paired = FALSE)+
  theme(axis.title.x = element_blank(),plot.margin = margin(3,3,3,3, "cm"),plot.title = element_text(hjust = 0.5,size = 11),
        legend.position = "none",axis.text = element_text(colour = "black"))
dev.off()

Da=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/MW/MW_Result_cls1.tab",row.names = 1)
head(Da)
Da$BH =p.adjust(Da$P_value,method = "BH")
write.table(Da,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/MW/MW_Result_cls1_adj.tab",sep="\t",col.names = NA,quote = FALSE)

###########################################################################################################

library(psych)
clini=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/Correlation/Clinical_Variables.txt",check.names=FALSE)
head(clini)
gene=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/Correlation/ART_BatchCorr.txt",check.names=FALSE)
cliniX=as.matrix(clini)
geneX=as.matrix(gene)
Res=corr.test(cliniX,geneX,use = "pairwise",method="spearman",adjust="none")

library(reshape2)
m=melt(Res$p,id.vars = "ID")
head(m)
write.table(m,file ="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/Correlation/Correlation_PvalueCorrRaw.txt",col.names = NA,quote = FALSE,sep="\t")

m1=melt(Res$r,id.vars = "ID")
head(m1)
write.table(m1,file ="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/Correlation/CorrelationCorrRaw.txt",col.names = NA,quote = FALSE,sep="\t")
##################################################################

clini=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/Correlation/Clinical_Variables.txt",check.names=FALSE)
head(clini)
gene=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/Correlation/GeneSets.txt",check.names=FALSE)
cliniX=as.matrix(clini)
geneX=as.matrix(gene)
Res=corr.test(cliniX,geneX,use = "pairwise",method="spearman",adjust="none")
warnings()
library(reshape2)
m=melt(Res$p,id.vars = "ID")
head(m)
write.table(m,file ="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/Correlation/Correlation_PvalueRaw.txt",col.names = NA,quote = FALSE,sep="\t")

m1=melt(Res$r,id.vars = "ID")
head(m1)
write.table(m1,file ="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/Correlation/CorrelationRaw.txt",col.names = NA,quote = FALSE,sep="\t")



YY=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/Correlation/RegressionResultBatchCorrected.txt")
head(YY)
YY$AdjPVal =p.adjust(YY$pvalue,method = "BH")
write.table(YY,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/Correlation/RegressionResultBatchCorrectedADJ.txt",sep="\t",col.names = NA,quote = FALSE)

#########################

DAT=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/Correlation/HeatMap/Melted.txt")
head(DAT)
Re <- cast(DAT)
head(Re)

write.table(Re,file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/Correlation/HeatMap/Merged.txt",sep="\t",col.names = NA,quote = FALSE)



col_fun_lfc = colorRamp2(c(-1,0, 1), c("#008000","grey","#e50000"))
corr1=read.table("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/Correlation/HeatMap/Merged.txt",row.names = 1,header = T)
head(corr1)
H3=Heatmap(as.matrix(t(corr1)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="Log2FoldChange",width  = unit(3, "cm"),show_row_names = TRUE,
           row_names_gp=gpar(fontsize = 10),column_names_gp =gpar(fontsize = 12),na_col = "#e0e0e0",height  = unit(35, "cm")) 


pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/Correlation/HeatMap/HeatMapNew.pdf",height = 20,width =15)
Heatmap(as.matrix(t(corr1)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="Correlation",width  = unit(10, "cm"),show_row_names = TRUE,
        row_names_gp=gpar(fontsize = 10),column_names_gp =gpar(fontsize = 12),na_col = "#e0e0e0",height  = unit(35, "cm"))
dev.off()

#####################################################################
library(umap)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/UMAP/ART_T.txt",header=T,row.names = 1)
head(data)
Art.umap = umap((data))
head(Art.umap$layout,3)
write.table(Art.umap$layout,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/Umap.txt",sep="\t",col.names = NA,quote = FALSE)
dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/Umap.txt",row.names = 1)
library(ggplot2)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/maike/UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=4,shape=21,aes(fill=Group))+
  scale_color_manual(values=c(VP="#b20000",HC="#0d3a1b",EC="#315b7d"))+
  scale_fill_manual(values=c(VP="#ff4c4c",HC="#96a94e",EC="#4682b4"))+
  theme(axis.title = element_text(size=9,hjust = 0.5),legend.position = c(0.91, 0.12),plot.margin = margin(0.7,0.5,0.7,0.5, "cm"))+labs(x = "UMAP1",y="UMAP2")
dev.off()


library(M3C)
?tsne

################################################ ART Only Heatmap

library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/ART_Only/ART_HC.txt",header = T,row.names = 1)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/ART_Only/HeatMap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

col_runif = colorRamp2(c(-2,-1,0,1, 2), c("#000099","#0000ff","black","#ffff4c" ,"yellow"))

sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/ART_Only/samples.txt",row.names = 1,header = T)
colours <- list("Cohort"=c("sART"="#ffdb99","L-ART"="#cc8400","HC"="#96a94e"))

H1=Heatmap(as.matrix(t(X$carpet)),col=col_runif,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = TRUE,row_dend_width = unit(3, "cm"),
           column_order = sort(colnames(t(X$carpet))),
           top_annotation =columnAnnotation(Cohort = sampleinfo$Group,col=colours,show_legend=FALSE,show_annotation_name=FALSE),row_gap = unit(2, "mm"),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 5),height  = unit(40, "cm"),width  = unit(20, "cm"),row_split = 4,
           column_split =c(rep("sART",9),rep("L-ART",10),rep("HC",19)))


col_fun_lfc = colorRamp2(c(-1,-0.5, 0,0.5, 1), c("#007300","#198119","#f2f2f2","#e71919","#e50000"))

corr1=read.table("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/ART_Only/Mean.txt",row.names = 1,header = T)
head(corr1)
meanClr = colorRamp2(c(-5,-2,0,2,5), c("#007300","#b2d8b2","white" ,"#ffcccc","#ff0000"))
H3=Heatmap(as.matrix((corr1)),col=meanClr,cluster_rows=FALSE,cluster_columns = FALSE,name="MeanExp",width  = unit(3, "cm"),show_row_names = TRUE,
           row_names_gp=gpar(fontsize = 5),column_names_gp =gpar(fontsize = 12),height  = unit(40, "cm")) 
row.names(corr1)
tt=H1+H3
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/ART_Only/HeatMap.pdf",height = 20,width =15)
draw(tt, merge_legend = TRUE)
dev.off()


library(umap)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/ART_Only/ART_HC.txt",header=T,row.names = 1)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/ART_Only/Umap.txt",sep="\t",col.names = NA,quote = FALSE)
dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/ART_Only/Umap.txt",row.names = 1)
library(ggplot2)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/ART_Only/UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=4,shape=21,aes(fill=Group))+
  scale_color_manual(values=c("sART"="#ffd27f","L-ART"="#b27300",HC="#0d3a1b"))+
  scale_fill_manual(values=c("sART"="#ffdb99","L-ART"="#cc8400",HC="#96a94e"))+
  theme(axis.title = element_text(size=9,hjust = 0.5),legend.position = c(0.1, 0.12),plot.margin = margin(0.7,0.5,0.7,0.5, "cm"))+labs(x = "UMAP1",y="UMAP2")
dev.off()



library(umap)

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/UMAP/Log.txt",header=T,row.names = 1)
head(data)
LData=log10(data)
Art.umap = umap((data))
write.table(LData,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/UMAP/Log.txt",sep="\t",col.names = NA,quote = FALSE)

head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/UMAP/Umap.txt",sep="\t",col.names = NA,quote = FALSE)
dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/ART_Only/Umap.txt",row.names = 1)




dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/UMAP/Umap.txt",row.names = 1)
library(ggplot2)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/ROS/UMAP/UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=4,shape=21,aes(fill=Group))+
  scale_color_manual(values=c(VP="#b20000",HC="#0d3a1b",EC="#315b7d",ART="#b27300"))+
  scale_fill_manual(values=c(VP="#ff4c4c",HC="#96a94e",EC="#4682b4",ART="#FFA500"))+
  theme(axis.title = element_text(size=9,hjust = 0.5),legend.position = c(0.07, 0.12),plot.margin = margin(0.7,0.5,0.7,0.5, "cm"))+labs(x = "UMAP1",y="UMAP2")
dev.off()

################################# HIF Targetted genes

library(gplots)
Dat=read.csv("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/HIF/GeneSets.txt",header = T,row.names = 1)
head(Dat)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/HIF/HeatMap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/HIF/Zscore.txt",sep="\t",col.names = NA, quote = FALSE)

library(ComplexHeatmap)
library(circlize)

col_runif = colorRamp2(c(-2,-1,0,1, 2), c("#000099","#0000ff","black","#ffff4c" ,"yellow"))

sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/HIF/Samples.txt",row.names = 1,header = T)
colours <- list("Cohort"=c("ART"="#FFA500","EC"="#4682b4","B_HC"="#96a94e","VP"="#ff4c4c"))
head(sampleinfo)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/HIF/HeatMap.pdf",height = 20,width =15)
Heatmap(as.matrix(t(X$carpet)),col=col_runif,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_dend_width = unit(4, "cm"),
        top_annotation =columnAnnotation(Cohort = sampleinfo$Group,col=colours,show_legend=FALSE,show_annotation_name=FALSE),row_gap = unit(2, "mm"),
        name = "Z-Score",show_row_names = FALSE,row_names_gp=gpar(fontsize = 5),height  = unit(40, "cm"),width  = unit(20, "cm"),row_split = 4,
        column_split =c(rep("ART",19),rep("EC",19),rep("B_HC",19),rep("VP",19)))

dev.off()

library(reshape2)
C1=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/HIF/cluster1.txt")
C=melt(C1)
head(C,20)
library(ggplot2)
library(ggridges)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/HIF/cluster1.pdf")
ggplot(C, aes(x=value,y=GeneName,fill=GeneName,alpha=0.1,color=GeneName))+scale_fill_manual(values = c("ART"="#FFA500","EC"="#4682b4","B_HC"="#96a94e","VP"="#ff4c4c"))+
  scale_color_manual(values = c("ART"="#FFA500","EC"="#4682b4","B_HC"="#96a94e","VP"="#ff4c4c"))+
  stat_density_ridges(aes(fill = GeneName),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+xlim(-4,4)+
  theme(plot.margin = margin(3.5,3.5,3.5,2.5, "cm"),legend.position = "none",axis.title.y = element_blank())+labs(x="Zscore")
dev.off()


library(reshape2)
C1=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/HIF/cluster2.txt")
C=melt(C1)
head(C,20)
library(ggplot2)
library(ggridges)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/HIF/cluster2.pdf")
ggplot(C, aes(x=value,y=GeneName,fill=GeneName,alpha=0.1,color=GeneName))+scale_fill_manual(values = c("ART"="#FFA500","EC"="#4682b4","B_HC"="#96a94e","VP"="#ff4c4c"))+
  scale_color_manual(values = c("ART"="#FFA500","EC"="#4682b4","B_HC"="#96a94e","VP"="#ff4c4c"))+
  stat_density_ridges(aes(fill = GeneName),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+xlim(-4,4)+
  theme(plot.margin = margin(3.5,3.5,3.5,2.5, "cm"),legend.position = "none",axis.title.y = element_blank())+labs(x="Zscore")
dev.off()

C1=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/HIF/cluster3.txt")
C=melt(C1)
head(C,20)
library(ggplot2)
library(ggridges)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/HIF/cluster3.pdf")
ggplot(C, aes(x=value,y=GeneName,fill=GeneName,alpha=0.1,color=GeneName))+scale_fill_manual(values = c("ART"="#FFA500","EC"="#4682b4","B_HC"="#96a94e","VP"="#ff4c4c"))+
  scale_color_manual(values = c("ART"="#FFA500","EC"="#4682b4","B_HC"="#96a94e","VP"="#ff4c4c"))+
  stat_density_ridges(aes(fill = GeneName),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+xlim(-4,4)+
  theme(plot.margin = margin(3.5,3.5,3.5,2.5, "cm"),legend.position = "none",axis.title.y = element_blank())+labs(x="Zscore")
dev.off()

C1=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/HIF/cluster4.txt")
C=melt(C1)
head(C,20)
library(ggplot2)
library(ggridges)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/HIF/cluster4.pdf")
ggplot(C, aes(x=value,y=GeneName,fill=GeneName,alpha=0.1,color=GeneName))+scale_fill_manual(values = c("ART"="#FFA500","EC"="#4682b4","B_HC"="#96a94e","VP"="#ff4c4c"))+
  scale_color_manual(values = c("ART"="#FFA500","EC"="#4682b4","B_HC"="#96a94e","VP"="#ff4c4c"))+
  stat_density_ridges(aes(fill = GeneName),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+xlim(-4,4)+
  theme(plot.margin = margin(3.5,3.5,3.5,2.5, "cm"),legend.position = "none",axis.title.y = element_blank())+labs(x="Zscore")
dev.off()

################################ NRF2 targetted genes

Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/NRF2/GeneSets.txt",header = T,row.names = 1)
head(Dat)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/NRF2/HeatMap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/NRF2/Zscore.txt",sep="\t",col.names = NA, quote = FALSE)

col_runif = colorRamp2(c(-2,-1,0,1, 2), c("#000099","#0000ff","black","#ffff4c" ,"yellow"))

sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/HIF/Samples.txt",row.names = 1,header = T)
colours <- list("Cohort"=c("ART"="#FFA500","EC"="#4682b4","B_HC"="#96a94e","VP"="#ff4c4c"))
head(sampleinfo)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/NRF2/HeatMap.pdf",height = 20,width =15)
Heatmap(as.matrix(t(X$carpet)),col=col_runif,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_dend_width = unit(4, "cm"),
        top_annotation =columnAnnotation(Cohort = sampleinfo$Group,col=colours,show_legend=FALSE,show_annotation_name=FALSE),row_gap = unit(2, "mm"),
        name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 8),height  = unit(43, "cm"),width  = unit(20, "cm"),row_split = 3,
        column_split =c(rep("ART",19),rep("EC",19),rep("B_HC",19),rep("VP",19)))

dev.off()


C1=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/NRF2/cluster1.txt")
C=melt(C1)
head(C,20)
library(ggplot2)
library(ggridges)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/NRF2/cluster1.pdf")
ggplot(C, aes(x=value,y=GeneName,fill=GeneName,alpha=0.1,color=GeneName))+scale_fill_manual(values = c("ART"="#FFA500","EC"="#4682b4","B_HC"="#96a94e","VP"="#ff4c4c"))+
  scale_color_manual(values = c("ART"="#FFA500","EC"="#4682b4","B_HC"="#96a94e","VP"="#ff4c4c"))+
  stat_density_ridges(aes(fill = GeneName),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+xlim(-4,4)+
  theme(plot.margin = margin(3.5,3.5,3.5,2.5, "cm"),legend.position = "none",axis.title.y = element_blank())+labs(x="Zscore")
dev.off()


C1=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/NRF2/cluster2.txt")
C=melt(C1)
head(C,20)
library(ggplot2)
library(ggridges)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/NRF2/cluster2.pdf")
ggplot(C, aes(x=value,y=GeneName,fill=GeneName,alpha=0.1,color=GeneName))+scale_fill_manual(values = c("ART"="#FFA500","EC"="#4682b4","B_HC"="#96a94e","VP"="#ff4c4c"))+
  scale_color_manual(values = c("ART"="#FFA500","EC"="#4682b4","B_HC"="#96a94e","VP"="#ff4c4c"))+
  stat_density_ridges(aes(fill = GeneName),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+xlim(-4,4)+
  theme(plot.margin = margin(3.5,3.5,3.5,2.5, "cm"),legend.position = "none",axis.title.y = element_blank())+labs(x="Zscore")
dev.off()


C1=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/NRF2/cluster3.txt")
C=melt(C1)
head(C,20)
library(ggplot2)
library(ggridges)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Hall_Mark_GeneSets/NRF2/cluster3.pdf")
ggplot(C, aes(x=value,y=GeneName,fill=GeneName,alpha=0.1,color=GeneName))+scale_fill_manual(values = c("ART"="#FFA500","EC"="#4682b4","B_HC"="#96a94e","VP"="#ff4c4c"))+
  scale_color_manual(values = c("ART"="#FFA500","EC"="#4682b4","B_HC"="#96a94e","VP"="#ff4c4c"))+
  stat_density_ridges(aes(fill = GeneName),quantile_lines = TRUE,quantiles = 0.5,linetype="dashed")+xlim(-4,4)+
  theme(plot.margin = margin(3.5,3.5,3.5,2.5, "cm"),legend.position = "none",axis.title.y = element_blank())+labs(x="Zscore")
dev.off()


library(countToFPKM)
cnt=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/Network/Coding_Count_Hg38.txt",row.names = 1)

anno=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/Network/GeneLength.txt")

sam=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/Network/FS.txt")
head(sam)
meanFragmentLength <- sam$FS
featureLength <- anno$length
fpkm_matrix <- fpkm (as.matrix(cnt), featureLength, meanFragmentLength)
head(fpkm_matrix)
write.table(fpkm_matrix,file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/Network/Coding_FPKM.txt",sep="\t",quote = FALSE,col.names = NA)


counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))
  
  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}

TPM=counts_to_tpm(cnt,featureLength, meanFragmentLength)
head(TPM)
write.table(TPM,file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/Network/Coding_TPM.txt",sep="\t",quote = FALSE,col.names = NA)

#####################

library(countToFPKM)
cnt=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/NonCoding_Count_Hg38.txt",row.names = 1)

anno=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/Network/GeneLength_NonCod.txt")

sam=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/Network/FS.txt")
head(sam)
meanFragmentLength <- sam$FS
featureLength <- anno$length
fpkm_matrix <- fpkm (as.matrix(cnt), featureLength, meanFragmentLength)
head(fpkm_matrix)
write.table(fpkm_matrix,file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/Network/Coding_FPKM.txt",sep="\t",quote = FALSE,col.names = NA)


counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))
  
  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}

TPM=counts_to_tpm(cnt,featureLength, meanFragmentLength)
head(TPM)
write.table(TPM,file="/home/anoop/Desktop/COVID_Omics/Transcriptomics/Network/Coding_TPM.txt",sep="\t",quote = FALSE,col.names = NA)


################ Nrew volcano EC -vs- ART

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/Volcano_Result.txt",row.names = 1,sep = ",")
head(data)
library("extrafont")
loadfonts()
font_import()
library(ggrepel)
fonts()
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/Volcano_Result.pdf",family = "Arial")
ggplot(data, aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(data=subset(data, padj>=0.05),aes(x=log2FoldChange, y=-log10(padj),size=-log10(padj)),color="#bfbfbf")+
  geom_label_repel(aes(label = Top),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,segment.alpha=0.75,
                   box.padding=0,nudge_x = 0.55,nudge_y = 1.2)+
  geom_point(data=subset(data, padj<0.05 & log2FoldChange < 0),aes(x=log2FoldChange, y=-log10(padj),color=abs(log2FoldChange),size=-log10(padj)))+scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, padj<0.05 & log2FoldChange > 0),shape=21,aes(x=log2FoldChange, y=-log10(padj),fill=abs(log2FoldChange),size=-log10(padj)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")
dev.off()

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/Volcano_ResultART.txt",row.names = 1,sep = ",")
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/Volcano_ResultART.pdf",family = "Arial")
ggplot(data, aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(data=subset(data, padj>=0.05),aes(x=log2FoldChange, y=-log10(padj),size=-log10(padj)),color="#bfbfbf")+
  geom_point(data=subset(data, padj<0.05 & log2FoldChange < 0 & ART==1),aes(x=log2FoldChange, y=-log10(padj),color=abs(log2FoldChange),size=-log10(padj)))+scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, padj<0.05 & log2FoldChange > 0 & ART==1),shape=21,aes(x=log2FoldChange, y=-log10(padj),fill=abs(log2FoldChange),size=-log10(padj)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")
dev.off()



library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/HeatMap/NEW/Signi.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/HeatMap/NEW/HeatMap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/HeatMap/NEW/ZScore.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/HeatMap/NEW/ZScore1.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/HeatMap/NEW/Info1.txt",row.names = 1)

col_fun1 = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))

library(ComplexHeatmap)
library(circlize)

colours <- list("Group"=c("HC"="#48690E","Conv"="#C6E2FF","Mild"="#ffd700","Severe"="#ffa500"))

ha = HeatmapAnnotation(df = sampleinfo,show_annotation_name = FALSE,annotation_name_side = "right",annotation_legend_param = list(Group = list(direction = "horizontal",grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 18))),
                       col = list(Group=c("EC"="#1919ff","ART"="#ffa500")))

col_age = colorRamp2(c(25,30,40,50,60,70), c("#6ed7c5","#3ecab2","#0ebd9f","#0b977f","#08715f","#054b3f"))
col_art=colorRamp2(c(6,8,10,12,15,20,25), c("#ffffff","#dcecaa","#cee588","#c0dd66","#bada55","#94ae44","#6f8233"))
col_hiv=colorRamp2(c(0,0.1,0.5,1,5,10,25), c("#b7e7ef","#aee4ed","#a6e2eb","#95cbd3","#84b4bc","#749ea4","#63878d"))

ann=read.delim("/home/anoop/Desktop/COVID_Omics/Transcriptomics/DiffExp/Deseq/Figures/Annotation.txt",row.names = 1,check.names = FALSE)
head(ann)

ha = HeatmapAnnotation(na_col = "#999999",df = sampleinfo,show_annotation_name = FALSE,annotation_name_side = "right",
                       annotation_legend_param = list(Age = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                 grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                 labels_gp = gpar(fontsize = 18)),
                                                      Duration_ART = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                labels_gp = gpar(fontsize = 18)),
                                                      Duration_HIVPos = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                 grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                 labels_gp = gpar(fontsize = 18)),
                                                      Gender = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                             grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                             labels_gp = gpar(fontsize = 18)),
                                                      Group = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                             grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                             labels_gp = gpar(fontsize = 18))),
                       col = list(Age=col_age,Duration_ART=col_art,Duration_HIVPos=col_hiv,Group=c("EC"="#1919ff","ART"="#ffa500"),Gender=c("Male"="#5e348a","Female"="#ffc0cb")))

colnames(sampleinfo)
col_fun1 = colorRamp2(c(4,1.5,1, 0,-1,-1.5,-4), c("#d67834","#da8548" ,"#de935c","white","#6590b4","#5282ab","#3f75a2"))

col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7F7F00","#B2B200" ,"#E5E500","white","#BF7FBF","#993299","#590059"))
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows=TRUE,cluster_columns = TRUE,show_column_names = TRUE,row_title_gp = gpar(fontsize=20),
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),
           top_annotation  =ha,heatmap_legend_param =list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 15), labels_gp = gpar(fontsize = 15)),
           name = "Z-Score",show_row_names = FALSE,row_names_gp=gpar(fontsize = 14),height  = unit(30, "cm"),width  = unit(22, "cm")
           )

LFC=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/HeatMap/NEW/LFC.txt",row.names = 1)

col_fun_lfc = colorRamp2(c(-5, -3,-2,-1, 0,1,2,3,5), c("#0000cc","#4c4cff","#5d5dff","#6f6fff" ,"white","#c1c132","#b9b919","#b2b200","#7f7f00"))
col_fun1 = colorRamp2(c(-2, -1,0, 1,2), c("#004C00","#008000","white","#E50000","#7F0000"))
H3=Heatmap(as.matrix((LFC)),col=col_fun1,cluster_rows=FALSE,cluster_columns = FALSE,name="LFC",width  = unit(0.75, "cm"),show_column_names = FALSE,
           show_row_names = FALSE,column_names_side = "top",heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                        title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           row_names_gp =gpar(fontsize = 7),height  = unit(30, "cm"),column_names_gp =gpar(fontsize = 20),na_col = "#e6e6e6") 

tt=H1+H3
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/HeatMap/NEW/HeatMap.pdf",height = 20,width =20)
draw(tt,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()




###############
y
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/GSEA/HallMark.txt",row.names = 1,sep = "\t")
head(data)
library("extrafont")
loadfonts()
font_import()
library(ggrepel)
fonts()
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/GSEA/HallMark1.pdf",family = "Arial")
ggplot(data, aes(x=Order, y=NES)) + theme_bw()+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,segment.alpha=1,
                   box.padding=0.95,nudge_x = 1.8,nudge_y = -0.01)+
  geom_point(data=subset(data, NES>0),aes(x=Order, y=NES, size=Log),color="#90001c",fill="#ba0024")+geom_hline(yintercept=0,size=0.3)+
  geom_point(data=subset(data, NES<0),aes(x=Order, y=NES, size=Log),color="#004000",fill="#008000")+
  scale_size(range = c(3, 10))+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        axis.title.y=element_text(size=15),legend.position = "bottom",axis.ticks.x = element_blank(),panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),axis.text.x =element_text(size=10,color="black",vjust =61),plot.margin = margin(3,0.5,3,0.5, "cm"))+
  guides(size=guide_legend(override.aes=list(color="grey",fill="grey")))+
  labs(y="NES(EC-vs-ART)",size="-Log10(FDR)")
dev.off()


data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/GSEA/GO.txt",row.names = 1,sep = "\t")
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/GSEA/GO.pdf",family = "Arial")
ggplot(data, aes(x=Order, y=NES)) + theme_bw()+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,segment.alpha=1,
                   box.padding=0.95,nudge_x = 1.8,nudge_y = -0.01)+
  geom_point(data=subset(data, NES>0),aes(x=Order, y=NES, size=Log),shape=21,color="#90001c",fill="#ba0024")+geom_hline(yintercept=0,size=0.3)+
  geom_point(data=subset(data, NES<0),aes(x=Order, y=NES, size=Log),shape=21,color="#004000",fill="#008000")+
  scale_size(range = c(3, 10))+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        axis.title.y=element_text(size=15),legend.position = "bottom",axis.ticks.x = element_blank(),panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),axis.text.x =element_text(size=10,color="black",vjust =67),plot.margin = margin(3,0.5,3,0.5, "cm"))+
  guides(size=guide_legend(override.aes=list(color="grey",fill="grey")))+
  labs(y="NES(EC-vs-ART)",size="-Log10(FDR)")
dev.off()




data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/Result.txt",row.names = 1)
head(data)
library(ggrepel)
svg("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/MAPlot.svg")
ggplot(data, aes(x=log2(baseMean), y=log2FoldChange)) + theme_bw()+
  geom_point(data=subset(data, padj>=0.05),aes(log2(baseMean),y=log2FoldChange),color="#bfbfbf",size=1)+
  geom_point(data=subset(data, padj<0.05 & log2FoldChange < 0),aes(x=log2(baseMean),y=log2FoldChange),color="#2e84d5",size=1)+
  geom_point(data=subset(data, padj<0.05 & log2FoldChange > 0),aes(x=log2(baseMean),y=log2FoldChange),color="#b20000",size=1)+
  guides(fill = guide_legend(title = "Regulation",override.aes = aes(label = "")))+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=25),legend.position = "none",
        axis.title.x=element_text(size=25),axis.text=element_text(size=20,color="black"),plot.margin = margin(2,1.5,2,1, "cm"))+
  labs(x="Log2 mean expression",y="Log2 fold change")
dev.off()



dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Umap/ART_Umap.txt",row.names = 1,check.names = FALSE)

library(ggrepel)
library(ggplot2)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ART-Specific/Umap/Figure_2_a_name.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=7,shape=21,aes(fill=Group))+theme_bw()+
  scale_color_manual(values=c(VP="#800000",HC="#666666",EC="#0000ff",ART="#e59400"))+
  scale_fill_manual(values=c(VP="#8c1919",HC="#808080",EC="#1919ff",ART="#ffa500"))+
  theme(axis.title = element_text(size=20,hjust = 0.5),legend.title = element_blank(),legend.text = element_text(size=20),
        plot.margin = margin(0.7,0.5,0.7,0.5, "cm"),
        legend.position = c(0.91, 0.12))+labs(x = "UMAP1",y="UMAP2")+
geom_text_repel(data=dat,aes(x=V1,y=V2,label = ID),size=1.9,box.padding=0.3,show.legend=FALSE,colour="black")
dev.off()



library(ggplot2)
library(dplyr)
library(vegan)
Genus=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EPIC/Res1.txt",check.names=FALSE,row.names=1)
meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EPIC/M1.txt",check.names=FALSE,row.names = 1)

head(t(Genus))
permanova <- adonis(t(Genus) ~ Group,data = meta, permutations=99, method = "bray")
print(as.data.frame(permanova$aov.tab)["Group", "Pr(>F)"])
dist <- vegdist(t(Genus))
anova(betadisper(dist, meta$Group))
coef <- coefficients(permanova)["group1",]
top.coef <- coef[rev(order(abs(coef)))[1:20]]
library(phyloseq)
data(dietswap)
pseq <- dietswap




##########

bar=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EPIC/Median.txt",header = TRUE)
bar$Cell <- factor(bar$Cell, levels = bar$Cell)
bar$row <- 1
bar
P1=ggplot(bar, aes(x = row,y = ART, fill = Cell)) +
  geom_bar(stat="identity",width = .5) + ylab("Mean cell fraction")+theme_bw()+xlab("ART")+
  scale_fill_manual(values = c("CD4_Tcells"="#e3b7d2","CD8_Tcells"="#f6cebf","Monocytes"="#bfe7f6","NKcells"="#d0ffcb","Bcells"="#ffdada","Neutrophils"="#feffde"))+
  geom_text(aes(label = ART), position = position_stack(vjust = 0.5))+theme(axis.title.x = element_text(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),
                                                                            plot.margin = margin(1.5,3.1,1.5,-0.9, "cm"),panel.grid.major.x = element_blank(),
                                                                             panel.grid.minor.x = element_blank(),axis.title.y = element_blank(),
                                                                             axis.ticks.x=element_blank(),legend.position = "none")

P2=ggplot(bar, aes(x = row,y = EC, fill = Cell)) +
  geom_bar(stat="identity",width = .5) + ylab("Percentage")+theme_bw()+xlab("EC")+
  scale_fill_manual(values = c("CD4_Tcells"="#e3b7d2","CD8_Tcells"="#f6cebf","Monocytes"="#bfe7f6","NKcells"="#d0ffcb","Bcells"="#ffdada","Neutrophils"="#feffde"))+
  geom_text(aes(label = EC), position = position_stack(vjust = 0.5))+theme(axis.title.x = element_text(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),
                                                                           plot.margin = margin(1.5,5.15,1.5,-3, "cm"),panel.grid.major.x = element_blank(),
                                                                            panel.grid.minor.x = element_blank(),axis.title.y = element_blank(),
                                                                            axis.ticks.x=element_blank(),legend.position = "none")

P3=ggplot(bar, aes(x = row,y = HC, fill = Cell)) +
  geom_bar(stat="identity",width = .5) + ylab("Percentage")+theme_bw()+xlab("HC")+
  scale_fill_manual(values = c("CD4_Tcells"="#e3b7d2","CD8_Tcells"="#f6cebf","Monocytes"="#bfe7f6","NKcells"="#d0ffcb","Bcells"="#ffdada","Neutrophils"="#feffde"))+
  geom_text(aes(label = HC), position = position_stack(vjust = 0.5))+theme(axis.title.x = element_text(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),
                                                                           plot.margin = margin(1.5,4,1.5,-5, "cm"),panel.grid.major.x = element_blank(),
                                                                           panel.grid.minor.x = element_blank(),axis.title.y = element_blank(),
                                                                           axis.ticks.x=element_blank(),legend.position = "right")

P4=ggplot(bar, aes(x = row,y = VP, fill = Cell)) +
  geom_bar(stat="identity",width = .5) + ylab("Median Cell Fraction")+theme_bw()+xlab("VP")+
  scale_fill_manual(values = c("CD4_Tcells"="#e3b7d2","CD8_Tcells"="#f6cebf","Monocytes"="#bfe7f6","NKcells"="#d0ffcb","Bcells"="#ffdada","Neutrophils"="#feffde"))+
  geom_text(aes(label = VP), position = position_stack(vjust = 0.5))+theme(axis.title.x = element_text(),axis.text.x = element_blank(),
                                                                           plot.margin = margin(1.5,1,1.5,0.1, "cm"),panel.grid.major.x = element_blank(),
                                                                           panel.grid.minor.x = element_blank(),
                                                                           axis.ticks.x=element_blank(),legend.position = "none")
library(ggpubr)



################# Metabolomics
library(umap)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/MERGED.txt",header=TRUE,row.names = 1,check.names = FALSE)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/UMAP.txt",sep="\t",col.names = NA,quote = FALSE)


library(PCAtools)
count=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/MERGED.txt",header=TRUE,row.names = 1,check.names = FALSE)
head(count)
meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/MetaData.txt",row.names = 1)
p <- pca(count, metadata = meta, removeVar = 0.1)
head(p$variance)

write.table(p$rotated,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/PCA.txt",sep="\t",col.names = NA,quote = FALSE)


dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/PCA.txt",row.names = 1)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/PCA.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Batch)) + geom_point(size=5,aes(fill=Batch))+
  theme(axis.title = element_text(size=12,hjust = 0.5),legend.title = element_blank(),legend.text = element_text(size=12),
        plot.margin = margin(0.7,0.5,0.7,0.5, "cm"),
        legend.position = c(0.1, 0.1))+labs(x = "PC1(27.33% variation)",y="PC2(21.91% variation)")#+
#geom_text_repel(data=dat,aes(x=PC1,y=PC2,label =Label),size=1.9,box.padding=0.3,show.legend=FALSE,colour="black")
dev.off()

dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/UMAP.txt",row.names = 1)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Batch)) + geom_point(size=5,shape=21,aes(fill=Batch))+
  theme(axis.title = element_text(size=12,hjust = 0.5),legend.title = element_blank(),legend.text = element_text(size=12),
        plot.margin = margin(0.7,0.5,0.7,0.5, "cm"),
        legend.position = c(0.89, 0.1))+labs(x = "UMAP1",y="UMAP2")+
  geom_text_repel(data=dat,aes(x=V1,y=V2,label = rownames(dat)),size=1.9,box.padding=0.3,show.legend=FALSE,colour="black")
dev.off()



count=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/HC.txt",header=TRUE,row.names = 1,check.names = FALSE)
head(count)
meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Meta.txt",row.names = 1)
p <- pca(count, metadata = meta, removeVar = 0.1)
head(p$variance)

write.table(p$rotated,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/PCA_HC.txt",sep="\t",col.names = NA,quote = FALSE)


dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/PCA_HC.txt",row.names = 1)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/PCA.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Group)) + geom_point(size=5,shape=21,aes(fill=Group))+
  theme(axis.title = element_text(size=12,hjust = 0.5),legend.title = element_blank(),legend.text = element_text(size=12),
        plot.margin = margin(0.7,0.5,0.7,0.5, "cm"),
        legend.position = c(0.89, 0.1))+labs(x = "PC1(22.3% variation)",y="PC2(21.03% variation)")+
  geom_text_repel(data=dat,aes(x=PC1,y=PC2,label =rownames(dat)),size=1.9,box.padding=0.3,show.legend=FALSE,colour="black")
dev.off()



count=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Filtered/Data.txt",header=TRUE,row.names = 1,check.names = FALSE)
head(count)
meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Filtered/Metadata.txt",row.names = 1)
p <- pca(count, metadata = meta, removeVar = 0.1)
head(p$variance)

write.table(p$rotated,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Filtered/PCA.txt",sep="\t",col.names = NA,quote = FALSE)


dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Filtered/PCA.txt",row.names = 1)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Filtered/PCA.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Batch)) + geom_point(size=5,shape=21,aes(fill=Batch))+
  theme(axis.title = element_text(size=12,hjust = 0.5),legend.title = element_blank(),legend.text = element_text(size=12),
        plot.margin = margin(0.7,0.5,0.7,0.5, "cm"),
        legend.position = c(0.1, 0.1))+labs(x = "PC1(33.45% variation)",y="PC2(27.74% variation)")
dev.off()

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Filtered/Data.txt",header=TRUE,row.names = 1,check.names = FALSE)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Filtered/UMAP.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Filtered/UMAP.txt",row.names = 1)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Filtered/UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=5,shape=21,aes(fill=Group))+
  theme(axis.title = element_text(size=12,hjust = 0.5),legend.title = element_blank(),legend.text = element_text(size=12),
        plot.margin = margin(0.7,0.5,0.7,0.5, "cm"),
        legend.position = c(0.1, 0.1))+labs(x = "UMAP1",y="UMAP2")+
  geom_text_repel(data=dat,aes(x=V1,y=V2,label =Label),size=1.9,box.padding=0.3,show.legend=FALSE,colour="black")
dev.off()


library(sva)
prot=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Filtered/Data.txt",row.names = 1,check.names = FALSE)
mat <- as.matrix(prot)
des=read.table("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Filtered/Metadata.txt",sep="\t",header=TRUE)
designCombat = model.matrix(~ des$Cohort)
rnaseqCombat = ComBat(mat, batch = des$Batch, mod = designCombat, par.prior = TRUE, prior.plots = TRUE)
write.table(rnaseqCombat,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Filtered/Combat_Result.txt", sep="\t",quote=FALSE,col.names = NA)


data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Filtered/Combat_Result.txt",header=TRUE,row.names = 1,check.names = FALSE)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Filtered/UMAP_combat.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Filtered/UMAP_combat.txt",row.names = 1)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Filtered/UMAP_combat.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=5,shape=21,aes(fill=Group))+
  theme(axis.title = element_text(size=12,hjust = 0.5),legend.title = element_blank(),legend.text = element_text(size=12),
        plot.margin = margin(0.7,0.5,0.7,0.5, "cm"),
        legend.position = c(0.1, 0.1))+labs(x = "UMAP1",y="UMAP2")+
  geom_text_repel(data=dat,aes(x=V1,y=V2,label =Label),size=1.9,box.padding=0.3,show.legend=FALSE,colour="black")
dev.off()



setwd("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/")
rm(list = ls())
library(DESeq2)
count=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Long_Short_ART/Input.txt",row.names = 1, check.names = FALSE)
Design=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Long_Short_ART/Info.txt")
ds=DESeqDataSetFromMatrix(countData = count,colData = Design,design = ~sex+batch+group)
ds=DESeq(ds)
res=results(ds,contrast = c("group","long","short"),independentFiltering = FALSE)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Long_Short_ART/Result.txt",sep="\t", quote=FALSE,col.names = NA)
?results

count=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Long_Short_ART/IP_EC_long.txt",row.names = 1, check.names = FALSE)
Design=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Long_Short_ART/Info_EC_Long.txt")
ds=DESeqDataSetFromMatrix(countData = count,colData = Design,design = ~ group + sex+group:nestedBatch+group:sex)
ds=DESeq(ds)
res=results(ds,contrast = c("group","Long","EC"),independentFiltering = FALSE)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Long_Short_ART/Result_EC_Long.txt",sep="\t", quote=FALSE,col.names = NA)


count=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Long_Short_ART/IP_EC_short.txt",row.names = 1, check.names = FALSE)
Design=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Long_Short_ART/Info_EC_short.txt")
ds=DESeqDataSetFromMatrix(countData = count,colData = Design,design = ~ group + sex+group:nestedBatch+group:sex)
ds=DESeq(ds)
res=results(ds,contrast = c("group","Short","EC"),independentFiltering = FALSE)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Long_Short_ART/Result_EC_short.txt",sep="\t", quote=FALSE,col.names = NA)



########################## redox genes heatmap

library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ReporterMets/RedoxIP.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ReporterMets/RedoxHeatmap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ReporterMets/RedoxZ.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ReporterMets/RedoxZ.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ReporterMets/Info1.txt",row.names = 1)

col_fun1 = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))

library(ComplexHeatmap)
library(circlize)

ha = HeatmapAnnotation(na_col = "#999999",df = sampleinfo,show_annotation_name = FALSE,annotation_name_side = "right",
                       annotation_legend_param = list(Group = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                   grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                   labels_gp = gpar(fontsize = 18))),
                       col = list(Group=c("a_EC"="#1919ff","b_ART"="#ffa500")))


col_fun1 = colorRamp2(c(4,1.5,1, 0,-1,-1.5,-4), c("#d67834","#da8548" ,"#de935c","white","#6590b4","#5282ab","#3f75a2"))

col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7F7F00","#B2B200" ,"#E5E500","white","#8caabe","#407294","#335b76"))

col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#987316","#daa520" ,"#e5c062","white","#8caabe","#407294","#335b76"))
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=20),
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),
           top_annotation  =ha,heatmap_legend_param =list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                          labels_gp = gpar(fontsize = 19)),
           name = "Z-Score",show_row_names = FALSE,row_names_gp=gpar(fontsize = 14),height  = unit(30, "cm"),width  = unit(22, "cm"),
           column_split =c(rep("a_EC",19),rep("b_ART",19))
)

LFC=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ReporterMets/LFC.txt",row.names = 1)

col_fun_lfc = colorRamp2(c(-5, -3,-2,-1, 0,1,2,3,5), c("#0000cc","#4c4cff","#5d5dff","#6f6fff" ,"white","#c1c132","#b9b919","#b2b200","#7f7f00"))
col_fun2 = colorRamp2(c(-2, -1,0, 1,2), c("#004C00","#008000","white","#E50000","#7F0000"))


ha3 = rowAnnotation(foo = anno_mark(at = c(36,25,160,139,48,20,168,92,3,96,80,37,155),
                                    labels_gp = gpar(fontsize=20),lines_gp = gpar(col="black"),link_height = unit(35, "mm"),link_width=unit(8, "mm"),
                                    labels = c("ALDOA","FASN","TXNDC12","GLS","ALDH6A1","IDH3B","PARK7","PRPS1","SOD1","PRDX4","TXN","TXNDC11","GSTM4")))

H3=Heatmap(as.matrix((LFC)),col=col_fun2,cluster_rows=FALSE,cluster_columns = FALSE,name="LFC",width  = unit(0.75, "cm"),show_column_names = FALSE,
           show_row_names = FALSE,column_names_side = "top",heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                        title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           row_names_gp =gpar(fontsize = 7),height  = unit(30, "cm"),column_names_gp =gpar(fontsize = 20),na_col = "#e6e6e6",right_annotation = ha3) 

tt=H1+H3
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ReporterMets/RedoxHeatmap.pdf",height = 20,width =20)
draw(tt,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()

?results
###########

library(umap)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ReporterMets/RedoxIP_2.txt",header=TRUE,row.names = 1,check.names = FALSE)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ReporterMets/RedoxUMAP.txt",sep="\t",col.names = NA,quote = FALSE)



library(ggrepel)
library(ggplot2)
dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ReporterMets/RedoxUMAP.txt",row.names = 1)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ReporterMets/RedoxUMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=7,shape=21,aes(fill=Group))+theme_bw()+
  scale_color_manual(values=c(EC="#0000ff",ART="#e59400"))+
  scale_fill_manual(values=c(EC="#1919ff",ART="#ffa500"))+
  theme(axis.title = element_text(size=20,hjust = 0.5),legend.title = element_blank(),legend.text = element_text(size=20),
        plot.margin = margin(0.7,0.5,0.7,0.5, "cm"),
        legend.position = c(0.1, 0.12))+labs(x = "UMAP1",y="UMAP2")
dev.off()


count=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/HC.txt",header=TRUE,row.names = 1,check.names = FALSE)
head(count)
meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Meta.txt",row.names = 1)
p <- pca(count, metadata = meta, removeVar = 0.1)
head(p$variance)

######################

library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ReporterMets/RedoxIP.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ReporterMets/RedoxHeatmapAll.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ReporterMets/RedoxZ-All.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ReporterMets/RedoxZ-All.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ReporterMets/Info2.txt",row.names = 1)


library(ComplexHeatmap)
library(circlize)

ha = HeatmapAnnotation(na_col = "#999999",df = sampleinfo,show_annotation_name = FALSE,annotation_name_side = "right",
                       annotation_legend_param = list(Group = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                   grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                   labels_gp = gpar(fontsize = 18))),
                       col = list(Group=c(b_VP="#8c1919",a_HC="#808080",c_EC="#1919ff",d_ART="#ffa500")))


col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#987316","#daa520" ,"#e5c062","white","#8caabe","#407294","#335b76"))
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=20),
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),
           top_annotation  =ha,heatmap_legend_param =list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                          labels_gp = gpar(fontsize = 19)),
           name = "Z-Score",show_row_names = FALSE,row_names_gp=gpar(fontsize = 14),height  = unit(30, "cm"),width  = unit(22, "cm"),
           column_split =c(rep("a_HC",19),rep("b_VP",19),rep("c_EC",19),rep("d_ART",19))
)

H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows=TRUE,cluster_columns = TRUE,show_column_names = FALSE,row_title_gp = gpar(fontsize=20),
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),
           top_annotation  =ha,heatmap_legend_param =list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                          labels_gp = gpar(fontsize = 19)),
           name = "Z-Score",show_row_names = FALSE,row_names_gp=gpar(fontsize = 14),height  = unit(30, "cm"),width  = unit(22, "cm"))

LFC=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ReporterMets/LFC.txt",row.names = 1)

col_fun_lfc = colorRamp2(c(-5, -3,-2,-1, 0,1,2,3,5), c("#0000cc","#4c4cff","#5d5dff","#6f6fff" ,"white","#c1c132","#b9b919","#b2b200","#7f7f00"))
col_fun2 = colorRamp2(c(-2, -1,0, 1,2), c("#004C00","#008000","white","#E50000","#7F0000"))


ha3 = rowAnnotation(foo = anno_mark(at = c(22,34,152,59,18,88,111,130,140,127,149,151,78),
                                    labels_gp = gpar(fontsize=20),lines_gp = gpar(col="black"),link_height = unit(35, "mm"),link_width=unit(8, "mm"),
                                    labels = c("ALDOA","FASN","TXNDC12","GLS","ALDH6A1","IDH3B","PARK7","PRPS1","SOD1","PRDX4","TXN","TXNDC11","GSTM4")))

H3=Heatmap(as.matrix((LFC)),col=col_fun2,cluster_rows=FALSE,cluster_columns = FALSE,name="LFC",width  = unit(0.75, "cm"),show_column_names = FALSE,
           show_row_names = FALSE,column_names_side = "top",heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                        title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           row_names_gp =gpar(fontsize = 7),height  = unit(30, "cm"),column_names_gp =gpar(fontsize = 20),na_col = "#e6e6e6",right_annotation = ha3) 

tt=H1+H3
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/ReporterMets/RedoxHeatmapAll2.pdf",height = 20,width =20)
draw(tt,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()


library(umap)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Test/Data.txt",header=TRUE,row.names = 1,check.names = FALSE)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Test/Umap.txt",sep="\t",col.names = NA,quote = FALSE)
library(ggplot2)
library(ggrepel)
dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Test/Umap.txt",row.names = 1)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Test/Umap.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=5,shape=21,aes(fill=Group))+
  theme(axis.title = element_text(size=12,hjust = 0.5),legend.title = element_blank(),legend.text = element_text(size=12),
        plot.margin = margin(0.7,0.5,0.7,0.5, "cm"),
        legend.position = c(0.1, 0.1))+labs(x = "UMAP1",y="UMAP2")+
  geom_text_repel(data=dat,aes(x=V1,y=V2,label =Label),size=1.9,box.padding=0.3,show.legend=FALSE,colour="black")
dev.off()


DAT=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Limma/Combat_Result.txt",row.names = 1,check.names = FALSE)
L=log2(DAT)
write.table(L,file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Limma/Log2Data.txt",sep="\t",col.names = NA,quote = FALSE)

library("limma")
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Limma/Info.txt")
group <- paste(sampleinfo$Group)
group <- factor(group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Limma/Log2Data.txt",sep="\t",row.names=1,check.names = FALSE)
mat <- as.matrix(data)
fit <- lmFit(mat,design)
cont.matrix <- makeContrasts(ECvsART=ART - EC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="ECvsART",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Limma/EC_ART.txt",sep="\t",quote=FALSE,col.names = NA)

cont.matrix <- makeContrasts(HCvsART=ART - HC,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,coef="HCvsART",sort.by="p",n="Inf")
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/Metabolomics/Limma/HC_ART.txt",sep="\t",quote=FALSE,col.names = NA)



library(umap)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/FBA/Results/Combat.txt",header=TRUE,row.names = 1,check.names = FALSE)
head(data)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/FBA/Results/Umap.txt",sep="\t",col.names = NA,quote = FALSE)
library(ggplot2)
library(ggrepel)
library(ggrepel)
library(ggplot2)
dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/FBA/Results/Umap.txt",row.names = 1)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/FBA/Results/Umap.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=7,shape=21,aes(fill=Group))+theme_bw()+
  scale_color_manual(values=c(EC="#0000ff",ART="#e59400"))+
  scale_fill_manual(values=c(EC="#1919ff",ART="#ffa500"))+
  theme(axis.title = element_text(size=20,hjust = 0.5),legend.title = element_blank(),legend.text = element_text(size=20),
        plot.margin = margin(0.7,0.5,0.7,0.5, "cm"),
        legend.position = c(0.1, 0.12))+labs(x = "UMAP1",y="UMAP2")
dev.off()


library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/FBA/Temp/Name_Input.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/FBA/Temp/RXN_genes.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/FBA/Temp/Name_InputZ.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/FBA/Temp/Name_InputZ.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/FBA/Temp/Info1.txt",row.names = 1)


library(ComplexHeatmap)
library(circlize)



ha = HeatmapAnnotation(df = sampleinfo,show_annotation_name = FALSE,annotation_name_side = "right",annotation_legend_param = list(Group = list(direction = "horizontal",grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 18))),
                       col = list(Group=c("EC"="#1919ff","ART"="#ffa500")))

col_age = colorRamp2(c(25,30,40,50,60,70), c("#6ed7c5","#3ecab2","#0ebd9f","#0b977f","#08715f","#054b3f"))
col_art=colorRamp2(c(6,8,10,12,15,20,25), c("#ffffff","#dcecaa","#cee588","#c0dd66","#bada55","#94ae44","#6f8233"))



ha = HeatmapAnnotation(na_col = "#999999",df = sampleinfo,show_annotation_name = FALSE,annotation_name_side = "right",
                       annotation_legend_param = list(Age = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                 grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                 labels_gp = gpar(fontsize = 18)),
                                                      Duration_ART = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                          grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                          labels_gp = gpar(fontsize = 18)),
                                                      Gender = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                    grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                    labels_gp = gpar(fontsize = 18)),
                                                      Group = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                   grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                   labels_gp = gpar(fontsize = 18))),
                       col = list(Age=col_age,Duration_ART=col_art,Group=c("EC"="#1919ff","ART"="#ffa500"),Gender=c("Male"="#5e348a","Female"="#ffc0cb")))

colnames(sampleinfo)
col_fun1 = colorRamp2(c(4,1.5,1, 0,-1,-1.5,-4), c("#d67834","#da8548" ,"#de935c","white","#6590b4","#5282ab","#3f75a2"))

col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#987316","#daa520" ,"#e5c062","white","#8caabe","#407294","#335b76"))
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows=TRUE,cluster_columns = TRUE,show_column_names = FALSE,row_title_gp = gpar(fontsize=20),
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),
           top_annotation  =ha,heatmap_legend_param =list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 15), 
                                                          labels_gp = gpar(fontsize = 15)),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 14),height  = unit(30, "cm"),width  = unit(22, "cm")
)

LFC=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/EC-ART/HeatMap/NEW/LFC.txt",row.names = 1)

col_fun_lfc = colorRamp2(c(-5, -3,-2,-1, 0,1,2,3,5), c("#0000cc","#4c4cff","#5d5dff","#6f6fff" ,"white","#c1c132","#b9b919","#b2b200","#7f7f00"))
col_fun1 = colorRamp2(c(-2, -1,0, 1,2), c("#004C00","#008000","white","#E50000","#7F0000"))
H3=Heatmap(as.matrix((LFC)),col=col_fun1,cluster_rows=FALSE,cluster_columns = FALSE,name="LFC",width  = unit(0.75, "cm"),show_column_names = FALSE,
           show_row_names = FALSE,column_names_side = "top",heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                        title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           row_names_gp =gpar(fontsize = 7),height  = unit(30, "cm"),column_names_gp =gpar(fontsize = 20),na_col = "#e6e6e6") 

tt=H1+H3
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/FBA/Temp/RXN_genes.pdf",height = 20,width =20)
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()



######################

library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/FBA/Temp/data/Reporter_IP.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/FBA/Temp/data/ReporterMetsGenes.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/FBA/Temp/data/Reporter_Z.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/FBA/Temp/data/Reporter_Z.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/FBA/Temp/data/Info1.txt",row.names = 1)


library(ComplexHeatmap)
library(circlize)


col_age = colorRamp2(c(25,30,40,50,60,70), c("#6ed7c5","#3ecab2","#0ebd9f","#0b977f","#08715f","#054b3f"))
col_art=colorRamp2(c(6,8,10,12,15,20,25), c("#ffffff","#dcecaa","#cee588","#c0dd66","#bada55","#94ae44","#6f8233"))



ha = HeatmapAnnotation(na_col = "#999999",df = sampleinfo,show_annotation_name = FALSE,annotation_name_side = "right",
                       annotation_legend_param = list(Age = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                 grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                 labels_gp = gpar(fontsize = 18)),
                                                      Duration_ART = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                          grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                          labels_gp = gpar(fontsize = 18)),
                                                      Gender = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                    grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                    labels_gp = gpar(fontsize = 18)),
                                                      Group = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                   grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                   labels_gp = gpar(fontsize = 18))),
                       col = list(Age=col_age,Duration_ART=col_art,Group=c("EC"="#1919ff","ART"="#ffa500"),Gender=c("Male"="#5e348a","Female"="#ffc0cb")))

colnames(sampleinfo)
col_fun1 = colorRamp2(c(4,1.5,1, 0,-1,-1.5,-4), c("#d67834","#da8548" ,"#de935c","white","#6590b4","#5282ab","#3f75a2"))

col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#987316","#daa520" ,"#e5c062","white","#797087","#60576d","#4a4355"))
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows=TRUE,cluster_columns = TRUE,show_column_names = FALSE,row_title_gp = gpar(fontsize=20),
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),
           top_annotation  =ha,heatmap_legend_param =list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 15), 
                                                          labels_gp = gpar(fontsize = 15)),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 14),height  = unit(14, "cm"),width  = unit(18, "cm")
)


pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/FBA/Temp/data/ReporterMetsGenes.pdf",height = 20,width =20)
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()


df=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/FBA/Test/artCount.txt")
head(df)
ggplot(data=df, aes(x=key, y=frac)) + geom_line()


#######################

library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/FBA/Temp/NewHeatMap/Filtered_CombatNames.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/FBA/Temp/NewHeatMap/RXNsAssocGenes.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/FBA/Temp/NewHeatMap/Filtered_CombatNamesZ.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/FBA/Temp/NewHeatMap/Filtered_CombatNamesZ.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/FBA/Temp/NewHeatMap/Info1.txt",row.names = 1)


library(ComplexHeatmap)
library(circlize)

col_age = colorRamp2(c(25,30,40,50,60,70), c("#6ed7c5","#3ecab2","#0ebd9f","#0b977f","#08715f","#054b3f"))
col_art=colorRamp2(c(6,8,10,12,15,20,25), c("#ffffff","#dcecaa","#cee588","#c0dd66","#bada55","#94ae44","#6f8233"))



ha = HeatmapAnnotation(na_col = "#999999",df = sampleinfo,show_annotation_name = FALSE,annotation_name_side = "right",
                       annotation_legend_param = list(Age = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                 grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                 labels_gp = gpar(fontsize = 18)),
                                                      Duration_ART = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                          grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                          labels_gp = gpar(fontsize = 18)),
                                                      Gender = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                    grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                    labels_gp = gpar(fontsize = 18)),
                                                      Group = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                   grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                   labels_gp = gpar(fontsize = 18))),
                       col = list(Age=col_age,Duration_ART=col_art,Group=c("EC"="#1919ff","ART"="#ffa500","HC"="#808080"),Gender=c("Male"="#5e348a","Female"="#ffc0cb")))

geneinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/FBA/Temp/NewHeatMap/GeneInfo.txt",row.names = 1)

ha2 = rowAnnotation(na_col = "#999999",df = geneinfo,show_annotation_name = FALSE,
                        annotation_legend_param = list(subsystem = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                  grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                  labels_gp = gpar(fontsize = 18))),
                                                       col = list(subsystem=c("TCA cycle"="#b5ddd1","Glycolysis"="#d3c0f9","Oxidative phosphorylation"="#f99a9c")))


col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#987316","#daa520" ,"#e5c062","white","#8caabe","#407294","#335b76"))
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows=TRUE,cluster_columns = TRUE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),right_annotation =ha2,
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),row_split = geneinfo$subsystem,
           top_annotation  =ha,heatmap_legend_param =list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 15), 
                                                          labels_gp = gpar(fontsize = 15)),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 14),height  = unit(30, "cm"),width  = unit(22, "cm")
)

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/FBA/Temp/NewHeatMap/RXNsAssocGenes.pdf",height = 20,width =20)
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()



################### for Ujjwal

count=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/For_Ujjwal/Names_Combat_Result.txt",row.names = 1,check.names = FALSE)
countX=as.matrix(count)
calc_coef_var <- function(x) sd(x) / mean(x)
coef_var <- apply(countX, 1, calc_coef_var)
HVG_5 <- countX[rank(coef_var) / length(coef_var) > 0.95, ]
write.table(HVG_5,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/For_Ujjwal/Combat_Top5per.txt",sep="\t",col.names = NA,quote = FALSE)



library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/For_Ujjwal/Combat_Top10per.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/For_Ujjwal/Top_10percHeatmap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/For_Ujjwal/Combat_Top10per_Z.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/For_Ujjwal/Combat_Top10per_Z.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/For_Ujjwal/MetaData.txt",row.names = 1)
library(ComplexHeatmap)
library(circlize)

col_age = colorRamp2(c(25,30,40,50,60,70), c("#6ed7c5","#3ecab2","#0ebd9f","#0b977f","#08715f","#054b3f"))

ha = HeatmapAnnotation(na_col = "#999999",df = sampleinfo,show_annotation_name = FALSE,annotation_name_side = "right",
                       annotation_legend_param = list(Age = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                 grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                 labels_gp = gpar(fontsize = 18)),
                                                      Sex = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                    grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                    labels_gp = gpar(fontsize = 18)),
                                                      Group = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                   grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                   labels_gp = gpar(fontsize = 18))),
                       col = list(Age=col_age,Group=c(VP="#8c1919",HC="#808080",EC="#1919ff"),Sex=c("Male"="#b266b2","Female"="#ffc0cb")))


col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#987316","#daa520" ,"#e5c062","white","#8caabe","#407294","#335b76"))
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows=TRUE,cluster_columns = TRUE,show_column_names = FALSE,row_title_gp = gpar(fontsize=20),
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),
           top_annotation  =ha,heatmap_legend_param =list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 15), 
                                                          labels_gp = gpar(fontsize = 15)),
           name = "Z-Score",show_row_names = FALSE,row_names_gp=gpar(fontsize = 14),height  = unit(30, "cm"),width  = unit(20, "cm")
)


pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/For_Ujjwal/Top_10percHeatmap.pdf",height = 20,width =20)
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()

###########

library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/For_Ujjwal/Combat_Top5per.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/For_Ujjwal/Top_5percHeatmap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/For_Ujjwal/Combat_Top5per_Z.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/For_Ujjwal/Combat_Top5per_Z.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/For_Ujjwal/MetaData.txt",row.names = 1)
library(ComplexHeatmap)
library(circlize)

col_age = colorRamp2(c(25,30,40,50,60,70), c("#6ed7c5","#3ecab2","#0ebd9f","#0b977f","#08715f","#054b3f"))

ha = HeatmapAnnotation(na_col = "#999999",df = sampleinfo,show_annotation_name = FALSE,annotation_name_side = "right",
                       annotation_legend_param = list(Age = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                 grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                 labels_gp = gpar(fontsize = 18)),
                                                      Sex = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                 grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                 labels_gp = gpar(fontsize = 18)),
                                                      Group = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                   grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                   labels_gp = gpar(fontsize = 18))),
                       col = list(Age=col_age,Group=c(VP="#8c1919",HC="#808080",EC="#1919ff"),Sex=c("Male"="#747e96","Female"="#ffc0cb")))


col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#987316","#daa520" ,"#e5c062","white","#8caabe","#407294","#335b76"))

col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7F7F00","#B2B200" ,"#E5E500","white","#BF7FBF","#993299","#590059"))
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows=TRUE,cluster_columns = TRUE,show_column_names = FALSE,row_title_gp = gpar(fontsize=20),
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),
           top_annotation  =ha,heatmap_legend_param =list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 15), 
                                                          labels_gp = gpar(fontsize = 15)),
           name = "Z-Score",show_row_names = FALSE,row_names_gp=gpar(fontsize = 14),height  = unit(25, "cm"),width  = unit(20, "cm")
)


pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/For_Ujjwal/Top_5percHeatmap.pdf",height = 20,width =20)
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()

############

library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/For_Ujjwal/Proteom/HeatMap.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/For_Ujjwal/Proteom/Heatamap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/For_Ujjwal/Proteom/heatmap_Z.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/For_Ujjwal/Proteom/heatmap_Z.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/For_Ujjwal/Proteom/Meta.txt",row.names = 1)
library(ComplexHeatmap)
library(circlize)

ha = HeatmapAnnotation(na_col = "#999999",df = sampleinfo,show_annotation_name = FALSE,annotation_name_side = "right",
                       annotation_legend_param = list(Group = list(direction = "horizontal",grid_width = unit(1, "cm"),
                                                                   grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 19), 
                                                                   labels_gp = gpar(fontsize = 18))),
                       col = list(Group=c(VP="#8c1919",HC="#808080",EC="#1919ff")))


col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#987316","#daa520" ,"#e5c062","white","#8caabe","#407294","#335b76"))

col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7F7F00","#B2B200" ,"#E5E500","white","#BF7FBF","#993299","#590059"))
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows=TRUE,cluster_columns = TRUE,show_column_names = FALSE,row_title_gp = gpar(fontsize=20),
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),
           top_annotation  =ha,heatmap_legend_param =list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 15), 
                                                          labels_gp = gpar(fontsize = 15)),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 14),height  = unit(20, "cm"),width  = unit(25, "cm")
)


pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/For_Ujjwal/Proteom/Heatamap.pdf",height = 20,width =20)
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()

################################# NEW EPIC 

library(EPIC)
tdata <- read.delim(file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NEW_EPIC/Coding_TPM2.txt", as.is=T, check.names=F)
bulk <- as.matrix(tdata[,-(1:2)])
rownames(bulk) <- tdata[,1]

rdata <- read.delim(file="/home/anoop/Desktop/COVID_Omics/EPIC/NEW/Blood_Cells_casted.txt", as.is=T, check.names=F)
REF <- as.matrix(rdata[,-(1:2)])
rownames(REF) <- rdata[,1]
head(REF)

sig=read.delim("/home/anoop/Desktop/COVID_Omics/EPIC/NEW/NewSigGenes.txt",header = FALSE)
SIG <- as.vector(sig[,1])

refList <- list(refProfiles=REF,sigGenes=SIG)

Res2 <- EPIC(bulk = bulk, reference = refList,scaleExprs=FALSE)

?EPIC

write.table(Res2$mRNAProportions,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NEW_EPIC/mRNAProp.txt",col.names = NA,quote = FALSE,sep="\t")
write.table(Res2$cellFractions,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NEW_EPIC/CellFrac.txt",col.names = NA,quote = FALSE,sep="\t")
write.table(Res2$fit.gof,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NEW_EPIC/OtherRes.txt",col.names = NA,quote = FALSE,sep="\t")


##########

Data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NEW_EPIC/KW/CellFrac.txt",header = TRUE,check.names = FALSE)
head(Data)
library(reshape2)
M=melt(Data)
head(M)
library(ggpubr)
library(ggplot2)
my_comparisons = list( c("HC", "ART"),c("HC", "EC"),c("EC", "ART"),c("VP", "ART"))
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NEW_EPIC/KW/CellFracKW.pdf",width = 17.5,height = 12)
ggplot(M,aes(x=factor(Group,levels = unique(Group)),y=value,fill=Group,color=Group))+geom_boxplot(outlier.shape = NA)+facet_wrap(~ variable, ncol = 6,scales ="free" )+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  scale_color_manual(values=c(VP="#B20000",HC="#0D3A1B",EC="#315B7D",ART="#B27300"))+
  scale_fill_manual(values=c(VP="#FF4C4C",HC="#96A94E",EC="#4682B4",ART="#FFA500"))+
  theme_bw()+stat_compare_means(label = "p.format",comparisons = my_comparisons)+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position ="bottom",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        plot.margin = margin(2,0.1,2,0.1, "cm"),strip.text.x = element_text(size = 12),
        axis.text.y = element_text(color="black",size=12),plot.title = element_text(hjust = 0.5,size = 15))+
  guides(color = guide_legend(nrow = 2),fill = guide_legend(nrow = 2))
dev.off()


##########################  Correlation

library(psych)
library(reshape2)
tpm=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NEW_EPIC/Network/ART.txt",check.names=FALSE,header = TRUE)

Corr=corr.test(as.matrix(tpm),use = "pairwise",method="spearman",adjust="BH")

write.table(Corr$r,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NEW_EPIC/Network/ART_R.txt",sep="\t",col.names = NA,quote = FALSE)
corr=melt(Corr$r)
pval=melt(Corr$p)

corr$BH=pval
write.table(corr,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NEW_EPIC/Network/ART_R.txt",sep="\t",col.names = NA,quote = FALSE)

############
tpm=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NEW_EPIC/Network/EC.txt",check.names=FALSE,header = TRUE)

Corr=corr.test(as.matrix(tpm),use = "pairwise",method="spearman",adjust="BH")
corr=melt(Corr$r)
pval=melt(Corr$p)

corr$BH=pval
write.table(corr,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NEW_EPIC/Network/EC_R.txt",sep="\t",col.names = NA,quote = FALSE)

###################

tpm=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NEW_EPIC/Network/HC.txt",check.names=FALSE,header = TRUE)

Corr=corr.test(as.matrix(tpm),use = "pairwise",method="spearman",adjust="BH")
corr=melt(Corr$r)
pval=melt(Corr$p)

corr$BH=pval
write.table(corr,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NEW_EPIC/Network/HC_R.txt",sep="\t",col.names = NA,quote = FALSE)

###########

tpm=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NEW_EPIC/Network/VP.txt",check.names=FALSE,header = TRUE)

Corr=corr.test(as.matrix(tpm),use = "pairwise",method="spearman",adjust="BH")
corr=melt(Corr$r)
pval=melt(Corr$p)

corr$BH=pval
write.table(corr,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NEW_EPIC/Network/VP_R.txt",sep="\t",col.names = NA,quote = FALSE)

packageVersion("DESeq2")

############################################# Flux Balance Analysis #############
library(ComplexHeatmap)
library(circlize)
FBA=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/FBA/Merged/Analyse/InputPlot.txt",row.names = 1)
sub=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/FBA/Merged/Analyse/SubSytem2.txt",row.names = 1)

sub_anno = rowAnnotation(df = sub, width =unit(0.2, "mm"),annotation_name_gp = gpar(fontsize = 0),
                         annotation_legend_param  = list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),
                                                         title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
                         col = list(SubSystem = c("Bile acid biosynthesis"="#303964","Bile acid recycling"="#596083",
                                                  "Carnitine shuttle (mitochondrial)"="#979cb1","Fatty acid oxidation"="#d5d7e0",
                                                  "Folate metabolism"="#9c1d5d","Glycolysis / Gluconeogenesis"="#df2a85",
                                                  "Nucleotide metabolism"="#eb7fb5","Oxidative phosphorylation"="#f5bfda",
                                                  "Pentose phosphate pathway"="#004c00","Purine metabolism"="#008000",
                                                  "Pyrimidine metabolism"="#4ca64c","Other"="#93A6A3","Starch and sucrose metabolism"="#99cc99","Transport reactions"="#E1CDB6")))


ha = rowAnnotation(foo = anno_mark(at = c(40,41,50,52,53,54,61,66,67,70,71,73,77,79,81,87,88,92,93,99,108,110,116,120,130,148,153,179,180,185,233,235,247,251,252),
                                   labels_gp = gpar(fontsize=15),lines_gp = gpar(col="black"),
                                   link_height = unit(13, "mm"),link_width=unit(10, "mm"),
                                   labels = c("LCAT1e","C181CRNt","HMR_4655","HMR_8144","HMR_4442","HMR_3923","HMR_7876",
                                              "HMR_8470","HMR_7870","HMR_7885","HMR_7891","HMR_4680","HMR_7882","HMR_8462","HMR_8461",
                                              "HMR_8453","HMR_7874","HMR_7888","HMR_7872","HMR_3827","HMR_7663","HMR_8667","HMR_2626","HMR_0226","HMR_1329",
                                              "HMR_4398","HMR_4710","HMR_8917","HMR_4914","HMR_5579","HMR_5585","HMR_4964","HMR_5008","HMR_5541","MLTHFte3")))

chrt <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/FBA/Merged/Analyse/SampleInfo.txt",row.names = 1,header = TRUE)

chrt_anno = HeatmapAnnotation(df = chrt, height =unit(0.1, "mm"),annotation_name_side = "left",annotation_name_gp = gpar(fontsize = 0),
                              annotation_legend_param  = list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),
                                                              title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
                              col = list(Cohort=c("HC"="#808080","EC"="#1919ff","ART"="#FFA500")))


col_fun1 = colorRamp2(c(1000,500,10, 0,-10,-500,-1000), c("#7F7F00","#B2B200" ,"#E5E500","white","#BF7FBF","#993299","#590059"))


H1=Heatmap(as.matrix((FBA)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,na_col = "white",
           row_title_gp =gpar(fontsize = 0),column_split = chrt$Cohort,left_annotation = sub_anno,right_annotation = ha,
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),column_names_gp =gpar(fontsize = 12),column_dend_height =unit(3, "cm"), 
           top_annotation = chrt_anno,clustering_distance_columns = "euclidean",row_split = sub$SubSystem,
           row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),heatmap_legend_param = list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),
                                                                                          title_gp = gpar(fontsize = 20), 
                                                                                          labels_gp = gpar(fontsize = 20)),
           name = "Flux",show_row_names = FALSE,row_names_gp=gpar(fontsize = 15),height  = unit(40, "cm"),width  = unit(25, "cm"))


pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/FBA/Merged/Analyse/FBA_Heatmap.pdf",height = 25,width =25)
draw(H1, merge_legend = TRUE)
dev.off()

#########################

FBA=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/GroupFBA/FBA/Merged/Figures/OtherComp.txt",row.names = 1)
sub=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/GroupFBA/FBA/Merged/Figures/OtherSubs.txt",row.names = 1)

sub_anno = rowAnnotation(df = sub, width =unit(0.2, "mm"),annotation_name_gp = gpar(fontsize = 0),
                         annotation_legend_param  = list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),
                                                         title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
                         col = list(SubSystem = c("Alanine, aspartate and glutamate metabolism"="#303964","Amino sugar and nucleotide sugar metabolism"="#596083",
                                                  "Carnitine shuttle (cytosolic)"="#979cb1","Fatty acid activation (cytosolic)"="#d5d7e0",
                                                  "Fatty acid oxidation"="#9c1d5d","Glycolysis / Gluconeogenesis"="#df2a85",
                                                  "Nucleotide metabolism"="#eb7fb5","Purine metabolism"="#f5bfda",
                                                  "Pyrimidine metabolism"="#004c00","Retinol metabolism"="#008000",
                                                  "Sphingolipid metabolism"="#4ca64c","TCA Cycle"="#99cc99","Tryptophan metabolism"="#E1CDB6")))


chrt <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/GroupFBA/FBA/Merged/Figures/SampleInfo.txt",row.names = 1,header = TRUE)

chrt_anno = HeatmapAnnotation(df = chrt, height =unit(0.1, "mm"),annotation_name_side = "left",annotation_name_gp = gpar(fontsize = 0),
                              annotation_legend_param  = list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),
                                                              title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
                              col = list(Cohort=c("HC"="#808080","EC"="#1919ff","ART"="#FFA500")))


col_fun1 = colorRamp2(c(1000,500,10, 0,-10,-500,-1000), c("#7F7F00","#B2B200" ,"#E5E500","white","#BF7FBF","#993299","#590059"))


H1=Heatmap(as.matrix((FBA)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,na_col = "white",
           row_title_gp =gpar(fontsize = 0),column_split = chrt$Cohort,left_annotation = sub_anno,
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),column_names_gp =gpar(fontsize = 12),column_dend_height =unit(3, "cm"), 
           top_annotation = chrt_anno,clustering_distance_columns = "euclidean",row_split = sub$SubSystem,
           row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),heatmap_legend_param = list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),
                                                                                          title_gp = gpar(fontsize = 20), 
                                                                                          labels_gp = gpar(fontsize = 20)),
           name = "Flux",show_row_names = TRUE,row_names_gp=gpar(fontsize = 20),height  = unit(40, "cm"),width  = unit(7, "cm"))


pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/GroupFBA/FBA/Merged/Figures/OtherCompartment.pdf",height = 25,width =35)
draw(H1, merge_legend = TRUE,annotation_legend_side = "bottom",heatmap_legend_side = "left")
dev.off()


###########


FBA=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/GroupFBA/FBA/Merged/Figures/Transport.txt",row.names = 1)


chrt <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/GroupFBA/FBA/Merged/Figures/SampleInfo.txt",row.names = 1,header = TRUE)

chrt_anno = HeatmapAnnotation(df = chrt, height =unit(0.1, "mm"),annotation_name_side = "left",annotation_name_gp = gpar(fontsize = 0),
                              annotation_legend_param  = list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),
                                                              title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
                              col = list(Cohort=c("HC"="#808080","EC"="#1919ff","ART"="#FFA500")))


col_fun1 = colorRamp2(c(1000,500,10, 0,-10,-500,-1000), c("#7F7F00","#B2B200" ,"#E5E500","white","#BF7FBF","#993299","#590059"))


H1=Heatmap(as.matrix((FBA)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,na_col = "white",
           row_title_gp =gpar(fontsize = 0),column_split = chrt$Cohort,
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),column_names_gp =gpar(fontsize = 12),column_dend_height =unit(3, "cm"), 
           top_annotation = chrt_anno,clustering_distance_columns = "euclidean",
           row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),heatmap_legend_param = list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),
                                                                                          title_gp = gpar(fontsize = 20), 
                                                                                          labels_gp = gpar(fontsize = 20)),
           name = "Flux",show_row_names = TRUE,row_names_gp=gpar(fontsize = 20),height  = unit(30, "cm"),width  = unit(7, "cm"))


pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/GroupFBA/FBA/Merged/Figures/Transport.pdf",height = 25,width =35)
draw(H1, merge_legend = TRUE,annotation_legend_side = "bottom",heatmap_legend_side = "left")
dev.off()
##############################


FBA=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/GroupFBA/FBA/Merged/Figures/Mitochondria.txt",row.names = 1)
sub=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/GroupFBA/FBA/Merged/Figures/MitoSubs.txt",row.names = 1)

sub_anno = rowAnnotation(df = sub, width =unit(0.2, "mm"),annotation_name_gp = gpar(fontsize = 0),
                         annotation_legend_param  = list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),
                                                         title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
                         col = list(SubSystem = c("Alanine, aspartate and glutamate metabolism"="#303964","Amino sugar and nucleotide sugar metabolism"="#596083",
                                                  "Arginine and proline metabolism"="#979cb1","Bile acid biosynthesis"="#d5d7e0",
                                                  "Carnitine shuttle (cytosolic)"="#9c1d5d","Carnitine shuttle (mitochondrial)"="#df2a85",
                                                  "Fatty acid activation (cytosolic)"="#eb7fb5","Fatty acid oxidation"="#f5bfda",
                                                  "Glycolysis / Gluconeogenesis"="#004c00","Nucleotide metabolism"="#008000",
                                                  "Oxidative phosphorylation"="#4ca64c","Purine metabolism"="#99cc99","Pyrimidine metabolism"="#9d8f7f",
                                                  "Pyruvate metabolism"="#e1cdb6","Retinol metabolism"="#f0e6da","Sphingolipid metabolism"="#586361","TCA Cycle"="#93a6a3",
                                                  "Tryptophan metabolism"="#bec9c7","Tyrosine metabolism"="#e9edec")))


chrt <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/GroupFBA/FBA/Merged/Figures/SampleInfo.txt",row.names = 1,header = TRUE)

chrt_anno = HeatmapAnnotation(df = chrt, height =unit(0.1, "mm"),annotation_name_side = "left",annotation_name_gp = gpar(fontsize = 0),
                              annotation_legend_param  = list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),
                                                              title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
                              col = list(Cohort=c("HC"="#808080","EC"="#1919ff","ART"="#FFA500")))


col_fun1 = colorRamp2(c(1000,500,10, 0,-10,-500,-1000), c("#7F7F00","#B2B200" ,"#E5E500","white","#BF7FBF","#993299","#590059"))


H1=Heatmap(as.matrix((FBA)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,na_col = "white",
           row_title_gp =gpar(fontsize = 0),column_split = chrt$Cohort,left_annotation = sub_anno,
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),column_names_gp =gpar(fontsize = 12),column_dend_height =unit(3, "cm"), 
           top_annotation = chrt_anno,clustering_distance_columns = "euclidean",row_split = sub$SubSystem,
           row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),heatmap_legend_param = list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),
                                                                                          title_gp = gpar(fontsize = 20), 
                                                                                          labels_gp = gpar(fontsize = 20)),
           name = "Flux",show_row_names = TRUE,row_names_gp=gpar(fontsize = 20),height  = unit(44, "cm"),width  = unit(7, "cm"))


pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/GEM/personalized_GEM/GroupFBA/FBA/Merged/Figures/Mitochondria.pdf",height = 25,width =45)
draw(H1, merge_legend = TRUE,annotation_legend_side = "bottom",heatmap_legend_side = "left")
dev.off()


#################### NEW EPIC bar graph
library(reshape2)
library(ggplot2)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NEW_EPIC/MeanFrac.txt")
MM=melt(data)

head(MM)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/NEW_EPIC/CellFrac.pdf")
ggplot(data=MM, aes(fill=Sample, y=value, x=variable)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = c("classical_monocyte"="#b29600","eosinophil"="#e5c100","gdT-cell"="#ffdf32","intermediate_monocyte"="#ffef99","MAIT_T-cell"="#005900","memory_B-cell"="#008000",
                               "memory_CD4_T-cell"="#329932","memory_CD8_T-cell"="#7fbf7f","myeloid_DC"="#009999","naive_B-cell"="#00cccc","naive_CD4_T-cell"="#00ffff",
                               "naive_CD8_T-cell"="#b2ffff","neutrophil"="#660066","NK-cell"="#843284","non-classical_monocyte"="#b27fb2","otherCells"="#e0cce0","plasmacytoid_DC"="#088da5","T-reg"="#6abac9"))+
  theme(axis.title.x = element_blank(),axis.title= element_blank(),legend.title = element_blank(),legend.text = element_text(size = 10),
        plot.margin = margin(4,3,4,3, "cm"),
        axis.text= element_text(color="black",size=10))
dev.off()

#############

library(circlize)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/BatchCorrected/NEW/sampleInfo.txt",row.names = 1,header = T)
geneinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/BatchCorrected/NEW/GeneInfo_Sorted.txt",header = T,row.names = 1)

input=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/BatchCorrected/NEW/LFC.txt",header = T,check.names = FALSE,row.names = 1)

DAT=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/BatchCorrected/NEW/Combat_Sorted.txt",check.names = FALSE,header = T,row.names = 1)
col_fun = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))
library(ComplexHeatmap)
nrow(t(X$carpet))
colours <- list("Cohort"=c("ART"="#FFA500","EC"="#1919ff"))
comp<- list("Complex"=c("Complex1"="#16605f","Complex2"="#891b69","Complex3"="#e5ba33","Complex4"="#332a6c","Complex5"="#d17294"))
col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7F7F00","#B2B200" ,"#E5E500","white","#BF7FBF","#993299","#590059"))
H1=Heatmap(as.matrix(t(DAT)),col=col_fun1,cluster_rows=FALSE,cluster_columns = FALSE,left_annotation = rowAnnotation(Cohort = sampleinfo$Group,col=colours,show_legend=FALSE,show_annotation_name=FALSE),
           top_annotation =columnAnnotation(Complex = geneinfo$Complex,col=comp,show_legend=FALSE,show_annotation_name=FALSE),
           name = "Z-Score",show_row_names = FALSE,row_names_gp=gpar(fontsize = 5),row_order = order(rownames(t(DAT))),
           row_split=c(rep("ART",19),rep("EC",19)),column_split = geneinfo$Complex,height=unit(7, "cm"),width=unit(23, "cm"))

ha = columnAnnotation(foo = anno_mark(at = c(108,28,121,39,101,59,27,40,69,100,95,71,94,93,67,18,60,56,62,73,79,126,15,41,83,92,4,77,96,8,53,75,5,6,80,84,2,1,7,23,3),
                                   labels_gp = gpar(fontsize=10),link_height=unit(5, "mm"),labels = c("ATP6V0C","NDUFB4","ATP6V1G1","NDUFS3","ATP5PF",
                                                                                                                            "UQCRFS1","NDUFB3","NDUFS4","COX5A","ATP5PD",
                                                                                                                            "ATP5MC3","COX6A1","ATP5MC2","ATP5MC1","COX4I1",
                                                                                                                            "NDUFA6","UQCRH","UQCRB","UQCRQ","COX6B1","COX7B",
                                                                                                                            "MT-ATP6","NDUFA4","NDUFS5","MT-CO1","ATP5F1E",
                                                                                                                            "MT-ND4","COX7A2","ATP5ME","NDUFA1","MT-CYB",
                                                                                                                            "COX6C","MT-ND4L","MT-ND5","COX7C","MT-CO2",
                                                                                                                            "MT-ND2","MT-ND1","MT-ND6","NDUFB1","MT-ND3")))
H2=Heatmap(as.matrix((input)),col=col_fun,cluster_rows=FALSE,cluster_columns = FALSE,name="Log2FoldChange",row_names_gp=gpar(fontsize = 12),show_row_names = FALSE,show_column_names = FALSE,
           height=unit(0.25, "cm"),column_names_gp =gpar(fontsize = 8),top_annotation  = ha) 

tt=H1 %v% H2
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/BatchCorrected/NEW/HeatMap.pdf",width = 20,height=7)
draw(tt, merge_legend = TRUE)
dev.off()

H2=Heatmap(as.matrix((input)),col=col_fun,cluster_rows=FALSE,cluster_columns = FALSE,name="Log2FoldChange",row_names_gp=gpar(fontsize = 12),show_row_names = FALSE,show_column_names = FALSE,
           height=unit(0.25, "cm"),column_names_gp =gpar(fontsize = 8),bottom_annotation = ha) 

tt=H1 %v% H2
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/OxPhos/BatchCorrected/NEW/HeatMap2.pdf",width = 20,height=7)
draw(tt, merge_legend = TRUE)
dev.off()

