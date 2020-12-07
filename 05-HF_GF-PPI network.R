##### 选取全部的DEG进行PPI作图######
###### 对照标准格式的ID转换获取 ######
rm(list=ls())
library(clusterProfiler)
library(org.Mm.eg.db)
options(stringsAsFactors = F)
dfHF <- read.table("./string_interactions HvsF.tsv", header=F, sep="\t", fileEncoding="windows-1252") #读取TSV格式的情况
dfGF <- read.table("./string_interactions GvsF.tsv", header=F, sep="\t", fileEncoding="windows-1252") 
dfHF=dfHF[,c(1:4,13)]
dfGF=dfGF[,c(1:4,13)]
Name=c("Node1","Node2","AltNode1Name","AltNode2Name","Combined Score")
colnames(dfHF)=Name
colnames(dfGF)=Name

###### 生成edge文件  ######
edges_GF=dfGF[,c(1:2)]
colnames(edges_GF)=c("FromNode","ToNode")
write.table(edges_GF,file = "edges_GF.txt",col.names=T,row.names = F,sep = "\t",quote = F)
geneGF=c(as.character(edges_GF$FromNode),as.character(edges_GF$ToNode)) 
View(geneGF)


##### #生成node文件 ######
num= table(geneGF)
View(num)
write.table(num,file = "./GF_num.txt",row.names = F,col.names = T,sep = "\t")
num= read.table("./GF_num.txt",header = T,sep = "\t")
colnames(num)=c("symbol","numInteractions")

typeGF =read.csv("~/Documents/MNC-ceRNA数据分析/MNC-ceRNA/venn/transcript/G vs F_DEgene_voom_exp.csv",sep = ",",header = T)
typeGF=typeGF[,c(3,9)]

#匹配对应的上下游符号
typeGFu=typeGF[with(typeGF,typeGF$logFC>1),]
typeGFu$regulation=c(rep("up",length(typeGFu$logFC)))
typeGFd=typeGF[with(typeGF,typeGF$logFC<(-1)),]
typeGFd$regulation=c(rep("down",length(typeGFd$logFC)))
typeGF=rbind(typeGFd,typeGFu)
write.table(typeGF,file="./typeGF.txt",sep = "\t",row.names = F)

num$type= typeGF$regulation[match(num$symbol,typeGF$external_gene_name)]
View(num)
num$type[1]="down"
write.table(num,file="nodes_GF.txt",col.names = T,row.names = F,sep = "\t",quote=F)


##### 选取combined score>0.5进行PPI的作图(HF系列) #####
###### 生成edge文件  ######
dfHF=dfHF[with(dfHF,dfHF$`Combined Score` > 0.8),]
edges_HF=dfHF[,c(1:2)]
colnames(edges_HF)=c("FromNode","ToNode")
write.table(edges_HF,file = "./edges_HF.txt",col.names=T,row.names = F,sep = "\t",quote = F)

##### #生成node文件 ######
geneHF=c(as.character(edges_HF$FromNode),as.character(edges_HF$ToNode)) 
View(geneHF)
numHF = table(geneHF)
View(numHF)
write.table(numHF,file = "./HF_num.txt",row.names = F,col.names = T,sep = "\t")
numHF= read.table("./HF_num.txt",header = T,sep = "\t")
colnames(numHF)= c("symbol","numInteractions")
numHF =numHF [with(numHF,numHF$numInteractions >0),]

### 整合合并所有的数据类型  ###
# voom类型没有显示相应的指标情况，需要先进行转换，其余的再进行相关指标的分析 
typeHF1 =read.csv("~/Documents/MNC-ceRNA数据分析/MNC-ceRNA/venn/transcript/H vs F_DEgene_voom_exp.csv",sep = ",",header = T)
typeHF1=typeHF1[,c(3,9)]
typeHF1u=typeHF1[with(typeHF1,typeHF1$logFC>1),]
typeHF1u$regulation=c(rep("up",length(typeHF1u$logFC)))
typeHF1d=typeHF1[with(typeHF1,typeHF1$logFC<(-1)),]
typeHF1d$regulation=c(rep("down",length(typeHF1d$logFC)))
typeHF1=rbind(typeHF1d,typeHF1u)
typeHF1=typeHF1[,c(3,2)]
colnames(typeHF1)= c("sig","external_gene_name")

typeHF2 =read.csv("~/Documents/MNC-ceRNA数据分析/MNC-ceRNA/venn/transcript/H vs F_DEgene_trend_exp.csv",sep = ",",header = T)
typeHF2=typeHF2[,c(9,10)]

typeHF3 =read.csv("~/Documents/MNC-ceRNA数据分析/MNC-ceRNA/venn/transcript/H vs F_DEgene_DESeq2_exp.csv",sep = ",",header = T)
typeHF3=typeHF3[,c(15,16)]

typeHF4 =read.csv("~/Documents/MNC-ceRNA数据分析/MNC-ceRNA/venn/transcript/H vs F_DEgene_edgeR_exp.csv",sep = ",",header = T)
typeHF4=typeHF4[,c(7,8)]

typeall=unique(rbind(typeHF1,typeHF2,typeHF3,typeHF4))
numHF$type= typeall$sig[match(numHF$symbol,typeall$external_gene_name)]
View(numHF)
numHFx=numHF[with(numHF,is.na(numHF$type)==F),]
numHFn=numHF[with(numHF,is.na(numHF$type)==T),]
numHFn$type= c(rep("up",length(numHFn$type)))
numHF=rbind(numHFn,numHFx)
write.table(numHF,file="./nodes_HF.txt",col.names = T,row.names = F,sep = "\t",quote=F)



## 统计相关的结果可以采用在线分析工具venn diagram 或者采用excel中的函数
###=IF(AND(COUNTIF(B:B,A2)>0,COUNTIF(C:C,A2)>0),A2,"")



