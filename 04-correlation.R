rm(list=ls())
options(stringsAsFactors = F)
library(readxl)
require(stats)
library(stringr)
library(Hmisc)

loc_cod="~/Documents/MNC-ceRNA数据分析/公司数据合集/one element/deg/"
loc_crc="~/Documents/MNC-ceRNA数据分析/公司数据合集/one element/deg_circ/"
loc="~/Documents/MNC-ceRNA数据分析/MNC-ceRNA/correlation/"

GFname=c("BPD1","BPD2","BPD3","MSC1","MSC2","MSC3")
HFname=c("MNC1","MNC2","MNC3","BPD1","BPD2","BPD3")

setwd("~/Documents/MNC-ceRNA数据分析/MNC-ceRNA/correlation/")


##### 读入相关的基因列表 #####
GF_mRNA=readxl::read_xlsx(path=paste0(loc_cod,"F_vs_G/deg.up_down.anno.mRNA.xlsx"),sheet = 1,col_names = T)
GF_mRNA=GF_mRNA[,c(1:5,13)]

GF_lncRNA=readxl::read_xlsx(path=paste0(loc_cod,"F_vs_G/deg.up_down.anno.lncRNA.xlsx"),sheet = 1,col_names = T)
GF_lncRNA=GF_lncRNA[,c(1:5,13)]

GF_circRNA=readxl::read_xlsx(path=paste0(loc_crc,"F_vs_G/deg.up_down.anno.xlsx"),sheet = 1,col_names = T)
GF_circRNA=GF_circRNA[,c(1:4,9)]

HF_mRNA=readxl::read_xlsx(path=paste0(loc_cod,"F_vs_H/deg.up_down.anno.mRNA.xlsx"),sheet = 1,col_names = T)
HF_mRNA=HF_mRNA[,c(1:5,13)]

HF_lncRNA=readxl::read_xlsx(path=paste0(loc_cod,"F_vs_H/deg.up_down.anno.lncRNA.xlsx"),sheet = 1,col_names = T)
HF_lncRNA=HF_lncRNA[,c(1:5,13)]

HF_circRNA=readxl::read_xlsx(path=paste0(loc_crc,"F_vs_H/deg.up_down.anno.xlsx"),sheet = 1,col_names = T)
HF_circRNA=HF_circRNA[,c(1:4,9)]



##### mRNA 数据导入 #####
mRNAgf=readxl::read_xlsx(path=paste0(loc_cod,"F_vs_G/deg.up_down.anno.mRNA.xlsx"),sheet = 1,col_names = T)
print(colnames(mRNAgf))
mRNAgf=mRNAgf[,c(1,15:20)]
mRNAgf=as.data.frame(mRNAgf)
mRNAgf=unique(mRNAgf)
mRNA_gf=mRNAgf[,c(2:7)]
dir.create("./CBC")
write.csv(mRNA_gf,file=paste0(loc,"CBC/mRNAgf.csv"),quote = F,row.names = F)
mRNA_gf=read.csv(file = "./CBC/mRNAgf.csv",sep = ",",row.names = mRNAgf$ID)
mRNA_gf=mRNA_gf[which(rowSums(mRNA_gf==0)==0),] #删除任意行出现一个数值为零的情况
#*** 去掉只要有一列为零的行
#data[which(rowSums(data==0)==0),]
#>str(data)
#data : int ...  ### data的输出结果为数值类型
#若要把持数据类型不变的话，修改如下：
#data = [which(rowSums(data) > 0),  ,drop=FALSE]
mRNA_gf_z<-t(scale(t(mRNA_gf)))

mRNAhf=readxl::read_xlsx(path=paste0(loc_cod,"F_vs_H/deg.up_down.anno.mRNA.xlsx"),sheet = 1,col_names = T)
mRNAhf=mRNAhf[,c(1,15:20)]
mRNAhf=as.data.frame(mRNAhf)
mRNAhf=unique(mRNAhf)
mRNA_hf=mRNAhf[,c(2:7)]
dir.create("./CBC")
write.csv(mRNA_hf,file=paste0(loc,"CBC/mRNAhf.csv"),quote = F,row.names = F)
mRNA_hf=read.csv(file = "./CBC/mRNAhf.csv",sep = ",",row.names = mRNAhf$ID)
mRNA_hf=mRNA_hf[which(rowSums(mRNA_hf==0)==0),] #删除任意行出现一个数值
mRNA_hf_z<-t(scale(t(mRNA_hf)))


##### lncRNA 数据导入 #####
lncRNAgf=readxl::read_xlsx(path=paste0(loc_cod,"F_vs_G/deg.up_down.anno.lncRNA.xlsx"),sheet = 1,col_names = T)
lncRNAgf=lncRNAgf[,c(1,15:20)]
lncRNAgf=as.data.frame(lncRNAgf)
lncRNAgf=unique(lncRNAgf)
lncRNA_gf=lncRNAgf[,c(2:7)]
dir.create("./CBC")
write.csv(lncRNA_gf,file=paste0(loc,"CBC/lncRNAgf.csv"),quote = F,row.names = F)
lncRNA_gf=read.csv(file = "./CBC/lncRNAgf.csv",sep = ",",row.names = lncRNAgf$ID)
lncRNA_gf=lncRNA_gf[which(rowSums(lncRNA_gf==0)==0),]
lncRNA_gf_z<-t(scale(t(lncRNA_gf)))

lncRNAhf=readxl::read_xlsx(path=paste0(loc_cod,"F_vs_H/deg.up_down.anno.lncRNA.xlsx"),sheet = 1,col_names = T)
#print(colnames(lncRNA_hf))
lncRNAhf=lncRNAhf[,c(1,15:20)]
lncRNAhf=as.data.frame(lncRNAhf)
lncRNAhf=unique(lncRNAhf)
lncRNA_hf=lncRNAhf[,c(2:7)]
dir.create("./CBC")
write.csv(lncRNA_hf,file=paste0(loc,"CBC/lncRNAhf.csv"),quote = F,row.names = F)
lncRNA_hf=read.csv(file = "./CBC/lncRNAhf.csv",sep = ",",row.names = lncRNAhf$ID)
lncRNA_hf=lncRNA_hf[which(rowSums(lncRNA_hf==0)==0),] #删除任意行出现一个数值
lncRNA_hf_z<-t(scale(t(lncRNA_hf)))

#统一更换colnames
colnames(mRNA_gf_z)=GFname
colnames(lncRNA_gf_z)=GFname
colnames(mRNA_hf_z)=HFname
colnames(lncRNA_hf_z)=HFname
save(mRNA_gf_z,mRNA_hf_z,lncRNA_gf_z,lncRNA_hf_z,file="./CBC/coding_zscore.Rdata")


##### circRNA 数据导入 #####
circRNAgf=readxl::read_xlsx(path=paste0(loc_crc,"F_vs_G/deg.up_down.anno.xlsx"),sheet = 1,col_names = T)
print(colnames(circRNAgf))
circRNAgf=circRNAgf[,c(1,9,20:25)]
circRNAgf=as.data.frame(circRNAgf)

circRNA_gf=circRNAgf[,c(3:8)]
write.csv(circRNA_gf,file=paste0(loc,"CBC/circRNAgf.csv"),quote = F,row.names = F)
circRNA_gf=read.csv(file = "./CBC/circRNAgf.csv",sep = ",",row.names = circRNAgf$running_number)
circRNA_gf_z<-t(scale(t(circRNA_gf)))



circRNAhf=readxl::read_xlsx(path=paste0(loc_crc,"F_vs_H/deg.up_down.anno.xlsx"),sheet = 1,col_names = T)
print(colnames(circRNAhf))
circRNAhf=circRNAhf[,c(1,9,20:25)]
circRNAhf=as.data.frame(circRNAhf)

circRNA_hf=circRNAhf[,c(3:8)]
write.csv(circRNA_hf,file=paste0(loc,"CBC/circRNAhf.csv"),quote = F,row.names = F)
circRNA_hf=read.csv(file = "./CBC/circRNAhf.csv",sep = ",",row.names = circRNAhf$running_number)
circRNA_hf_z<-t(scale(t(circRNA_hf)))

circRNAhf2=circRNAhf[circRNAhf$circBase!="-",]
circRNA_hf2=circRNAhf2[,c(3:8)]
write.csv(circRNA_hf2,file=paste0(loc,"CBC/circRNAhf2.csv"),quote = F,row.names = F)
circRNA_hf2=read.csv(file = "./CBC/circRNAhf2.csv",sep = ",",row.names = circRNAhf2$circBase)
circRNA_hf2_z<-t(scale(t(circRNA_hf2)))
rownames(circRNA_hf2_z)=str_sub(rownames(circRNA_hf2_z),start = 1,end = -10L)



#删除没有靶向到circbase的情况
circRNAhf2=circRNAhf[circRNAhf$circBase!="-",]
circRNA_hf2=circRNAhf2[,c(3:8)]
write.csv(circRNA_hf2,file=paste0(loc,"CBC/circRNAhf2.csv"),quote = F,row.names = F)
circRNA_hf2=read.csv(file = "./CBC/circRNAhf2.csv",sep = ",",row.names = circRNAhf2$circBase)
circRNA_hf2_z<-t(scale(t(circRNA_hf2)))
rownames(circRNA_hf2_z)=str_sub(rownames(circRNA_hf2_z),start = 1,end = -10L)

#统一的更换colnames
colnames(circRNA_gf_z)=GFname
colnames(circRNA_gf2_z)=GFname
colnames(circRNA_hf_z)=HFname
colnames(circRNA_hf2_z)=HFname
save(circRNA_gf_z,circRNA_gf2_z,circRNA_hf_z,circRNA_hf2_z,file="./CBC/circRNA_zscore.Rdata")

save(mRNA_gf,lncRNA_gf,mRNA_hf,lncRNA_hf,circRNA_gf,circRNA_hf,circRNA_gf2,circRNAhf2,file = "expression relative quant.Rdata")


#### mRNA_lncRNA network building ####
load("coding_zscore.Rdata")
load("circRNA_zscore.Rdata")
##### r corr GF #####
lmGF=rbind(mRNA_gf_z,lncRNA_gf_z)
lmGF_z=t(lmGF)
result <- cor(lmGF_z, method = "pearson")#参数也可以改成pearson
result[1:8,1:8]
round(result, 8) #保留8位小数，可根据实际的情况调整
result2 <- rcorr(result, type = "pearson")
View(head(result2))
cormat = result2$r
round(cormat,8)
View(head(cormat))
pmat= result2$P
round(pmat,10)
View(head(pmat))
combind_Correlation<-function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame( fromNode = rownames(cormat)[col(cormat)[ut]],
              toNode = rownames(cormat)[row(cormat)[ut]],
              cor =(cormat)[ut],
              p = pmat[ut] )
}
GFl_m = combind_Correlation(result2$r,result2$P)
View(head(GFl_m))

GFl_m=GFl_m[is.na(GFl_m$p) == FALSE,]
index01=GFl_m$fromNode %in% rownames(mRNA_gf_z)==TRUE & GFl_m$toNode %in% rownames(mRNA_gf_z)==TRUE
index02=GFl_m$fromNode %in% rownames(lncRNA_gf_z)==TRUE & GFl_m$toNode %in% rownames(lncRNA_gf_z)==TRUE
## 判断元素是否存在于指定的向量中
GFl_m1=GFl_m[!index01,]
GFl_m2=GFl_m1[!index02,]
GFl_mp=GFl_m2[GFl_m2$p < 1.000000e-03, ]

##### r corr HF #####
lmHF=rbind(mRNA_hf_z,lncRNA_hf_z)
lmHF_z=t(lmHF)
result <- cor(lmHF_z, method = "pearson")#参数也可以改成pearson
result[1:8,1:8]
round(result, 8) #保留8位小数，可根据实际的情况调整library(Hmisc)#加载包
result2 <- rcorr(result, type = "pearson")
View(head(result2))
cormat = result2$r
round(cormat,8)
View(head(cormat))
pmat= result2$P
round(pmat,10)
View(head(pmat))
combind_Correlation<-function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame( fromNode = rownames(cormat)[col(cormat)[ut]],
              toNode = rownames(cormat)[row(cormat)[ut]],
              cor =(cormat)[ut],
              p = pmat[ut] )
}
HFl_m = combind_Correlation(result2$r,result2$P)
View(head(HFl_m))

HFl_m=HFl_m[is.na(HFl_m$p) == FALSE,]
index01=HFl_m$fromNode %in% rownames(mRNA_hf_z)==TRUE & HFl_m$toNode %in% rownames(mRNA_hf_z)==TRUE
index02=HFl_m$fromNode %in% rownames(lncRNA_hf_z)==TRUE & HFl_m$toNode %in% rownames(lncRNA_hf_z)==TRUE
## 判断元素是否存在于指定的向量中
HFl_m1=HFl_m[!index01,]
HFl_m2=HFl_m1[!index02,]
HFl_mp=HFl_m2[HFl_m2$p < 1.000000e-03, ]

#### 保存不同阈值差异的显著相关性分析结果 ####
HFl_m_sig=HFl_mp[abs(HFl_mp$cor) > 0.95 & abs(HFl_mp$cor) < 1, ]
save(HFl_m_sig,file = "HFl_m_sig0.95.Rdata")
HFl_m_sig1= HFl_m_sig[order(abs(HFl_m_sig$cor)[1:1000],decreasing = T),]
save(HFl_m_sig1,file = "HFl_m_sig0.95(1000).Rdata")

##### GF  cor.test sep  #####
dir.create("lncRNA-mRNA/separated",recursive = T)
setwd("./lncRNA-mRNA/separated")
gene_list=rownames(mRNA_gf_z)
sep_Correlation_GF<-function(x) {
  getCorrelationPvalue1 <- function(x){
    pvalue = cor.test(x,gene.wanted.exp,method = "spearman")$p.value
    return(pvalue)
  }
  getCorrelationPvalue2 <- function(x){
    pvalue = cor.test(x,gene.wanted.exp,method = "pearson")$p.value
    return(pvalue)
  }
  getCorrelationValue1 <- function(x){
    cor.value = cor.test(x,gene.wanted.exp,method = "spearman")$estimate
    return(cor.value)
  }
  getCorrelationValue2 <- function(x){
    cor.value = cor.test(x,gene.wanted.exp,method = "pearson")$estimate
    return(cor.value)
  }
  p.value.spearman = apply(lncRNA_gf_z,1,getCorrelationPvalue1)
  p.value.pearson = apply(lncRNA_gf_z,1,getCorrelationPvalue2)
  cor.value.spearman = apply(lncRNA_gf_z,1,getCorrelationValue1)
  cor.value.pearson = apply(lncRNA_gf_z,1,getCorrelationValue2)
  result= data.frame(Gene.cor= gene.wanted,
                     GeneName = names(p.value.spearman),
                     P.value.Spearman = unlist(p.value.spearman),
                     Cor.value.Spearman = unlist(cor.value.spearman),
                     P.value.Pearson = unlist(p.value.pearson),
                     Cor.value.Pearson = unlist(cor.value.pearson))
  #result00 = result[order(abs(result$Cor.value.Spearman)[1:100],decreasing = T),]
  result01=result[result$Cor.value.Spearman> 0.9 & result$Cor.value.Spearman != 1 & result$P.value.Spearman<0.05,]
  result02=result01[result01$Cor.value.Pearson> 0.9 & result01$Cor.value.Pearson != 1 & result01$P.value.Pearson<0.05,]
  write.table(result02,paste(as.character(gene_list[i]),"_","cor",".txt",sep=""),col.names = T,row.names = F,quote=F,sep="\t")
}
for (i in 1:length(rownames(mRNA_gf_z)))
{
  gene.wanted =  gene_list[i]
  gene.wanted.exp = unlist(mRNA_gf_z[rownames(mRNA_gf_z) == gene.wanted, ])
  sep_Correlation(x)
}


setwd("~/Documents/MNC-ceRNA数据分析/MNC-ceRNA/correlation/CBC/lncRNA-mRNA/GF/separated")
filenames <- list.files()
for (file in filenames) {
  if (!exists("all.data")) {
    all.data <- read.table(file, header = TRUE,sep = "\t")
  }
  if (exists("all.data")) {
    new.data <- read.table(file, header = TRUE,sep = "\t")
    all.data <- rbind(all.data, new.data)
  }
}

write.table(all.data,file = paste0(loc,"GF_mRNA_lncRNA_0.9.txt"),row.names=F,sep = "\t",quote = F)

# remain the unduplicated levels
GF_mRNA_lncRNA=read.table(file=paste0(loc,"GF_mRNA_lncRNA_0.9.txt"),sep = "\t",header = T)
GF_mRNA_lncRNA=GF_mRNA_lncRNA[GF_mRNA_lncRNA$Gene.cor != "" & GF_mRNA_lncRNA$GeneName != "", ]
GF_mRNA_lncRNA$Gene.cor.sym=GF_mRNA$symbol[match(GF_mRNA_lncRNA$Gene.cor,GF_mRNA$ID)]
GF_mRNA_lncRNA=GF_mRNA_lncRNA[GF_mRNA_lncRNA$Cor.value.Spearman > 0 & GF_mRNA_lncRNA$Cor.value.Pearson > 0,]
#此步骤可省，将abs参数改正

GF_mRNA_lncRNA$Gene.cor.reg=GF_mRNA$regulation[match(GF_mRNA_lncRNA$Gene.cor,GF_mRNA$ID)]
GF_mRNA_lncRNA$GeneName.reg=GF_mRNA$regulation[match(GF_mRNA_lncRNA$GeneName,GF_lncRNA$ID)]

# 进一步匹配到上下调信息
GF_mRNA_lncRNA=GF_mRNA_lncRNA[is.na(GF_mRNA_lncRNA$Gene.cor.sym)==F &
                                is.na(GF_mRNA_lncRNA$Gene.cor.reg)==F &
                                is.na(GF_mRNA_lncRNA$GeneName.reg)==F ,]


GF_mRNA_lncRNA=GF_mRNA_lncRNA[str_sub(GF_mRNA_lncRNA$GeneName,start=1,end=4)!="ENSM",]
GF_mRNA_lncRNA= GF_mRNA_lncRNA[str_sub(GF_mRNA_lncRNA$GeneName,start  = 1,end = 5) !="MERGE", ]
GF_mRNA_lncRNA= GF_mRNA_lncRNA[str_sub(GF_mRNA_lncRNA$Gene.cor.sym,start  = -3L,end = -1L) !="Rik", ]

#结合生物学相关性进行相关信息的进一步滤除
GF_mRNA_lncRNA=GF_mRNA_lncRNA[GF_mRNA_lncRNA$Gene.cor.reg == GF_mRNA_lncRNA$GeneName.reg, ]
GF_ml_edges=GF_mRNA_lncRNA[,c(1,2,5:9)]
write.table(GF_ml_edges,file=paste0(loc,"GF_ml_edges.txt"),row.names = F,quote = F,sep = "\t")


GF_fromnode=data.frame(table(GF_mRNA_lncRNA$Gene.cor))
GF_tonode=data.frame(table(GF_mRNA_lncRNA$GeneName))
colnames(GF_fromnode)=colnames(GF_tonode)=c("symbol","num")

GF_fromnode$reg= GF_mRNA$regulation[match(GF_fromnode$symbol,GF_mRNA$ID)]
GF_tonode$reg= GF_lncRNA$regulation[match(GF_tonode$symbol,GF_lncRNA$ID)]
GF_node=rbind(GF_fromnode,GF_tonode)
GF_node$type=c(rep("pc",length(GF_fromnode$symbol)),rep("lnc",length(GF_tonode$symbol)))
write.table(GF_node,file=paste0(loc,"GF_ml_nodes.txt"),row.names = F,quote = F,sep = "\t")

## 适当根据网络情况调整
GF_mRNA_lncRNA01=GF_mRNA_lncRNA[GF_mRNA_lncRNA$Cor.value.Pearson>0.97 &GF_mRNA_lncRNA$P.value.Pearson <0.05 ,]
write.table(GF_mRNA_lncRNA01,file=paste0(loc,"GF_ml_edges01.txt"),row.names = F,quote = F,sep = "\t")

GF_fromnode01=data.frame(table(GF_mRNA_lncRNA01$Gene.cor))
GF_tonode01=data.frame(table(GF_mRNA_lncRNA01$GeneName))
colnames(GF_fromnode01)=colnames(GF_tonode01)=c("symbol","num")

GF_fromnode01$reg= GF_mRNA$regulation[match(GF_fromnode01$symbol,GF_mRNA$ID)]
GF_tonode01$reg= GF_lncRNA$regulation[match(GF_tonode01$symbol,GF_lncRNA$ID)]
GF_node01=rbind(GF_fromnode01,GF_tonode01)
GF_node01$type=c(rep("pc",length(GF_fromnode01$symbol)),rep("lnc",length(GF_tonode01$symbol)))
write.table(GF_node01,file=paste0(loc,"GF_ml_nodes01.txt"),row.names = F,quote = F,sep = "\t")



##### HF  cor.test sep  #####
dir.create("./CBC/lncRNA-mRNA/HF/separated",recursive = T)
setwd("./CBC/lncRNA-mRNA/HF/separated")
gene_list=rownames(mRNA_hf_z)
sep_Correlation_HF<-function(x) {
  getCorrelationPvalue1 <- function(x){
    pvalue = cor.test(x,gene.wanted.exp,method = "spearman")$p.value
    return(pvalue)
  }
  getCorrelationPvalue2 <- function(x){
    pvalue = cor.test(x,gene.wanted.exp,method = "pearson")$p.value
    return(pvalue)
  }
  getCorrelationValue1 <- function(x){
    cor.value = cor.test(x,gene.wanted.exp,method = "spearman")$estimate
    return(cor.value)
  }
  getCorrelationValue2 <- function(x){
    cor.value = cor.test(x,gene.wanted.exp,method = "pearson")$estimate
    return(cor.value)
  }
  p.value.spearman = apply(lncRNA_hf_z,1,getCorrelationPvalue1)
  p.value.pearson = apply(lncRNA_hf_z,1,getCorrelationPvalue2)
  cor.value.spearman = apply(lncRNA_hf_z,1,getCorrelationValue1)
  cor.value.pearson = apply(lncRNA_hf_z,1,getCorrelationValue2)
  result= data.frame(Gene.cor= gene.wanted,
                     GeneName = names(p.value.spearman),
                     P.value.Spearman = unlist(p.value.spearman),
                     Cor.value.Spearman = unlist(cor.value.spearman),
                     P.value.Pearson = unlist(p.value.pearson),
                     Cor.value.Pearson = unlist(cor.value.pearson))
  #result00 = result[order(abs(result$Cor.value.Spearman)[1:100],decreasing = T),]
  result01=result[result$Cor.value.Spearman> 0.9 & result$Cor.value.Spearman != 1 & result$P.value.Spearman<0.05,]
  result02=result01[result01$Cor.value.Pearson> 0.9 & result01$Cor.value.Pearson != 1 & result01$P.value.Pearson<0.05,]
  write.table(result02,paste(as.character(gene_list[i]),"_","cor",".txt",sep=""),col.names = T,row.names = F,quote=F,sep="\t")
}
for (i in 1:length(rownames(mRNA_hf_z)))
{
  gene.wanted =  gene_list[i]
  gene.wanted.exp = unlist(mRNA_hf_z[rownames(mRNA_hf_z) == gene.wanted, ])
  sep_Correlation(x)
}


filenames <- list.files()
for (file in filenames) {
  if (!exists("all.data")) {
    all.data <- read.table(file, header = TRUE,sep = "\t")
  }
  if (exists("all.data")) {
    new.data <- read.table(file, header = TRUE,sep = "\t")
    all.data <- rbind(all.data, new.data)
  }
}

write.table(all.data,file = paste0(loc,"HF_mRNA_lncRNA_0.9.txt"),row.names=F,sep = "\t",quote = F)


HF_mRNA_lncRNA=all.data
HF_mRNA_lncRNA$Gene.cor.sym=HF_mRNA$symbol[match(HF_mRNA_lncRNA$Gene.cor,HF_mRNA$ID)]

HF_mRNA_lncRNA=HF_mRNA_lncRNA[HF_mRNA_lncRNA$Cor.value.Spearman > 0 & HF_mRNA_lncRNA$Cor.value.Pearson > 0,]
#此步骤可省，将abs参数改正

HF_mRNA_lncRNA$Gene.cor.reg=HF_mRNA$regulation[match(HF_mRNA_lncRNA$Gene.cor,HF_mRNA$ID)]
HF_mRNA_lncRNA$GeneName.reg=HF_mRNA$regulation[match(HF_mRNA_lncRNA$GeneName,HF_lncRNA$ID)]
# 进一步匹配到上下调信息
HF_mRNA_lncRNA=HF_mRNA_lncRNA[is.na(HF_mRNA_lncRNA$Gene.cor.sym)==F &
                              is.na(HF_mRNA_lncRNA$Gene.cor.reg)==F &
                              is.na(HF_mRNA_lncRNA$GeneName.reg)==F ,]


HF_mRNA_lncRNA=HF_mRNA_lncRNA[str_sub(HF_mRNA_lncRNA$GeneName,start=1,end=4)!="ENSM",]
HF_mRNA_lncRNA= HF_mRNA_lncRNA[str_sub(HF_mRNA_lncRNA$GeneName,start  = 1,end = 5) !="MERGE", ]
HF_mRNA_lncRNA= HF_mRNA_lncRNA[str_sub(HF_mRNA_lncRNA$Gene.cor.sym,start  = -3L,end = -1L) !="Rik", ]
#结合生物学背景进行相关信息的进一步滤除
HF_mRNA_lncRNA=HF_mRNA_lncRNA[HF_mRNA_lncRNA$Gene.cor.reg == HF_mRNA_lncRNA$GeneName.reg, ]
HF_ml_edges=HF_mRNA_lncRNA[,c(1,2,5:9)]
write.table(HF_ml_edges,file=paste0(loc,"HF_ml_edges.txt"),row.names = F,quote = F,sep = "\t")

HF_fromnode=data.frame(table(HF_mRNA_lncRNA$Gene.cor))
HF_tonode=data.frame(table(HF_mRNA_lncRNA$GeneName))
colnames(HF_fromnode)=colnames(HF_tonode)=c("symbol","num")

HF_fromnode$reg= HF_mRNA$regulation[match(HF_fromnode$symbol,HF_mRNA$ID)]
HF_tonode$reg= HF_lncRNA$regulation[match(HF_tonode$symbol,HF_lncRNA$ID)]
HF_node=rbind(HF_fromnode,HF_tonode)
HF_node$type=c(rep("pc",length(HF_fromnode$symbol)),rep("lnc",length(HF_tonode$symbol)))
write.table(HF_node,file=paste0(loc,"HF_ml_nodes.txt"),row.names = F,quote = F,sep = "\t")

## 适当根据网络情况调整
HF_mRNA_lncRNA01=HF_mRNA_lncRNA[HF_mRNA_lncRNA$Cor.value.Pearson>0.95 &HF_mRNA_lncRNA$P.value.Pearson <0.001 ,]
write.table(HF_mRNA_lncRNA01,file=paste0(loc,"HF_ml_edges01.txt"),row.names = F,quote = F,sep = "\t")

HF_fromnode01=data.frame(table(HF_mRNA_lncRNA01$Gene.cor))
HF_tonode01=data.frame(table(HF_mRNA_lncRNA01$GeneName))
colnames(HF_fromnode01)=colnames(HF_tonode01)=c("symbol","num")

HF_fromnode01$reg= HF_mRNA$regulation[match(HF_fromnode01$symbol,HF_mRNA$ID)]
HF_tonode01$reg= HF_lncRNA$regulation[match(HF_tonode01$symbol,HF_lncRNA$ID)]
HF_node01=rbind(HF_fromnode01,HF_tonode01)
HF_node01$type=c(rep("pc",length(HF_fromnode01$symbol)),rep("lnc",length(HF_tonode01$symbol)))
write.table(HF_node01,file=paste0(loc,"HF_ml_nodes01.txt"),row.names = F,quote = F,sep = "\t")




#### mRNA-circRNA network building ####
##### GF  cor.test sep  #####
dir.create("circRNA-mRNA/separated/GF",recursive = T)
setwd("./circRNA-mRNA/separated/GF")
gene_list=rownames(mRNA_gf_z)
sep_Correlation_GF<-function(x) {
  getCorrelationPvalue1 <- function(x){
    pvalue = cor.test(x,gene.wanted.exp,method = "spearman")$p.value
    return(pvalue)
  }
  getCorrelationPvalue2 <- function(x){
    pvalue = cor.test(x,gene.wanted.exp,method = "pearson")$p.value
    return(pvalue)
  }
  getCorrelationValue1 <- function(x){
    cor.value = cor.test(x,gene.wanted.exp,method = "spearman")$estimate
    return(cor.value)
  }
  getCorrelationValue2 <- function(x){
    cor.value = cor.test(x,gene.wanted.exp,method = "pearson")$estimate
    return(cor.value)
  }
  p.value.spearman = apply(circRNA_gf2_z,1,getCorrelationPvalue1)
  p.value.pearson = apply(circRNA_gf2_z,1,getCorrelationPvalue2)
  cor.value.spearman = apply(circRNA_gf2_z,1,getCorrelationValue1)
  cor.value.pearson = apply(circRNA_gf2_z,1,getCorrelationValue2)
  result= data.frame(Gene.cor= gene.wanted,
                     GeneName = names(p.value.spearman),
                     P.value.Spearman = unlist(p.value.spearman),
                     Cor.value.Spearman = unlist(cor.value.spearman),
                     P.value.Pearson = unlist(p.value.pearson),
                     Cor.value.Pearson = unlist(cor.value.pearson))
  #result00 = result[order(abs(result$Cor.value.Spearman)[1:100],decreasing = T),]
  result01=result[result$Cor.value.Spearman> 0.6 & result$Cor.value.Spearman != 1 & result$P.value.Spearman<0.05,]
  result02=result01[result01$Cor.value.Pearson> 0.6 & result01$Cor.value.Pearson != 1 & result01$P.value.Pearson<0.05,]
  write.table(result02,paste(as.character(gene_list[i]),"_","cor",".txt",sep=""),col.names = T,row.names = F,quote=F,sep="\t")
}
for (i in 1:length(rownames(mRNA_gf_z)))
{
  gene.wanted =  gene_list[i]
  gene.wanted.exp = unlist(mRNA_gf_z[rownames(mRNA_gf_z) == gene.wanted, ])
  sep_Correlation_GF(x)
}

filenames <- list.files()
for (file in filenames) {
  if (!exists("all.data")) {
    all.data <- read.table(file, header = TRUE,sep = "\t")
  }
  if (exists("all.data")) {
    new.data <- read.table(file, header = TRUE,sep = "\t")
    all.data <- rbind(all.data, new.data)
  }
}

write.table(all.data,file = paste0(loc,"GF_mRNA_circRNA.txt"),row.names=F,sep = "\t",quote = F)

# remain the unduplicated levels
GF_mRNA_circRNA=read.table(file=paste0(loc,"GF_mRNA_circRNA.txt"),sep = "\t",header = T)
GF_mRNA_circRNA=GF_mRNA_circRNA[GF_mRNA_circRNA$Gene.cor != "" & GF_mRNA_circRNA$GeneName != "", ]
GF_mRNA_circRNA$Gene.cor.sym=GF_mRNA$symbol[match(GF_mRNA_circRNA$Gene.cor,GF_mRNA$ID)]
GF_mRNA_circRNA=GF_mRNA_circRNA[GF_mRNA_circRNA$Cor.value.Spearman > 0 & GF_mRNA_circRNA$Cor.value.Pearson > 0,]
#此步骤可省，将abs参数改正
GF_circRNA$circBase=str_sub(GF_circRNA$circBase,start= 1L,end = -10L)
GF_mRNA_circRNA$Gene.cor.reg=GF_mRNA$regulation[match(GF_mRNA_circRNA$Gene.cor,GF_mRNA$ID)]
GF_mRNA_circRNA$GeneName.reg=GF_mRNA$regulation[match(GF_mRNA_circRNA$GeneName,GF_circRNA$circBase)]

# 进一步匹配到上下调信息
GF_mRNA_circRNA=GF_mRNA_circRNA[is.na(GF_mRNA_circRNA$Gene.cor.sym)==F &
                                is.na(GF_mRNA_circRNA$Gene.cor.reg)==F &
                                is.na(GF_mRNA_circRNA$GeneName.reg)==F ,]


#结合生物学相关性进行相关信息的进一步滤除
GF_mRNA_circRNA=GF_mRNA_circRNA[GF_mRNA_circRNA$Gene.cor.reg == GF_mRNA_circRNA$GeneName.reg, ]
GF_mc_edges=GF_mRNA_circRNA[,c(1,2,5:9)]
write.table(GF_mc_edges,file=paste0(loc,"GF_mc_edges.txt"),row.names = F,quote = F,sep = "\t")


GF_mc_fromnode=data.frame(table(GF_mRNA_circRNA$Gene.cor))
GF_mc_tonode=data.frame(table(GF_mRNA_circRNA$GeneName))
colnames(GF_mc_fromnode)=colnames(GF_mc_tonode)=c("symbol","num")

GF_mc_fromnode$reg= GF_mRNA$regulation[match(GF_mc_fromnode$symbol,GF_mRNA$ID)]
GF_mc_tonode$reg= GF_circRNA$regulation[match(GF_mc_tonode$symbol,GF_circRNA$circBase)]
GF_mc_node=rbind(GF_mc_fromnode,GF_mc_tonode)
GF_mc_node$type=c(rep("pc",length(GF_mc_fromnode$symbol)),rep("circ",length(GF_mc_tonode$symbol)))
write.table(GF_mc_node,file=paste0(loc,"GF_mc_nodes.txt"),row.names = F,quote = F,sep = "\t")

#### 适当根据网络情况调整 ####
GF_mRNA_circRNA01=GF_mRNA_circRNA[GF_mRNA_circRNA$Cor.value.Pearson>0.97 &GF_mRNA_circRNA$P.value.Pearson <0.05 ,]
write.table(GF_mRNA_circRNA01,file=paste0(loc,"GF_ml_edges01.txt"),row.names = F,quote = F,sep = "\t")

GF_fromnode01=data.frame(table(GF_mRNA_circRNA01$Gene.cor))
GF_tonode01=data.frame(table(GF_mRNA_circRNA01$GeneName))
colnames(GF_fromnode01)=colnames(GF_tonode01)=c("symbol","num")

GF_fromnode01$reg= GF_mRNA$regulation[match(GF_fromnode01$symbol,GF_mRNA$ID)]
GF_tonode01$reg= GF_circRNA$regulation[match(GF_tonode01$symbol,GF_circRNA$ID)]
GF_node01=rbind(GF_fromnode01,GF_tonode01)
GF_node01$type=c(rep("pc",length(GF_fromnode01$symbol)),rep("lnc",length(GF_tonode01$symbol)))
write.table(GF_node01,file=paste0(loc,"GF_ml_nodes01.txt"),row.names = F,quote = F,sep = "\t")



##### HF  cor.test sep  #####
setwd("~/Documents/MNC-ceRNA数据分析/MNC-ceRNA/correlation/CBC/")
dir.create("circRNA-mRNA/separated/HF",recursive = T)
setwd("circRNA-mRNA/separated/HF")
gene_list=rownames(mRNA_hf_z)
sep_Correlation_HF<-function(x) {
  getCorrelationPvalue1 <- function(x){
    pvalue = cor.test(x,gene.wanted.exp,method = "spearman")$p.value
    return(pvalue)
  }
  getCorrelationPvalue2 <- function(x){
    pvalue = cor.test(x,gene.wanted.exp,method = "pearson")$p.value
    return(pvalue)
  }
  getCorrelationValue1 <- function(x){
    cor.value = cor.test(x,gene.wanted.exp,method = "spearman")$estimate
    return(cor.value)
  }
  getCorrelationValue2 <- function(x){
    cor.value = cor.test(x,gene.wanted.exp,method = "pearson")$estimate
    return(cor.value)
  }
  p.value.spearman = apply(circRNA_hf2_z,1,getCorrelationPvalue1)
  p.value.pearson = apply(circRNA_hf2_z,1,getCorrelationPvalue2)
  cor.value.spearman = apply(circRNA_hf2_z,1,getCorrelationValue1)
  cor.value.pearson = apply(circRNA_hf2_z,1,getCorrelationValue2)
  result= data.frame(Gene.cor= gene.wanted,
                     GeneName = names(p.value.spearman),
                     P.value.Spearman = unlist(p.value.spearman),
                     Cor.value.Spearman = unlist(cor.value.spearman),
                     P.value.Pearson = unlist(p.value.pearson),
                     Cor.value.Pearson = unlist(cor.value.pearson))
  #result00 = result[order(abs(result$Cor.value.Spearman)[1:100],decreasing = T),]
  result01=result[result$Cor.value.Spearman> 0.6 & result$Cor.value.Spearman != 1 & result$P.value.Spearman<0.05,]
  result02=result01[result01$Cor.value.Pearson> 0.6 & result01$Cor.value.Pearson != 1 & result01$P.value.Pearson<0.05,]
  write.table(result02,paste(as.character(gene_list[i]),"_","cor",".txt",sep=""),col.names = T,row.names = F,quote=F,sep="\t")
}
for (i in 1:length(rownames(mRNA_hf_z)))
{
  gene.wanted =  gene_list[i]
  gene.wanted.exp = unlist(mRNA_hf_z[rownames(mRNA_hf_z) == gene.wanted, ])
  sep_Correlation_HF(x)
}


filenames <- list.files()
for (file in filenames) {
  if (!exists("all.data")) {
    all.data <- read.table(file, header = TRUE,sep = "\t")
  }
  if (exists("all.data")) {
    new.data <- read.table(file, header = TRUE,sep = "\t")
    all.data <- rbind(all.data, new.data)
  }
}

write.table(all.data,file = paste0(loc,"HF_mRNA_circRNA.txt"),row.names=F,sep = "\t",quote = F)


HF_mRNA_circRNA=all.data
HF_mRNA_circRNA$Gene.cor.sym=HF_mRNA$symbol[match(HF_mRNA_circRNA$Gene.cor,HF_mRNA$ID)]
#此步骤可省，将abs参数改正
HF_circRNA$circBase=str_sub(HF_circRNA$circBase,start= 1L,end = -10L)
HF_mRNA_circRNA$Gene.cor.reg=HF_mRNA$regulation[match(HF_mRNA_circRNA$Gene.cor,HF_mRNA$ID)]
HF_mRNA_circRNA$GeneName.reg=HF_circRNA$regulation[match(HF_mRNA_circRNA$GeneName,HF_circRNA$circBase)]
# 进一步匹配到上下调信息
HF_mRNA_circRNA=HF_mRNA_circRNA[is.na(HF_mRNA_circRNA$Gene.cor.sym)==F &
                                is.na(HF_mRNA_circRNA$Gene.cor.reg)==F &
                                is.na(HF_mRNA_circRNA$GeneName.reg)==F ,]



#结合生物学背景进行相关信息的进一步滤除
HF_mc_edges=HF_mRNA_circRNA[,c(1,2,5:9)]
HF_mc_edges=unique(HF_mc_edges)
write.table(HF_mc_edges,file=paste0(loc,"HF_mc_edges.txt"),row.names = F,quote = F,sep = "\t")

HF_mc_fromnode=data.frame(table(HF_mc_edges$Gene.cor))
HF_mc_tonode=data.frame(table(HF_mc_edges$GeneName))
colnames(HF_mc_fromnode)=colnames(HF_mc_tonode)=c("symbol","num")

HF_mc_fromnode$reg= HF_mRNA$regulation[match(HF_mc_fromnode$symbol,HF_mRNA$ID)]
HF_mc_tonode$reg= HF_circRNA$regulation[match(HF_mc_tonode$symbol,HF_circRNA$circBase)]
HF_mc_node=rbind(HF_mc_fromnode,HF_mc_tonode)
HF_mc_node$type=c(rep("pc",length(HF_mc_fromnode$symbol)),rep("circ",length(HF_mc_tonode$symbol)))
write.table(HF_mc_node,file=paste0(loc,"HF_mc_nodes.txt"),row.names = F,quote = F,sep = "\t")

## 适当根据网络情况调整

HF_mc_edges01=HF_mc_edges[HF_mc_edges$Cor.value.Pearson> 0.9 &HF_mc_edges$P.value.Pearson <0.05 ,]
write.table(HF_mc_edges01,file=paste0(loc,"HF_mc_edges01.txt"),row.names = F,quote = F,sep = "\t")

HF_mc_fromnode01=data.frame(table(HF_mc_edges01$Gene.cor))
HF_mc_tonode01=data.frame(table(HF_mc_edges01$GeneName))
colnames(HF_mc_fromnode01)=colnames(HF_mc_tonode01)=c("symbol","num")

HF_mc_fromnode01$reg= HF_mRNA$regulation[match(HF_mc_fromnode01$symbol,HF_mRNA$ID)]
HF_mc_tonode01$reg= HF_circRNA$regulation[match(HF_mc_tonode01$symbol,HF_circRNA$circBase)]
HF_mc_node01=rbind(HF_mc_fromnode01,HF_mc_tonode01)
HF_mc_node01$type=c(rep("pc",length(HF_mc_fromnode01$symbol)),rep("circ",length(HF_mc_tonode01$symbol)))
write.table(HF_mc_node01,file=paste0(loc,"HF_mc_nodes01.txt"),row.names = F,quote = F,sep = "\t")



### ceRNA network building ####
