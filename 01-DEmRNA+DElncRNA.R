#### environment settings ####
rm(list=ls())
library("stringr")
library("readxl")
library("DESeq2")
library("limma")
library("edgeR")

#创建保存数据的子文件夹 
dir.create(path="./mRNA+lncRNA/pc",recursive = T)
dir.create(path="./mRNA+lncRNA/lnc",recursive = T)

#设定筛选的阈值 
F.change = 1.00
p.adj = 0.05
P.Val = 0.05

#设定需要通用的函数
expo<-function(x){2^(x)}


#### read in the files ####
trans_fpkm <- readxl::read_xlsx(path = "./count/trans.fpkm.anno.xlsx",sheet = 1)
genes_fpkm <- readxl::read_xlsx(path = "./count/genes.fpkm.anno.xlsx",sheet = 1)
trans_count <- readxl::read_xlsx(path = "./count/trans.rawcount.anno.xlsx",sheet = 1)
genes_count <- readxl::read_xlsx(path = "./count/genes.rawcount.anno.xlsx",sheet = 1)

##### gene-lnc 系列的录入 #####
colnames(genes_fpkm)
lnc_fpkm_7days= genes_fpkm[genes_fpkm$Description =="lncRNA",c(1,14:25,29,30)]
lnc_fpkm_7days=as.data.frame(lnc_fpkm_7days)
rownames(lnc_fpkm_7days)=lnc_fpkm_7days$`#id`
lnc_fpkm_7days=lnc_fpkm_7days[,-1]
View(head(lnc_fpkm_7days))

lnc_count_7days= genes_count[genes_count$Description =="lncRNA",c(1,14:25,29,30)]
lnc_count_7days=as.data.frame(lnc_count_7days)
rownames(lnc_count_7days)=lnc_count_7days$`#id`
lnc_count_7days=lnc_count_7days[,-1]
View(head(lnc_count_7days))

lnc_fpkm_28days= genes_fpkm[genes_fpkm$Description =="lncRNA",c(1:13,29,30)]
lnc_fpkm_28days=as.data.frame(lnc_fpkm_28days)
rownames(lnc_fpkm_28days)=lnc_fpkm_28days$`#id`
lnc_fpkm_28days=lnc_fpkm_28days[,-1]
View(head(lnc_fpkm_28days))


lnc_count_28days= genes_count[genes_count$Description =="lncRNA",c(1:13,29,30)]
lnc_count_28days=as.data.frame(lnc_count_28days)
rownames(lnc_count_28days)=lnc_count_28days$`#id`
lnc_count_28days=lnc_count_28days[,-1]
View(head(lnc_count_28days))

##### gene-pc系列的录入 #####
mRNA_fpkm_7days= genes_fpkm[genes_fpkm$Description !="lncRNA",c(1,14:25,29,30)]
mRNA_fpkm_7days=as.data.frame(mRNA_fpkm_7days)
rownames(mRNA_fpkm_7days)=mRNA_fpkm_7days$`#id`
mRNA_fpkm_7days=mRNA_fpkm_7days[,-1]
View(head(mRNA_fpkm_7days))

mRNA_count_7days= genes_count[genes_count$Description !="lncRNA",c(1,14:25,29,30)]
mRNA_count_7days=as.data.frame(mRNA_count_7days)
rownames(mRNA_count_7days)=mRNA_count_7days$`#id`
mRNA_count_7days=mRNA_count_7days[,-1]
View(head(mRNA_count_7days))

mRNA_fpkm_28days= genes_fpkm[genes_fpkm$Description !="lncRNA",c(1:13,29,30)]
mRNA_fpkm_28days=as.data.frame(mRNA_fpkm_28days)
rownames(mRNA_fpkm_28days)=mRNA_fpkm_28days$`#id`
mRNA_fpkm_28days=mRNA_fpkm_28days[,-1]
View(head(mRNA_fpkm_28days))


mRNA_count_28days= genes_count[genes_count$Description !="lncRNA",c(1:13,29,30)]
mRNA_count_28days=as.data.frame(mRNA_count_28days)
rownames(mRNA_count_28days)=mRNA_count_28days$`#id`
mRNA_count_28days=mRNA_count_28days[,-1]
View(head(mRNA_count_28days))







load("./mRNA+lncRNA/mRNA_7days.Rdata")
load("./mRNA+lncRNA/mRNA_28days.Rdata")



 
save(list=ls(),file = "./mRNA+lncRNA/mlnc.Rdata")


colnames(trans_fpkm)
##### trans-lnc 系列的录入 #####

Tlnc_fpkm_7days= trans_fpkm[trans_fpkm$Description =="lncRNA",c(1,14:25,30,31)]
Tlnc_fpkm_7days=as.data.frame(Tlnc_fpkm_7days)
rownames(Tlnc_fpkm_7days)=Tlnc_fpkm_7days$`#id`
Tlnc_fpkm_7days=Tlnc_fpkm_7days[,-1]
View(head(Tlnc_fpkm_7days))

Tlnc_count_7days= trans_count[trans_count$Description =="lncRNA",c(1,14:25,30,31)]
Tlnc_count_7days=as.data.frame(Tlnc_count_7days)
rownames(Tlnc_count_7days)=Tlnc_count_7days$`#id`
Tlnc_count_7days=Tlnc_count_7days[,-1]
View(head(Tlnc_count_7days))

Tlnc_fpkm_28days= trans_fpkm[trans_fpkm$Description =="lncRNA",c(1:13,30,31)]
Tlnc_fpkm_28days=as.data.frame(Tlnc_fpkm_28days)
rownames(Tlnc_fpkm_28days)=Tlnc_fpkm_28days$`#id`
Tlnc_fpkm_28days=Tlnc_fpkm_28days[,-1]
View(head(Tlnc_fpkm_28days))


Tlnc_count_28days= trans_count[trans_count$Description =="lncRNA",c(1:13,30,31)]
Tlnc_count_28days=as.data.frame(Tlnc_count_28days)
rownames(Tlnc_count_28days)=Tlnc_count_28days$`#id`
Tlnc_count_28days= Tlnc_count_28days[,-1]
View(head(Tlnc_count_28days))

##### trans-pc系列的录入 #####
TmRNA_fpkm_7days= trans_fpkm[trans_fpkm$Description !="lncRNA",c(1,14:25,30,31)]
TmRNA_fpkm_7days=as.data.frame(TmRNA_fpkm_7days)
rownames(TmRNA_fpkm_7days)=TmRNA_fpkm_7days$`#id`
TmRNA_fpkm_7days=TmRNA_fpkm_7days[,-1]
View(head(TmRNA_fpkm_7days))

TmRNA_count_7days= trans_count[trans_count$Description !="lncRNA",c(1,14:25,30,31)]
TmRNA_count_7days=as.data.frame(TmRNA_count_7days)
rownames(TmRNA_count_7days)=TmRNA_count_7days$`#id`
TmRNA_count_7days=TmRNA_count_7days[,-1]
View(head(TmRNA_count_7days))

TmRNA_fpkm_28days= trans_fpkm[trans_fpkm$Description !="lncRNA",c(1:13,30,31)]
TmRNA_fpkm_28days=as.data.frame(TmRNA_fpkm_28days)
rownames(TmRNA_fpkm_28days)=TmRNA_fpkm_28days$`#id`
TmRNA_fpkm_28days=TmRNA_fpkm_28days[,-1]
View(head(TmRNA_fpkm_28days))


TmRNA_count_28days= trans_count[trans_count$Description !="lncRNA",c(1:13,30,31)]
TmRNA_count_28days=as.data.frame(TmRNA_count_28days)
rownames(TmRNA_count_28days)=TmRNA_count_28days$`#id`
TmRNA_count_28days=TmRNA_count_28days[,-1]
View(head(TmRNA_count_28days))





#### DE_mRNA analysis ####
#####  7days  #####
# 创建批量保存文件时的文件参数
Gname<-c("BPD","BPDM","BPDN","Con")

G7cp<-list(compare1="B vs A",
           compare2="C vs A",
           compare3="D vs A"
)

# 设置子文件夹路径
dir.create(path="./mRNA+lncRNA/pc/7days",recursive = T)

#### 3 bpd vs 3 control ####

####  DESeq2 ####
library(DESeq2)
datcount= TmRNA_count_7days
datcount=datcount[str_sub(rownames(datcount),start = 1,end = 5)!= "MERGE",]
datcount<-datcount[,-c(13:14)]
View(head(datcount))

loc1<- "./mRNA+lncRNA/pc/7days/"

for (i in 1:3){
datacount <- cbind(datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==Gname[4]],datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==Gname[i]])
##构建dds对象
datacount <- round(as.matrix(datacount))
##分两组
groups<- c1
Data <- data.frame(condition = as.factor(groups),
                   row.names = colnames(datacount))
dds <-DESeqDataSetFromMatrix(countData = datacount, 
                             colData = Data, design= ~condition)
dds$condition <- relevel(dds$condition, ref = "control")
## 去掉所有条件都没有read的基因
dds <- dds[rowSums(counts(dds))>1,]
##使用DESeq函数预估离散度
dds <-DESeq(dds)
res <- results(dds)
##设定阈值筛选差异基因
res <- res[order(res$padj), ]
resData <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
resData <-resData[order(resData$padj),]
resData[which(res$log2FoldChange > F.change & resData$padj < p.adj),'sig'] <- 'up'
resData[which(resData$log2FoldChange < -F.change & resData$padj < p.adj),'sig'] <- 'down'
resData[which(abs(resData$log2FoldChange) <= F.change| resData$padj >= p.adj),'sig'] <- 'none'

diffData <- subset(resData, sig %in% c('up', 'down'))
diffData = na.omit(diffData)

diffDataUP = diffData[(diffData$padj < p.adj & (diffData$log2FoldChange > F.change)),]
diffDataDOWN = diffData[(diffData$padj < p.adj & (diffData$log2FoldChange < (-F.change))),]

##保存差异转录本
DEG_DESeq2 <- diffData$Row.names
exp_res_Data=resData[,c(8:length(colnames(resData)))]
row.names(exp_res_Data)=resData$Row.names

# 写入DE分析文件
write.csv(exp_res_Data, file =paste0(loc1,G7cp[[i]],"_mRNA_DESeq2_exp.csv"), row.names = T)
write.csv(resData, file =paste0(loc1,G7cp[[i]],"_mRNA_DESeq2.csv"), row.names = T)
write.csv(diffData, file =paste0(loc1,G7cp[[i]],"_mRNA_Diff_DESeq2_ALL.csv"), row.names = T)
write.csv(diffDataUP, file=paste0(loc1,G7cp[[i]],"_mRNA_Diff_DESeq2_UP.csv"), row.names = T)
write.csv(diffDataDOWN, file=paste0(loc1,G7cp[[i]],"_mRNA_Diff_DESeq2_DOWN.csv"), row.names = T)
write.table(DEG_DESeq2,file=paste0(loc1,G7cp[[i]],"_DEmRNA_DESeq2.txt"),sep = "\t",quote = F,row.names = F,col.names = F)
}

#### edgeR ####
library(edgeR)
#读取数据
for (i in 1:3){
exprSet <- cbind(datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==Gname[i]],datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==Gname[4]])
##分为两组
groups <- factor(c(rep("control", 3),rep("bpd", 3)))
#exprSet <- exprSet[rowMeans(exprSet) > 1,]
##构建DGEList对象
genelist <- DGEList(counts = exprSet, group = groups)
##分组
design <- model.matrix(~0+groups)
colnames(design) <- levels(groups)
design
## 数据过滤
###  基因至少在某一些文库的count超过10 ~ 15 才被认为是表达。这一步全靠尝试， 剔除太多就缓缓，剔除太少就严格点。 我们可以简单的对每个基因的raw count进行比较，但是建议用CPM（count-per-million)标准化 后再比较，避免了文库大小的影响。
keep <- rowSums(cpm(genelist) > 0.5 ) >=2
#table(keep)
genelist.filter <- genelist[keep,keep.lib.sizes=FALSE]
##数据normalization
genelist.norm <- calcNormFactors(genelist.filter)
##估计离散度

#当不存在实验设计矩阵(design matrix)的时候，estimateDisp 等价于 estimateCommonDisp 和 estimateTagwiseDisp 
#### x <- estimateCommonDispx
#### x <- estimateTagwiseDisp(x)


#而当给定实验设计矩阵(design matrix)时， estimateDisp 等价于 estimateGLMCommonDisp, estimateGLMTrendedDisp 和 estimateGLMTagwiseDisp。 其中tag与gene同义

##找差异基因，并按阈值进行筛选
genelist.Disp <- estimateDisp(genelist.norm, design, robust = TRUE)
#plotBCV(genelist.Disp)

###### 较为宽松的检验方法#####
list <- exactTest(genelist.Disp) #负二项分布的常规方法
Gene_diff <- topTags(list, n=nrow(exprSet))
Gene_edgeR <- Gene_diff$table
#在总表中对应的标记出来基因具体的上下调
Gene_edgeR[which(Gene_edgeR$logFC > F.change & Gene_edgeR$PValue < P.Val),'sig'] <- 'up'
Gene_edgeR[which(Gene_edgeR$logFC < -F.change & Gene_edgeR$PValue < P.Val),'sig'] <- 'down'
Gene_edgeR[which(abs(Gene_edgeR$logFC) <= F.change| Gene_edgeR$PValue >= P.Val),'sig'] <- 'none'
Gene_diff_edgeR <- subset(Gene_edgeR, sig %in% c('up', 'down'))
DEG_edgeR <-row.names(Gene_diff_edgeR)
edgeR_UP <- Gene_diff_edgeR[Gene_diff_edgeR$logFC > F.change & Gene_diff_edgeR$PValue < P.Val,]
edgeR_DOWN <- Gene_diff_edgeR[Gene_diff_edgeR$logFC < -(F.change) & Gene_diff_edgeR$P.Value < P.Val,]

exp_edgeR <- genelist.norm$counts

# 写入DE分析文件
write.csv(exp_edgeR, file = paste0(loc1,G7cp[[i]],"_mRNA_edgeR_exp.csv"))
write.csv(Gene_edgeR, file = paste0(loc1,G7cp[[i]],"_mRNA_edgeR.csv"))
write.csv(Gene_diff_edgeR, file = paste0(loc1,G7cp[[i]],"_mRNA_Diff_edgeR_ALL.csv"))
write.csv(edgeR_UP, file = paste0(loc1,G7cp[[i]],"_mRNA_Diff_edgeR_UP.csv"))
write.csv(edgeR_DOWN, file = paste0(loc1,G7cp[[i]],"_mRNA_Diff_edgeR_DOWN.csv"))
write.table(DEG_edgeR,file=paste0(loc1,G7cp[[i]],"_DEmRNA_edgeR.txt"),sep = "\t",quote = F,row.names = F,col.names = F)


}



###### 较为严格的检验方法 ######
fit <- glmFit(genelist.Disp, design, robust=TRUE) #拟合NB模型，用于解释生物学和技术性导致的基因特异性变异,

head(fit$coefficients)
cntr.vs.KD <- makeContrasts(bpd - control, levels=design)
res <- glmLRT(fit, contrast=cntr.vs.KD)#这里用的是glmQLFTest而不是glmLRT是因为前面用了glmQLTFit进行拟合，所以需要用QL F-test进行检验。如果前面用的是glmFit，那么对应的就是glmLRT. 作者称QL F-test更加严格。多重试验矫正用的也是BH方法。
ig.edger <- res$table[p.adjust(res$table$PValue, method = "BH") < 0.05, ]
topTags(res,n=10)
is.de <- decideTestsDGE(res)
summary(is.de)
plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")






#### limma ####
library(limma)

#在处理RNA-Seq数据时，raw read count先被转成log2-counts-per-million (logCPM)，然后对mean-variance关系建模。建模有两种方法：

#精确权重法（precision weights）也就是“voom"
#经验贝叶斯先验趋势（empirical Bayes prior trend），也就是”limma-trend“

for (i in 1:3){
 
exprSet=cbind(datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==Gname[4]],datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==Gname[i]]) # 表达数据

#exprSet=exprSet[rowMeans(exprSet)>1,] 

group_list = factor(c(rep("control",3),rep("disease",3)))
design <- model.matrix(~group_list)
colnames(design) <- levels(group)
rownames(design) <- colnames(exprSet)

dge <- DGEList(counts = exprSet,group = group_list)
### filter base  use CPM
#keep <- rowSums(cpm(dge) > 0.5 ) >=2
table(keep)
dge.filter <- dge[keep,keep.lib.sizes=FALSE]
dge.norm <- calcNormFactors(dge.filter,method = "TMM")

##### limma-trend #####    
logCPM <- cpm(dge.norm, log=TRUE, prior.count=3)
fit <- lmFit(logCPM, design)
fit <- eBayes(fit,trend = T)
diff_trend <- topTable(fit, coef=ncol(design),n=Inf) ## 可以取出不显著的值


diff_trend[which(diff_trend$logFC > F.change & diff_trend$P.Value < 0.05),'sig'] <- 'up'
diff_trend[which(diff_trend$logFC < -F.change & diff_trend$P.Value < 0.05),'sig'] <- 'down'
diff_trend[which(abs(diff_trend$logFC) <= F.change| diff_trend$P.Value >= 0.05),'sig'] <- 'none'
DEG_trend<- subset(diff_trend, sig %in% c('up', 'down'))

limma_trend <- row.names(DEG_trend)
DEG_trend_UP <- DEG_trend[with(DEG_trend, logFC > F.change), ]
DEG_trend_DOWN <- DEG_trend[with(DEG_trend, logFC < -F.change), ]

trend_exp<- apply(logCPM,2,expo)



## 保存DE分析结果
write.csv(trend_exp,file=paste0(loc1,G7cp[i],"_mRNA_trend_exp.csv"),row.names=T)
write.csv(diff_trend,file=paste0(loc1,G7cp[i],"_mRNA_trend.csv"),row.names=T)
#deg_all_none = topTable(fit2,adjust.method="none",coef=1,adj.P.Val=0.05,
#lfc=log(2,2),number=5000,sort.by = 'logFC')
write.csv(DEG_trend,file=paste0(loc1,G7cp[i],"_mRNA_Diff_trend_ALL.csv"),row.names = T)
write.csv(DEG_trend_UP,file=paste0(loc1,G7cp[i],"_mRNA_Diff_trend_UP.csv"),row.names = T)
write.csv(DEG_trend_DOWN,file=paste0(loc1,G7cp[i],"_mRNA_Diff_trend_DOWN.csv"),row.names = T)
write.table(limma_trend,file=paste0(loc1,G7cp[[i]],"_DEmRNA_trend.txt"),sep = "\t",quote = F,row.names = F,col.names = F)






##### limma-voom #####
   v <- voom(dge.norm, design,plot=T)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  output <- topTable(fit, coef=2,n=Inf,adjust.method = "none")
  diff_voom = na.omit(output)
  
  diff_voom[which(diff_voom$logFC > F.change & diff_voom$P.Value < P.Val),'sig'] <- 'up'
  diff_voom[which(diff_voom$logFC < -F.change & diff_voom$P.Value < P.Val),'sig'] <- 'down'
  diff_voom[which(abs(diff_voom$logFC) <= F.change| diff_voom$P.Value >= P.Val),'sig'] <- 'none'
  DEG_voom<- subset(diff_voom, sig %in% c('up', 'down'))
  
  DEG_voom_UP<- DEG_voom[with(DEG_voom, logFC >1), ]
  DEG_voom_DOWN<- DEG_voom[with(DEG_voom, logFC < -1), ]
  limma_voom<-row.names(DEG_voom)
  
  voom_exp <- v$E
  
  
  
  # 写入DE分析结果
  write.csv(voom_exp,file=paste0(loc1,G7cp[i],"_mRNA_voom_exp.csv"),row.names = T)
  write.csv(DEG_voom,file=paste0(loc1,G7cp[i],"_mRNA_voom.csv"),row.names = T)
  write.csv(DEG_voom_UP,file=paste0(loc1,G7cp[i],"_mRNA_Diff_voom_UP.csv"),row.names = T)
  write.csv(DEG_voom_DOWN,file=paste0(loc1,G7cp[i],"_mRNA_Diff_voom_DOWN.csv"),row.names = T)
  write.table(limma_voom,file=paste0(loc1,G7cp[[i]],"_DEmRNA_voom.txt"),sep = "\t",quote = F,row.names = F,col.names = F)
  
  
}


##### 9 bpd VS 3 con 进行的分析 #####
loc39<-"./mRNA+lncRNA/pc/7days/bpd_vs_ctrl_7days_39"
Ggroup2=c(c(rep("bpd",9)),c(rep("control",3)))
###### DESeq2 ######
library(DESeq2)
datcount= TmRNA_count_7days
datcount=datcount[str_sub(rownames(datcount),start = 1,end = 5)!= "MERGE",]
datcount<-datcount[,-c(13:14)]
datcount[1:4,1:4]
datacount=datcount[,c(1:12)]
##构建dds对象
  datacount <- round(as.matrix(datacount))
  ##分为control和cold两组
  groups<-Ggroup2
  Data <- data.frame(condition = as.factor(groups))
  rownames(Data) <- colnames(datacount)
  dds <-DESeqDataSetFromMatrix(countData = datacount, 
                               colData = Data, design= ~condition)
  dds$condition <- relevel(dds$condition, ref = "control")
  ## 去掉所有条件都没有read的基因
  dds <- dds[rowSums(counts(dds))>1,]
  ##使用DESeq函数预估离散度
  dds <-DESeq(dds)
  res <- results(dds)
  ##设定阈值筛选差异基因
  res <- res[order(res$padj), ]
 resData <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  resData <-resData[order(resData$padj),]
  resData[which(res$log2FoldChange > F.change & resData$padj < P.Val),'sig'] <- 'up'
  resData[which(resData$log2FoldChange < -F.change & resData$padj < P.Val),'sig'] <- 'down'
  resData[which(abs(resData$log2FoldChange) <= F.change| resData$padj >= P.Val),'sig'] <- 'none'
  
  diffData <- subset(resData, sig %in% c('up', 'down'))
  diffData=na.omit(diffData)
  
  diffDataUP = diffData[(diffData$padj < p.adj & (diffData$log2FoldChange > F.change)),]
  diffDataDOWN = diffData[(diffData$padj < p.adj & (diffData$log2FoldChange < (-F.change))),]
  
  DEG_DESeq2 <- diffData$Row.names
  exp_res_Data=resData[,c(8:length(colnames(resData)))]
  row.names(exp_res_Data)=resData$Row.names
  
  # 写入DE文件
  write.csv(exp_res_Data, file = paste0(loc39,"_mRNA_DESeq2_exp.csv"), row.names = T)
  write.csv(resData, file =paste0(loc39,"_mRNA_DESeq2.csv"), row.names = T)
  write.csv(diffData, file =paste0(loc39,"_Diff_DESeq2_ALL.csv"), row.names = T)
  write.csv(diffDataUP, file =paste0(loc39,"._Diff_DESeq2_UP.csv"), row.names = T)
  write.csv(diffDataDOWN, file =paste0(loc39,"_Diff_DESeq2_DOWN.csv"), row.names = T)
  write.table(DEG_DESeq2,file=paste0(loc39,"_DEmRNA_DESeq2.txt"),sep = "\t",quote = F,row.names = F,col.names = F)
  
 
  
###### edgeR ######
library(edgeR)
#读取数据

##构建DGEList对象
  exprSet <- datacount[,c(10:12,1:9)]
  ##分为两组
  groups <- c(rep("control",3),rep("disease",9))
  ##分为两组
  ##构建DGEList对象
  genelist <- DGEList(counts = exprSet, group = groups)
  ##分组
  design <- model.matrix(~0+groups)
  colnames(design) <- levels(groups)
  design

  keep <- rowSums(cpm(genelist) > 0.5 ) >=2
  table(keep)
  genelist.filter <- genelist[keep,keep.lib.sizes=FALSE]
  ##数据normalization
  genelist.norm <- calcNormFactors(genelist.filter)
  ##估计离散度
  genelist.Disp <- estimateDisp(genelist.norm, design, robust = TRUE)
  plotBCV(genelist.Disp)
  
  ###### 较为宽松的检验方法#####
  list <- exactTest(genelist.Disp) #负二项分布的常规方法
  Gene_diff <- topTags(list, n=nrow(exprSet))
  Gene_edgeR <- Gene_diff$table
  #在总表中对应的标记出来基因具体的上下调
  Gene_edgeR[which(Gene_edgeR$logFC > F.change & Gene_edgeR$PValue < P.Val),'sig'] <- 'up'
  Gene_edgeR[which(Gene_edgeR$logFC < -F.change & Gene_edgeR$PValue < P.Val),'sig'] <- 'down'
  Gene_edgeR[which(abs(Gene_edgeR$logFC) <= F.change| Gene_edgeR$PValue >= P.Val),'sig'] <- 'none'
  Gene_diff_edgeR <- subset(Gene_edgeR, sig %in% c('up', 'down'))
  DEG_edgeR <-row.names(Gene_diff_edgeR)
  edgeR_UP <- Gene_diff_edgeR[Gene_diff_edgeR$logFC > F.change & Gene_diff_edgeR $ PValue < P.Val,]
  edgeR_DOWN <- Gene_diff_edgeR[Gene_diff_edgeR$logFC < -(F.change) & Gene_diff_edgeR$PValue < P.Val,]
  exp_edgeR <- genelist.norm$counts
  
  # 写入DE分析文件
  write.csv(exp_edgeR, file = paste0(loc39,"_mRNA_edgeR_exp.csv"))
  write.csv(Gene_edgeR, file = paste0(loc39,"_mRNA_edgeR.csv"))
  write.csv(Gene_diff_edgeR, file = paste0(loc39,"_mRNA_Diff_edgeR_ALL.csv"))
  write.csv(edgeR_UP, file = paste0(loc39,"_mRNA_Diff_edgeR_UP.csv"))
  write.csv(edgeR_DOWN, file = paste0(loc39,"_mRNA_Diff_edgeR_DOWN.csv"))
  write.table(DEG_edgeR,file=paste0(loc39,"_DEmRNA_edgeR.txt"),sep = "\t",quote = F,row.names = F,col.names = F)
  

###### limma ######
  
  exprSet <- datcount[,c(10:12,1:9)]
  group_list = factor(c(rep("control",3),rep("disease",9)))
  design <- model.matrix(~group_list)
  colnames(design) <- levels(group_list)
  rownames(design) <- colnames(exprSet)
  
  dge <- DGEList(counts = exprSet,group = group_list)
  ### filter base  use CPM
  keep <- rowSums(cpm(dge) > 0.5 ) >=2
  table(keep)
  dge.filter <- dge[keep,keep.lib.sizes=FALSE]
  dge.norm <- calcNormFactors(dge.filter,method = "TMM")
  
  ##### limma-trend #####    
  logCPM <- cpm(dge.norm, log=TRUE, prior.count=3)
  fit <- lmFit(logCPM, design)
  fit <- eBayes(fit,trend = T)
  diff_trend <- topTable(fit, coef=ncol(design),n=Inf) ## 可以取出不显著的值
  
  
  diff_trend[which(diff_trend$logFC > F.change & diff_trend$P.Value < P.Val),'sig'] <- 'up'
  diff_trend[which(diff_trend$logFC < -F.change & diff_trend$P.Value < P.Val),'sig'] <- 'down'
  diff_trend[which(abs(diff_trend$logFC) <= F.change| diff_trend$P.Value >= P.Val),'sig'] <- 'none'
  DEG_trend<- subset(diff_trend, sig %in% c('up', 'down'))
  
  limma_trend <- row.names(DEG_trend)
  DEG_trend_UP <- DEG_trend[with(DEG_trend, logFC > F.change), ]
  DEG_trend_DOWN <- DEG_trend[with(DEG_trend, logFC < -F.change), ]
  
  trend_exp<- apply(logCPM,2,expo)
  
  
  
  ## 保存DE分析结果
  write.csv(trend_exp,file=paste0(loc39,"_mRNA_trend_exp.csv"),row.names=T)
  write.csv(diff_trend,file=paste0(loc39,"_mRNA_trend.csv"),row.names=T)
  #deg_all_none = topTable(fit2,adjust.method="none",coef=1,adj.P.Val=0.05,
  #lfc=log(2,2),number=5000,sort.by = 'logFC')
  write.csv(DEG_trend,file=paste0(loc39,"_mRNA_Diff_trend_ALL.csv"),row.names = T)
  write.csv(DEG_trend_UP,file=paste0(loc39,"_mRNA_Diff_trend_UP.csv"),row.names = T)
  write.csv(DEG_trend_DOWN,file=paste0(loc39,"_mRNA_Diff_trend_DOWN.csv"),row.names = T)
  write.table(limma_trend,file=paste0(loc39,"_DEmRNA_trend.txt"),sep = "\t",quote = F,row.names = F,col.names = F)
  
  
  

  ##### limma-voom #####
  v <- voom(dge.norm, design,plot=T)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  output <- topTable(fit, coef=2,n=Inf,adjust.method = "none")
  diff_voom = na.omit(output)
  
  diff_voom[which(diff_voom$logFC > F.change & diff_voom$P.Value < P.Val),'sig'] <- 'up'
  diff_voom[which(diff_voom$logFC < -F.change & diff_voom$P.Value < P.Val),'sig'] <- 'down'
  diff_voom[which(abs(diff_voom$logFC) <= F.change| diff_voom$P.Value >= P.Val),'sig'] <- 'none'
  DEG_voom<- subset(diff_voom, sig %in% c('up', 'down'))
  
  DEG_voom_UP<- DEG_voom[with(DEG_voom, logFC >1), ]
  DEG_voom_DOWN<- DEG_voom[with(DEG_voom, logFC < -1), ]
  limma_voom<-row.names(DEG_voom)
  
  voom_exp <- v$E
  
  
  
  # 写入DE分析结果
  write.csv(voom_exp,file=paste0(loc39,"_mRNA_voom_exp.csv"),row.names = T)
  write.csv(DEG_voom,file=paste0(loc39,"_mRNA_voom.csv"),row.names = T)
  write.csv(DEG_voom_UP,file=paste0(loc39,"_mRNA_Diff_voom_UP.csv"),row.names = T)
  write.csv(DEG_voom_DOWN,file=paste0(loc39,"_mRNA_Diff_voom_DOWN.csv"),row.names = T)
  write.table(limma_voom,file=paste0(loc39,"_DEmRNA_voom.txt"),sep = "\t",quote = F,row.names = F,col.names = F)
  
#### 28days ####
# 设置批量保存文件的参数 
    name28 =c("control","control","control","msc","msc","msc","mnc","mnc","mnc","bpd","bpd","bpd")
    d28=list(
      name1= c("control","control","control","bpd","bpd","bpd"),
      name2= c("bpd","bpd","bpd","msc","msc","msc"),
      name3= c("bpd","bpd","bpd","mnc","mnc","mnc"),
      name4 = c("bpd","bpd","bpd","control","control","control"),
      name5= c("msc","msc","msc","bpd","bpd","bpd"),
      name6= c("mnc","mnc","mnc","bpd","bpd","bpd")
    )
    G28cp<-list(compare1="F vs E",
                compare2="G vs F",
                compare3="H vs F"
    )
    G28name<-c("Con","BPDM","BPDN","BPD")
#  设置文件保存的位置

#因为在实际筛选的过程中发现由于R2较小，得到的DEG总量明显小于7days 分组，因此采用p.adj/FDR和pValue两种矫正方法.基于G系列的案例，选择将pValue作为参照的对象,可以调整参数

dir.create(path = "./mRNA+lncRNA/pc/28days")
loc2<-("./mRNA+lncRNA/pc/28days")

##### DESeq2 #####
library(DESeq2)
datcount= TmRNA_count_28days
datcount=datcount[str_sub(rownames(datcount),start = 1,end = 5)!= "MERGE",]
datcount=datcount[,c(10:12,4:9,1:3)]
View(head(datcount))

loc28<-"./mRNA+lncRNA/pc/28days/"
for (i in 2:3){
  i=3
  datacount <- cbind(datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==G28name[4]],datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==G28name[i]])
  ##构建dds对象
  datacount <- round(as.matrix(datacount))
  ##分两组
  groups<-d28[[i]]
  Data <- data.frame(condition = as.factor(groups),
                     row.names = colnames(datacount))
  dds <-DESeqDataSetFromMatrix(countData = datacount, 
                               colData = Data, design= ~condition)
  dds$condition <- relevel(dds$condition, ref = "bpd")
  ## 去掉所有条件都没有read的基因
  dds <- dds[rowSums(counts(dds))>1,]
  ##使用DESeq函数预估离散度
  dds <-DESeq(dds)
  res <- results(dds)
  ##设定阈值筛选差异基因
  res <- res[order(res$padj), ]
  resData <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  resData <-resData[order(resData$padj),]
  resData[which(res$log2FoldChange > F.change & resData$pvalue < P.Val),'sig'] <- 'up'
  resData[which(resData$log2FoldChange < -F.change & resData$pvalue < P.Val),'sig'] <- 'down'
  resData[which(abs(resData$log2FoldChange) <= F.change| resData$pvalue >= P.Val),'sig'] <- 'none'
  
  diffData <- subset(resData, sig %in% c('up', 'down'))
  diffData = na.omit(diffData)
  
  diffDataUP = diffData[(diffData$pvalue < P.Val & (diffData$log2FoldChange > F.change)),]
  diffDataDOWN = diffData[(diffData$pvalue < P.Val & (diffData$log2FoldChange < (-F.change))),]
  
  ##保存差异转录本
  DEG_DESeq2 <- diffData$Row.names
  exp_res_Data=resData[,c(8:length(colnames(resData)))]
  row.names(exp_res_Data)=resData$Row.names
  
# 写入DE文件
write.csv(exp_res_Data,file =paste0(loc28,G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_mRNA_DESeq2_exp.csv"), row.names = T)
write.csv(resData, file =paste0(loc28,G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_mRNA_DESeq2.csv"), row.names = T)
write.csv(diffData, file =paste0(loc28,G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_mRNA_Diff_DESeq2_ALL.csv"), row.names = T)
write.csv(diffDataUP, file=paste0(loc28,G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_mRNA_Diff_DESeq2_UP.csv"), row.names = T)
write.csv(diffDataDOWN, file=paste0(loc28,G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_mRNA_Diff_DESeq2_DOWN.csv"), row.names = T)
write.table(DEG_DESeq2,file=paste0(loc28,G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_DEmRNA_DEseq2.txt"),sep = "\t",quote = F,row.names = F,col.names = F)
}

# bpd vs con

i=1
datacount <- cbind(datcount[,str_sub(colnames(datcount),start = 3,end = -2L)== G28name[i]],datcount[,str_sub(colnames(datcount),start = 3,end = -2L)== G28name[4]])
##构建dds对象
datacount <- round(as.matrix(datacount))
##分为control和cold两组
groups<- d28[[i]]
Data <- data.frame(condition = as.factor(groups))
rownames(Data) <- colnames(datacount)
dds <-DESeqDataSetFromMatrix(countData = datacount, 
                             colData = Data, design= ~condition)
dds$condition <- relevel(dds$condition, ref = "control")

## 去掉所有条件都没有read的基因
dds <- dds[rowSums(counts(dds))>1,]
##使用DESeq函数预估离散度
dds <-DESeq(dds)
res <- results(dds)
##设定阈值筛选差异基因
res <- res[order(res$pvalue), ]
resData <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
resData <-resData[order(resData$pvalue),]
resData[which(res$log2FoldChange > F.change & resData$pvalue < P.Val),'sig'] <- 'up'
resData[which(resData$log2FoldChange < -F.change & resData$pvalue < P.Val),'sig'] <- 'down'
resData[which(abs(resData$log2FoldChange) <= F.change| resData$pvalue >= P.Val),'sig'] <- 'none'

diffData <- subset(resData, sig %in% c('up', 'down'))
diffData = na.omit(diffData)

diffDataUP = diffData[(diffData$pvalue < P.Val & (diffData$log2FoldChange > F.change)),]
diffDataDOWN = diffData[(diffData$pvalue < P.Val & (diffData$log2FoldChange < (-F.change))),]

DEG_DESeq2 <- diffData$Row.names
exp_res_Data=resData[,c(8:length(colnames(resData)))]
row.names(exp_res_Data)=resData$Row.names

# 写入DE文件
write.csv(exp_res_Data, file =paste0(loc28,G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_mRNA_DESeq2_exp.csv"), row.names = T)
write.csv(resData, file =paste0(loc28,G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_mRNA_DESeq2.csv"), row.names = T)
write.csv(diffData, file =paste0(loc28,G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_mRNA_Diff_DESeq2_ALL.csv"), row.names = T)
write.csv(diffDataUP, file=paste0(loc28,G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_mRNA_Diff_DESeq2_UP.csv"), row.names = T)
write.csv(diffDataDOWN, file=paste0(loc28,G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_mRNA_Diff_DESeq2_DOWN.csv"), row.names = T)
write.table(DEG_DESeq2,file=paste0(loc28,G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_DEmRNA_DEseq2.txt"),sep = "\t",quote = F,row.names = F,col.names = F)



##### edgeR #####
library(edgeR)
#读取数据
for (i in 2:3){ 
exprSet <- cbind(datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==G28name[4]],datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==G28name[i]])
##分为两组
groups<-d28[[i]]
##分为两组
##构建DGEList对象
genelist <- DGEList(counts = exprSet, group = groups)
##分组
design <- model.matrix(~0+groups)
colnames(design) <- levels(groups)
design

keep <- rowSums(cpm(genelist) > 0.5 ) >=2
table(keep)
genelist.filter <- genelist[keep,keep.lib.sizes=FALSE]
##数据normalization
genelist.norm <- calcNormFactors(genelist.filter)
##估计离散度
genelist.Disp <- estimateDisp(genelist.norm, design, robust = TRUE)
plotBCV(genelist.Disp)

###### 较为宽松的检验方法#####
list <- exactTest(genelist.Disp) #负二项分布的常规方法
Gene_diff <- topTags(list, n=nrow(exprSet))
Gene_edgeR <- Gene_diff$table
#在总表中对应的标记出来基因具体的上下调
Gene_edgeR[which(Gene_edgeR$logFC > F.change & Gene_edgeR$PValue < P.Val),'sig'] <- 'up'
Gene_edgeR[which(Gene_edgeR$logFC < -F.change & Gene_edgeR$PValue < P.Val),'sig'] <- 'down'
Gene_edgeR[which(abs(Gene_edgeR$logFC) <= F.change| Gene_edgeR$PValue >= P.Val),'sig'] <- 'none'
Gene_diff_edgeR <- subset(Gene_edgeR, sig %in% c('up', 'down'))
DEG_edgeR <-row.names(Gene_diff_edgeR)
edgeR_UP <- Gene_diff_edgeR[Gene_diff_edgeR$logFC > F.change & Gene_diff_edgeR $ PValue < P.Val,]
edgeR_DOWN <- Gene_diff_edgeR[Gene_diff_edgeR$logFC < -(F.change) & Gene_diff_edgeR$PValue < P.Val,]
exp_edgeR <- genelist.norm$counts

  write.csv(exp_edgeR, file =paste0(loc28,G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_mRNA_edgeR_exp.csv"), row.names = T)
  write.csv(Gene_edgeR, file =paste0(loc28,G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_mRNA_edgeR.csv"), row.names = T)
  write.csv(Gene_diff_edgeR,file = paste0(loc28,G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_mRNA_Diff_edgeR_ALL.csv"), row.names = T)
  write.csv(edgeR_UP, file =paste0(loc28,G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_mRNA_Diff_edgeR_UP.csv"), row.names = T)
  write.csv(edgeR_DOWN, file =paste0(loc28,G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_mRNA_Diff_edgeR_DOWN.csv"), row.names = T)
  write.table(DEG_edgeR,file =paste0(loc28,G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_DEmRNA_edgeR.txt"), row.names = T)
  

}

# ----- bpd vs control 
i = 1
exprSet <- cbind(datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==G28name[i]],datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==G28name[4]])
##分为两组
groups <- factor(d28[[i]])
groups<-d28[[i]]
##分为两组
##构建DGEList对象
genelist <- DGEList(counts = exprSet, group = groups)
##分组
design <- model.matrix(~0+groups)
colnames(design) <- levels(groups)
design

keep <- rowSums(cpm(genelist) > 0.5 ) >=2
table(keep)
genelist.filter <- genelist[keep,keep.lib.sizes=FALSE]
##数据normalization
genelist.norm <- calcNormFactors(genelist.filter)
##估计离散度
genelist.Disp <- estimateDisp(genelist.norm, design, robust = TRUE)
plotBCV(genelist.Disp)

###### 较为宽松的检验方法#####
list <- exactTest(genelist.Disp) #负二项分布的常规方法
Gene_diff <- topTags(list, n=nrow(exprSet))
Gene_edgeR <- Gene_diff$table
#在总表中对应的标记出来基因具体的上下调
Gene_edgeR[which(Gene_edgeR$logFC > F.change & Gene_edgeR$PValue < P.Val),'sig'] <- 'up'
Gene_edgeR[which(Gene_edgeR$logFC < -F.change & Gene_edgeR$PValue < P.Val),'sig'] <- 'down'
Gene_edgeR[which(abs(Gene_edgeR$logFC) <= F.change| Gene_edgeR$PValue >= P.Val),'sig'] <- 'none'
Gene_diff_edgeR <- subset(Gene_edgeR, sig %in% c('up', 'down'))
DEG_edgeR <-row.names(Gene_diff_edgeR)
edgeR_UP <- Gene_diff_edgeR[Gene_diff_edgeR$logFC > F.change & Gene_diff_edgeR $ PValue < P.Val,]
edgeR_DOWN <- Gene_diff_edgeR[Gene_diff_edgeR$logFC < -(F.change) & Gene_diff_edgeR$PValue < P.Val,]
exp_edgeR <- genelist.norm$counts


# 写入DE文件
write.csv(exp_edgeR, file = paste0(loc28,G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_mRNA_edgeR_exp.csv"), row.names = T)
write.csv(Gene_edgeR, file = paste0(loc28,G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_mRNA_edgeR.csv"), row.names = T)
write.csv(Gene_diff_edgeR, file = paste0(loc28,G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_mRNA_Diff_edgeR.csv"), row.names = T)
write.csv(edgeR_UP, file = paste0(loc28,G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_mRNA_Diff_edgeR_UP.csv"), row.names = T)
write.csv(edgeR_DOWN, file = paste0(loc28,G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_mRNA_Diff_edgeR_DOWN.csv"), row.names = T)
write.table(DEG_edgeR,file = paste0(loc28,G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_DEmRNA_edgeR.txt"),sep = "\t",quote = F,row.names = F,col.names = F)



##### #limma ######
##### limma-trend #####
for (i in 2:3){
  
  exprSet=cbind(datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==Gname[4]],datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==Gname[i]]) # 表达数据
  
  #exprSet=exprSet[rowMeans(exprSet)>1,] 
  
  group_list = factor(d28[[i]])
  design <- model.matrix(~group_list)
  colnames(design) <- levels(group)
  rownames(design) <- colnames(exprSet)
  
  dge <- DGEList(counts = exprSet,group = group_list)
  ### filter base  use CPM
  #keep <- rowSums(cpm(dge) > 0.5 ) >=2
  table(keep)
  dge.filter <- dge[keep,keep.lib.sizes=FALSE]
  dge.norm <- calcNormFactors(dge.filter,method = "TMM")
  
  ##### limma-trend #####    
  logCPM <- cpm(dge.norm, log=TRUE, prior.count=3)
  fit <- lmFit(logCPM, design)
  fit <- eBayes(fit,trend = T)
  diff_trend <- topTable(fit, coef=ncol(design),n=Inf) ## 可以取出不显著的值
  
  
  diff_trend[which(diff_trend$logFC > F.change & diff_trend$P.Value < 0.05),'sig'] <- 'up'
  diff_trend[which(diff_trend$logFC < -F.change & diff_trend$P.Value < 0.05),'sig'] <- 'down'
  diff_trend[which(abs(diff_trend$logFC) <= F.change| diff_trend$P.Value >= 0.05),'sig'] <- 'none'
  DEG_trend<- subset(diff_trend, sig %in% c('up', 'down'))
  
  limma_trend <- row.names(DEG_trend)
  DEG_trend_UP <- DEG_trend[with(DEG_trend, logFC > F.change), ]
  DEG_trend_DOWN <- DEG_trend[with(DEG_trend, logFC < -F.change), ]
  
  trend_exp<- apply(logCPM,2,expo)
  
  
  
  ## 保存DE分析结果
  write.csv(trend_exp,file=paste0(loc1,G7cp[i],"_mRNA_trend_exp.csv"),row.names=T)
  write.csv(diff_trend,file=paste0(loc1,G7cp[i],"_mRNA_trend.csv"),row.names=T)
  #deg_all_none = topTable(fit2,adjust.method="none",coef=1,adj.P.Val=0.05,
  #lfc=log(2,2),number=5000,sort.by = 'logFC')
  write.csv(DEG_trend,file=paste0(loc1,G7cp[i],"_mRNA_Diff_trend_ALL.csv"),row.names = T)
  write.csv(DEG_trend_UP,file=paste0(loc1,G7cp[i],"_mRNA_Diff_trend_UP.csv"),row.names = T)
  write.csv(DEG_trend_DOWN,file=paste0(loc1,G7cp[i],"_mRNA_Diff_trend_DOWN.csv"),row.names = T)
  write.table(limma_trend,file=paste0(loc1,G7cp[[i]],"_DEmRNA_trend.txt"),sep = "\t",quote = F,row.names = F,col.names = F)
  
  
  
  
  
  
  ##### limma-voom #####
  v <- voom(dge.norm, design,plot=T)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  output <- topTable(fit, coef=2,n=Inf,adjust.method = "none")
  diff_voom = na.omit(output)
  
  diff_voom[which(diff_voom$logFC > F.change & diff_voom$P.Value < P.Val),'sig'] <- 'up'
  diff_voom[which(diff_voom$logFC < -F.change & diff_voom$P.Value < P.Val),'sig'] <- 'down'
  diff_voom[which(abs(diff_voom$logFC) <= F.change| diff_voom$P.Value >= P.Val),'sig'] <- 'none'
  DEG_voom<- subset(diff_voom, sig %in% c('up', 'down'))
  
  DEG_voom_UP<- DEG_voom[with(DEG_voom, logFC >1), ]
  DEG_voom_DOWN<- DEG_voom[with(DEG_voom, logFC < -1), ]
  limma_voom<-row.names(DEG_voom)
  
  voom_exp <- v$E
  
  
  
  # 写入DE分析结果
  write.csv(voom_exp,file=paste0(loc1,G7cp[i],"_mRNA_voom_exp.csv"),row.names = T)
  write.csv(DEG_voom,file=paste0(loc1,G7cp[i],"_mRNA_voom.csv"),row.names = T)
  write.csv(DEG_voom_UP,file=paste0(loc1,G7cp[i],"_mRNA_Diff_voom_UP.csv"),row.names = T)
  write.csv(DEG_voom_DOWN,file=paste0(loc1,G7cp[i],"_mRNA_Diff_voom_DOWN.csv"),row.names = T)
  write.table(limma_voom,file=paste0(loc1,G7cp[[i]],"_DEmRNA_voom.txt"),sep = "\t",quote = F,row.names = F,col.names = F)
  
  
}


for (i in 2:3) { 
  exprSet=cbind(datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==G28name[4]],datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==G28name[i]]) # 表达数据
  exprSet=exprSet[rowMeans(exprSet)>1,] 
  
  dge <- DGEList(counts = exprSet)
  dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log=TRUE, prior.count=3)
  
  group_list=d28[[i]]
  group_list=factor(group_list)
  design <- model.matrix(~group_list)
  colnames(design) <- levels(group_list)
  rownames(design) <- colnames(exprSet)
  
  fit <- lmFit(logCPM, design)
  fit <- eBayes(fit)
  diff_trend <- topTable(fit, coef=2,n=Inf) ## 可以取出不显著的值
  diff_trend[which(diff_trend$logFC > F.change & diff_trend$P.Value < P.Val),'sig'] <- 'up'
  diff_trend[which(diff_trend$logFC < -F.change & diff_trend$P.Value < P.Val),'sig'] <- 'down'
  diff_trend[which(abs(diff_trend$logFC) <= F.change| diff_trend$P.Value >= P.Val),'sig'] <- 'none'
  DEG_trend<- subset(diff_trend, sig %in% c('up', 'down'))
  
  limma_trend <- row.names(DEG_trend)
  DEG_trend_UP <- DEG_trend[with(DEG_trend, logFC > F.change), ]
  DEG_trend_DOWN <- DEG_trend[with(DEG_trend, logFC < -F.change), ]
  
  trend_exp<- apply(logCPM,2,expo)
  
  ## 写入DE文件
  write.csv(trend_exp,file=paste0("./mRNA+lncRNA/Tseries/pc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_mRNA_trend_exp.csv"), row.names = T)
  write.csv(diff_trend, file =paste0("./mRNA+lncRNA/Tseries/pc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_mRNA_trend.csv"), row.names = T)
  write.csv(DEG_trend, file =paste0("./mRNA+lncRNA/Tseries/pc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_mRNA_Diff_trend_ALL.csv"), row.names = T)
  write.csv(DEG_trend_UP, file=paste0("./mRNA+lncRNA/Tseries/pc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_mRNA_Diff_trend_UP.csv"), row.names = T)
  write.csv(DEG_trend_DOWN, file=paste0("./mRNA+lncRNA/Tseries/pc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_mRNA_Diff_trend_DOWN.csv"), row.names = T)
  write.table(limma_trend,file=paste0("./mRNA+lncRNA/Tseries/pc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_DEmRNA_trend.txt"),sep = "\t",quote = F,row.names = F,col.names = F)
  
  ## ID 转换
  gene<- as.data.frame(limma_trend)
  colnames(gene)="ensembl_transcript_id"
  gene_name<-getBM(attributes=c("ensembl_transcript_id","external_gene_name","ensembl_gene_id"),filters = "ensembl_transcript_id",values = gene, mart = mart)
  
  DEG_trends=DEG_trend
  DEG_trends$ensembl_transcript_id=rownames(DEG_trends)
  
  DEG_trends<-DEG_trends[match(gene_name$ensembl_transcript_id,DEG_trends$ensembl_transcript_id),]
  
  deg_trend<-merge(DEG_trends,gene_name,by="ensembl_transcript_id")
  deglist_trend<-unique(deg_trend$external_gene_name)
  
  ## 保存Venn结果
  write.csv(deg_trend, file =paste0(loc2,G28cp[[i]],"_DEgene_trend_exp.csv"), row.names = T)
  write.table(deglist_trend,file=paste0(loc2,G28cp[[i]],"_DEgene_trend.txt"),sep = "\t",quote = F,row.names = F,col.names = "symbol")
  
}

# bpd VS con 
i=1
exprSet=cbind(datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==G28name[4]],datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==G28name[i]]) # 表达数据
exprSet=exprSet[rowMeans(exprSet)>1,] 

dge <- DGEList(counts = exprSet)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3)

group_list = d28[[i]]
group_list = factor(group_list)
design <- model.matrix(~ group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(exprSet)

fit <- lmFit(logCPM, design)
fit <- eBayes(fit)
diff_trend <- topTable(fit, coef=2,n=Inf) ## 可以取出不显著的值
diff_trend[which(diff_trend$logFC > F.change & diff_trend$P.Value < 0.05),'sig'] <- 'up'
diff_trend[which(diff_trend$logFC < -F.change & diff_trend$P.Value < 0.05),'sig'] <- 'down'
diff_trend[which(abs(diff_trend$logFC) <= F.change| diff_trend$P.Value >= 0.05),'sig'] <- 'none'
DEG_trend<- subset(diff_trend, sig %in% c('up', 'down'))

limma_trend <- row.names(DEG_trend)
DEG_trend_UP <- DEG_trend[with(DEG_trend, logFC > F.change), ]
DEG_trend_DOWN <- DEG_trend[with(DEG_trend, logFC < -(F.change)), ]

trend_exp<- apply(logCPM,2,expo)

# 写入DE文件
write.csv(logCPM,file=paste0("./mRNA+lncRNA/Tseries/pc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_mRNA_trend_exp.csv"), row.names = T)
write.csv(diff_trend, file =paste0("./mRNA+lncRNA/Tseries/pc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_mRNA_trend.csv"), row.names = T)
write.csv(DEG_trend, file =paste0("./mRNA+lncRNA/Tseries/pc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_mRNA_Diff_trend_ALL.csv"), row.names = T)
write.csv(DEG_trend_UP, file=paste0("./mRNA+lncRNA/Tseries/pc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_mRNA_Diff_trend_UP.csv"), row.names = T)
write.csv(DEG_trend_DOWN, file=paste0("./mRNA+lncRNA/Tseries/pc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_mRNA_Diff_trend_DOWN.csv"), row.names = T)
write.table(limma_trend,file=paste0("./mRNA+lncRNA/Tseries/pc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_DEmRNA_trend.txt"),sep = "\t",quote = F,row.names = F,col.names = F)


# ID转换
gene<- as.data.frame(limma_trend)
colnames(gene)="ensembl_transcript_id"
gene_name<-getBM(attributes=c("ensembl_transcript_id","external_gene_name","ensembl_gene_id"),filters = "ensembl_transcript_id",values = gene, mart = mart)

DEG_trends=DEG_trend
DEG_trends$ensembl_transcript_id=rownames(DEG_trends)

DEG_trends<-DEG_trends[match(gene_name$ensembl_transcript_id,DEG_trends$ensembl_transcript_id),]

deg_trend<-merge(DEG_trends,gene_name,by="ensembl_transcript_id")
deglist_trend<-unique(deg_trend$external_gene_name)

## 保存Venn结果
write.csv(deg_trend, file =paste0(loc2,G28cp[[i]],"_DEgene_trend_exp.csv"), row.names = T)
write.table(deglist_trend,file=paste0(loc2,G28cp[[i]],"_DEgene_trend.txt"),sep = "\t",quote = F,row.names = F,col.names = "symbol")

##### limma-voom #####
##### pValue #####
for (i in 2:3){
  exprSet=cbind(datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==Gname[4]],datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==Gname[i]]) # 表达数据
  exprSet=exprSet[rowMeans(exprSet)>1,] 
  
  dge <- DGEList(counts = exprSet)
  dge <- calcNormFactors(dge)
  group_list=factor(d28[[i]])
  design <- model.matrix(~group_list)
  colnames(design) <- levels(group_list)
  rownames(design) <- colnames(exprSet)
  
  v <- voom(dge, design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  output <- topTable(fit, coef=2,n=Inf)
  diff_voom = na.omit(output)
  DEG_voom <- diff_voom[with(diff_voom, ((logFC > F.change| logFC < (-F.change)) & P.Value < P.Val)), ]
  limma_voom <- row.names(DEG_voom)
  DEG_voom_UP<- DEG_voom[with(DEG_voom, ((logFC > F.change) & P.Value < P.Val )), ]
  DEG_voom_DOWN<- DEG_voom[with(DEG_voom, ((logFC < - F.change) & P.Value < P.Val)), ]   
  voom_exp <- v$E
  voom_exp=apply(voom_exp,2,expo)
  
  write.csv(voom_exp, file =paste0("./mRNA+lncRNA/Tseries/pc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_mRNA_voom_exp.csv"), row.names = T)
  write.csv(diff_voom, file =paste0("./mRNA+lncRNA/Tseries/pc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_mRNA_voom.csv"), row.names = T)
  write.csv(DEG_voom, file =paste0("./mRNA+lncRNA/Tseries/pc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_mRNA_Diff_voom_ALL.csv"), row.names = T)
  write.csv(DEG_voom_UP, file=paste0("./mRNA+lncRNA/Tseries/pc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_mRNA_Diff_voom_UP.csv"), row.names = T)
  write.csv(DEG_voom_DOWN, file=paste0("./mRNA+lncRNA/Tseries/pc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_mRNA_Diff_voom_DOWN.csv"), row.names = T)
  write.table(limma_voom,file=paste0("./mRNA+lncRNA/Tseries/pc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_DEmRNA_voom.txt"),sep = "\t",quote = F,row.names = F,col.names = F)
  
  #ID转换
  gene<- as.data.frame(limma_voom)
  colnames(gene)="ensembl_transcript_id"
  gene_name<-getBM(attributes=c("ensembl_transcript_id","external_gene_name","ensembl_gene_id"),filters = "ensembl_transcript_id",values = gene, mart = mart)
  DEG_vooms=DEG_voom
  DEG_vooms$ensembl_transcript_id=rownames(DEG_vooms)
  
  DEG_vooms<-DEG_vooms[match(gene_name$ensembl_transcript_id,DEG_vooms$ensembl_transcript_id),]
  
  deg_voom<-merge(DEG_vooms,gene_name,by="ensembl_transcript_id")
  deglist_voom<-unique(deg_voom$external_gene_name)
  
  # 写入Venn分析结果
  write.csv(deg_voom, file =paste0(loc2,G28cp[[i]],"_DEgene_voom_exp.csv"), row.names = T)
  write.table(deglist_voom,file=paste0(loc2,G28cp[[i]],"_DEgene_voom.txt"),sep = "\t",quote = F,row.names = F,col.names = "symbol")
}   

# bpd VS control 
i=1
exprSet=cbind(datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==Gname[i]],datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==Gname[4]]) # 表达数据
exprSet=exprSet[rowMeans(exprSet)>1,] 

dge <- DGEList(counts = exprSet)
dge <- calcNormFactors(dge)
group_list=factor(d28[[i]])
design <- model.matrix(~ group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(exprSet)

v <- voom(dge, design)
fit <- lmFit(v, design)
fit <- eBayes(fit)
output <- topTable(fit, coef=2,n=Inf)
diff_voom = na.omit(output)
DEG_voom <- diff_voom[with(diff_voom, ((logFC > F.change| logFC< (-F.change)) & P.Value < P.Val )), ]
limma_voom <- row.names(DEG_voom)
DEG_voom_UP<- DEG_voom[with(DEG_voom, ((logFC > F.change) & P.Value < P.Val )), ]

voom_exp <- v$E
voom_exp=apply(voom_exp,2,expo)

write.csv(voom_exp,file =paste0("./mRNA+lncRNA/Tseries/pc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_mRNA_voom_exp.csv"), row.names = T)
write.csv(diff_voom, file =paste0("./mRNA+lncRNA/Tseries/pc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_mRNA_voom.csv"), row.names = T)
write.csv(DEG_voom, file =paste0("./mRNA+lncRNA/Tseries/pc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_mRNA_Diff_voom_ALL.csv"), row.names = T)
write.csv(DEG_voom_UP, file=paste0("./mRNA+lncRNA/Tseries/pc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_mRNA_Diff_voom_UP.csv"), row.names = T)
write.csv(DEG_voom_DOWN, file=paste0("./mRNA+lncRNA/Tseries/pc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_mRNA_Diff_voom_DOWN.csv"), row.names = T)
write.table(limma_voom,file=paste0("./mRNA+lncRNA/Tseries/pc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_DEmRNA_voom.txt"),sep = "\t",quote = F,row.names = F,col.names = F)

#ID转换
gene<- as.data.frame(limma_voom)
colnames(gene)="ensembl_transcript_id"
gene_name<-getBM(attributes=c("ensembl_transcript_id","external_gene_name","ensembl_gene_id"),filters = "ensembl_transcript_id",values = gene, mart = mart)
DEG_vooms=DEG_voom
DEG_vooms$ensembl_transcript_id=rownames(DEG_vooms)

DEG_vooms<-DEG_vooms[match(gene_name$ensembl_transcript_id,DEG_vooms$ensembl_transcript_id),]

deg_voom<-merge(DEG_vooms,gene_name,by="ensembl_transcript_id")
deglist_voom<-unique(deg_voom$external_gene_name)

# 写入Venn分析结果
write.csv(deg_voom, file =paste0(loc2,G28cp[[i]],"_DEgene_voom_exp.csv"), row.names = T)
write.table(deglist_voom,file=paste0(loc2,G28cp[[i]],"_DEgene_voom.txt"),sep = "\t",quote = F,row.names = F,col.names = "symbol")


