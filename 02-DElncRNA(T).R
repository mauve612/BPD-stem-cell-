library("biomaRt")
mart <- useMart("ensembl","mmusculus_gene_ensembl")##人类选择hsapiens_gene_ensembl


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
write.csv(deg_trend, file =paste0(loc2,G7cp[[i]],"_DEgene_trend_exp.csv"), row.names = T)
write.table(deglist_trend,file=paste0(loc2,G7cp[[i]],"_DEgene_trend.txt"),sep = "\t",quote = F,row.names = F,col.names = "symbol")

##ID转换提取差异基因
gene<- as.data.frame(DEG_DESeq2)
colnames(gene)="ensembl_transcript_id"
gene_name<-getBM(attributes=c("ensembl_transcript_id","external_gene_name","ensembl_gene_id"),filters = "ensembl_transcript_id",values = gene, mart = mart)
diffDatas=diffData
colnames(diffDatas)[1]="ensembl_transcript_id"

diffDatas<-diffDatas[match(gene_name$ensembl_transcript_id,diffDatas$ensembl_transcript_id),]

deg_DESeq2<-merge(diffDatas,gene_name,by="ensembl_transcript_id")
deglist_DESeq2<-unique(deg_DESeq2$external_gene_name)

#写入Venn分析文件
write.csv(deg_DESeq2, file =paste0(loc2,G7cp[[i]],"_DEgene_DESeq2_exp.csv"), row.names = T)
write.table(deglist_DESeq2,file=paste0(loc2,G7cp[[i]],"_DEgene_DESeq2.txt"),sep = "\t",quote = F,row.names = F,col.names = "symbol")


##ID转换提取差异基因
gene<- as.data.frame(DEG_DESeq2)
colnames(gene)="ensembl_transcript_id"
gene_name<-getBM(attributes=c("ensembl_transcript_id","external_gene_name","ensembl_gene_id"),filters = "ensembl_transcript_id",values = gene, mart = mart)
diffDatas=diffData
colnames(diffDatas)[1]="ensembl_transcript_id"

diffDatas<-diffDatas[match(gene_name$ensembl_transcript_id,diffDatas$ensembl_transcript_id),]

deg_DESeq2<-merge(diffDatas,gene_name,by="ensembl_transcript_id")
deglist_DESeq2<-unique(deg_DESeq2$external_gene_name)

##写入Venn文件
write.csv(deg_DESeq2, file =paste0(loc2,"39_DEgene_DESeq2.csv"), row.names = T)
write.table(deglist_DESeq2,file=paste0(loc2,"39_DEgene_DESeq2.txt"),sep = "\t",quote = F,row.names = F,col.names = "symbol")

## ID转化
gene<- as.data.frame(DEG_edgeR)
colnames(gene)="ensembl_transcript_id"
gene_name<-getBM(attributes=c("ensembl_transcript_id","external_gene_name","ensembl_gene_id"),filters = "ensembl_transcript_id",values = gene, mart = mart)
Gene_diff_edgeRs=Gene_diff_edgeR
Gene_diff_edgeRs$ensembl_transcript_id=row.names(Gene_diff_edgeR)

Gene_diff_edgeRs<-Gene_diff_edgeRs[match(gene_name$ensembl_transcript_id,Gene_diff_edgeRs$ensembl_transcript_id),]

deg_edgeR<-merge(Gene_diff_edgeRs,gene_name,by="ensembl_transcript_id")
deglist_edgeR<-unique(deg_edgeR$external_gene_name)

## 写入venn
write.csv(deg_edgeR, file =paste0(loc2,"39_DEgene_edgeR_exp.csv"), row.names = T)
write.table(deglist_edgeR,file=paste0(loc2,"39_DEgene_edgeR.txt"),sep = "\t",quote = F,row.names = F,col.names = "symbol")





rm(list=ls())
library("stringr")
library("DESeq2")
library("edgeR")
library("limma")

#### 导入数据  ####
load("./mRNA+lncRNA/lnc_7days.Rdata")
load("./mRNA+lncRNA/lnc_28days.Rdata")

#设定筛选的阈值 
F.change = 1.00
p.adj = 0.05
P.Val = 0.05

#### 7days  ####
#### 批量运算参数设置  ####
#数据提取 
c1=c("control","control","control","bpd","bpd","bpd")

# 文件命名
Gname<-c("BPD","BPDM","BPDN","Con")

G7cp<-list(compare1="B vs A",
           compare2="C vs A",
           compare3="D vs A"
)

# 文件路径
dir.create(path="./mRNA+lncRNA/Tseries/lnc/",recursive = T)
dir.create(path="./mRNA+lncRNA/Tseries/lnc/7days",recursive = T)

#### 3 bpd vs 3 control ####
#####  DESeq2 #####
datcount=Tlnc7d[,-1]
rownames(datcount)=Tlnc7d$`#id`
datcount=datcount[,c(1:12)] 
datcount[1:8,1:4]
for (i in 1:3){
  datacount <- cbind(datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==Gname[4]],datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==Gname[i]])
  ##构建dds对象
  datacount <- round(as.matrix(datacount))
  groups<- c1
  Data <- data.frame(condition = factor(groups),
                     row.names = colnames(datacount))
  
  dds <-DESeqDataSetFromMatrix(countData = datacount, 
                               colData = Data, design= ~condition)
  
  dds$condition <- relevel(dds$condition, ref = "control")
  #levels(dds$condition)
   
  ## 去掉所有条件都没有read的基因
  dds <- dds[rowSums(counts(dds))>1,]
  ##使用DESeq函数预估离散度
  dds <-DESeq(dds)
  res <- results(dds)
  ##设定阈值筛选差异基因
  res <- res[order(res$padj), ]
  resData <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  resData <-resData[order(resData$padj < p.adj),]
  resData[which(res$log2FoldChange > F.change & resData$pvalue < P.Val),'sig'] <- 'up'
  resData[which(resData$log2FoldChange < -F.change & resData$pvalue < P.Val),'sig'] <- 'down'
  resData[which(abs(resData$log2FoldChange) <= F.change| resData$pvalue >= P.Val),'sig'] <- 'none'
  diffData <- subset(resData, sig %in% c('up', 'down'))
  diffData = na.omit(diffData)
  
  diffDataUP = diffData[(diffData$pvalue < P.Val & (diffData$log2FoldChange > F.change)),]
  diffDataDOWN = diffData[(diffData$pvalue < P.Val & (diffData$log2FoldChange < (-F.change))),]
 
  diffData=na.omit(diffData)
  diffDataUP = diffData[(diffData$pvalue < P.Val & (diffData$log2FoldChange > F.change)),]
  diffDataDOWN = diffData[(diffData$pvalue < P.Val & (diffData$log2FoldChange < (-F.change))),]
  
  DEG_DESeq2 <- row.names(diffData)
  
  exp_res_Data=resData[,c(8:length(colnames(resData)))]
  row.names(exp_res_Data)=resData$Row.names
  
  # 写入完整文件 
  write.csv(exp_res_Data, file =paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_",G7cp[[i]],"_lnc_DESeq2_exp.csv"), row.names = T)
  write.csv(resData, file =paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_",G7cp[[i]],"_lnc_DESeq2.csv"), row.names = T)
  write.csv(diffData, file =paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_",G7cp[[i]],"_lnc_Diff_DESeq2_ALL.csv"), row.names = T)
  write.csv(diffDataUP, file=paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_",G7cp[[i]],"_lnc_Diff_DESeq2_UP.csv"), row.names = T)
  write.csv(diffDataDOWN, file=paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_",G7cp[[i]],"_lnc_Diff_DESeq2_DOWN.csv"), row.names = T)
  #写入DEG 列表
  write.table(DEG_DESeq2,file=paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_",G7cp[[i]],"_DElnc_DESeq2.txt"),sep = "\t",quote = F,row.names = F,col.names = F)
  
}

##### edgeR #####
#读取数据
datcount=Tlnc7d[,-1]
rownames(datcount)=Tlnc7d$`#id`
datcount=datcount[,c(1:12)] 
exprSet <- datcount
for (i in 1:3){
  exprSet <- cbind(datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==Gname[i]],datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==Gname[4]])
  ##分为两组
  groups <- factor(c(rep("control", 3),rep("bpd", 3)), levels = c("bpd", "control"))
  exprSet <- exprSet[rowMeans(exprSet) > 1,]
  # exprSet <- exprSet[rowSums(cpm(exprSet) > 1) >= 2, ]
  ##构建DGEList对象
  exprSet <- DGEList(counts = exprSet, group = groups)
  exprSet <- calcNormFactors(exprSet)
  ##估计离散度
  exprSet <- estimateCommonDisp(exprSet)
  exprSet <- estimateTagwiseDisp(exprSet)
  ##找差异基因，并按阈值进行筛选
  list <- exactTest(exprSet,pair = c("control","bpd"))
  Gene_diff <- topTags(list, n=nrow(exprSet))
  ##筛选差异表达基因
  Gene_edgeR <- Gene_diff$table
 
  #在总表中对应的标记出来基因具体的上下调
  Gene_edgeR[which(Gene_edgeR$logFC > F.change & Gene_edgeR$FDR < p.adj),'sig'] <- 'up'
  Gene_edgeR[which(Gene_edgeR$logFC < -F.change & Gene_edgeR$FDR < p.adj),'sig'] <- 'down'
  Gene_edgeR[which(abs(Gene_edgeR$logFC) <= F.change| Gene_edgeR$FDR >= p.adj),'sig'] <- 'none'
  Gene_diff_edgeR <- subset(Gene_edgeR, sig %in% c('up', 'down'))
  
  DEG_edgeR <-row.names(Gene_diff_edgeR)
  edgeR_UP <- Gene_diff_edgeR[Gene_diff_edgeR$logFC > F.change & Gene_diff_edgeR$FDR < p.adj,]
  edgeR_DOWN <- Gene_diff_edgeR[Gene_diff_edgeR$logFC < -(F.change) & Gene_diff_edgeR$FDR < p.adj,]
  exp_edgeR <- exprSet$pseudo.counts
  
  # 写入差异基因数据
  write.csv(exp_edgeR, file = paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_",G7cp[[i]],"_lnc_edgeR_exp.csv"))
  write.csv(Gene_edgeR, file = paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_",G7cp[[i]],"_lnc_edgeR.csv"))
  write.csv(Gene_diff_edgeR, file = paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_",G7cp[[i]],"_lnc_Diff_edgeR_ALL.csv"))
  write.csv(edgeR_UP, file = paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_",G7cp[[i]],"_lnc_Diff_edgeR_UP.csv"))
  write.csv(edgeR_DOWN, file = paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_",G7cp[[i]],"_lnc_Diff_edgeR_DOWN.csv"))
  
  # 写入DEG 列表
  write.table(DEG_edgeR,file=paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_",G7cp[[i]],"_DElnc_edgeR.txt"),sep = "\t",quote = F,row.names = F,col.names = F)
}

##### limma #####
##### limma-trend #####    
for (i in 1:3){
  exprSet=cbind(datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==Gname[4]],datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==Gname[i]]) # 表达数据
  exprSet = as.matrix(exprSet)  #这一步的前处理非常关键，决定是否为合适的表达矩阵
  exprSet = exprSet[rowMeans(exprSet) > 1,] 
  
  dge <- DGEList(counts = exprSet)
  dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log=TRUE, prior.count=3)
  
 
  group_list = factor(c(rep("control",3),rep("disease",3))) ## factor 自带参数的bug?
  design <- model.matrix(~group_list)
  colnames(design) <- levels(group_list)
  rownames(design) <- colnames(counts)
  
  fit <- lmFit(logCPM, design)
  fit <- eBayes(fit,trend = T)
  diff_trend <- topTable(fit, coef=2,n=Inf) ## 可以取出不显著的值
  
  diff_trend[which(diff_trend$logFC > F.change & diff_trend$P.Value < P.Val),'sig'] <- 'up'
  diff_trend[which(diff_trend$logFC < -F.change & diff_trend$P.Value < P.Val),'sig'] <- 'down'
  diff_trend[which(abs(diff_trend$logFC) <= F.change| diff_trend$P.Value >= P.Val),'sig'] <- 'none'
  DEG_trend <- subset(diff_trend, sig %in% c('up', 'down'))
  
  limma_trend <- row.names(DEG_trend)
  DEG_trend_UP <- DEG_trend[with(DEG_trend, (logFC > F.change & P.Value < P.Val )), ]
  DEG_trend_DOWN <- DEG_trend[with(DEG_trend, (logFC < -F.change & P.Value < P.Val )), ]
  
  # save the diff and DEG of the trend analysis pathway
  write.csv(logCPM,file=paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_",G7cp[i],"lnc_trend_exp.csv"),row.names=T)
  write.csv(diff_trend,file=paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_",G7cp[i],"_lnc_trend.csv"),row.names=T)
  write.csv(DEG_trend,file=paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_",G7cp[i],"_lnc_Diff_trend_ALL.csv"),row.names = T)
  write.csv(DEG_trend_UP,file=paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_",G7cp[i],"_lnc_Diff_trend_UP.csv"),row.names = T)
  write.csv(DEG_trend_DOWN,file=paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_",G7cp[i],"_lnc_Diff_trend_DOWN.csv"),row.names = T)
  write.table(limma_trend,file=paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_",G7cp[[i]],"_DElnc_trend.txt"),sep = "\t",quote = F,row.names = F,col.names = F)
}
# count 转化为 log2的值进行标准化，标准化方法为TMM
#使用edgeR中的calcNormFactor进行计算即可,edgeR+limma联合分析
# !! 该方法适用于组见样本测序深度差异不大的情况

##### limma-voom #####
datcount=Tlnc7d[,-1]
rownames(datcount)=Tlnc7d$`#id`
datcount=datcount[,c(1:12)] 
for (i in 1:3){
  exprSet=cbind(datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==Gname[4]],datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==Gname[i]]) # 表达数据
  exprSet=exprSet[rowMeans(exprSet)>1,] 
  
  dge <- DGEList(counts = exprSet)
  dge <- calcNormFactors(dge)
  group_list=factor(c(rep("control",3),rep("disease",3)))
  design <- model.matrix(~group_list)
  colnames(design) <- levels(group_list)
  rownames(design) <- colnames(exprSet)
  
  v <- voom(dge, design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  output <- topTable(fit, coef=2,n=Inf)
  diff_voom = na.omit(output)
  
  DEG_voom <- diff_voom[with(diff_voom, ((logFC > F.change| logFC < (-F.change)) & adj.P.Val < p.adj )), ]
  limma_voom <- row.names(DEG_voom)
  DEG_voom_UP<- DEG_voom[with(DEG_voom, ((logFC > F.change) & adj.P.Val < p.adj)), ]
  DEG_voom_DOWN<- DEG_voom[with(DEG_voom, ((logFC < -F.change) & adj.P.Val < p.adj)), ]
  voom_exp=v$E
 
   # 写入差异基因数据
  write.csv(voom_exp,file=paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_",G7cp[i],"_lnc","voom_exp.csv"),row.names = T)
 write.csv(diff_voom,file=paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_",G7cp[i],"_lnc","_limma_voom.csv"),row.names = T)
 write.csv(DEG_voom,file=paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_",G7cp[i],"_lnc","_voom_ALL.csv"),row.names = T)
 write.csv(DEG_voom_UP,file=paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_",G7cp[i],"_lnc","_voom_UP.csv"),row.names = T)
  write.csv(DEG_voom_DOWN,file=paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_",G7cp[i],"_lnc","_voom_DOWN.csv"),row.names = T)
  
  # 写入DEG 列表
  write.table(limma_voom,file=paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_",G7cp[[i]],"_DElnc_limma(voom).txt"),sep = "\t",quote = F,row.names = F,col.names = F)
  
}



##### 9 bpd VS 3 con 进行的分析 #####
###### DESeq2 ######
datacount=Tlnc7d[,-1]
rownames(datacount)=Tlnc7d$`#id`
datacount=datacount[,c(1:12)]
##构建dds对象
datacount <- round(as.matrix(datacount))
##分为control和bpd两组
groups<- c(rep("bpd",9),rep("control",3))
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
Gene_diff <- subset(res, padj < p.adj & (log2FoldChange > F.change | log2FoldChange < -F.change))
resData <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
resData <-resData[order(resData$padj),]
diffData<-as.data.frame(Gene_diff)
diffData = na.omit(diffData)
diffDataUP<-diffData[with(diffData, log2FoldChange > F.change ),]
diffDataDown<-diffData[with(diffData, log2FoldChange < -F.change ),]

Gene_diff_DESeq2 <- row.names(diffData)

exp_res_Data=resData[,c(8:length(colnames(resData)))]
row.names(exp_res_Data)=resData$Row.names

# 写入完整文件 
write.csv(exp_res_Data, file =paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_39_lnc_DESeq2_exp.csv"), row.names = T)
write.csv(resData, file ="./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_39_lnc_DESeq2.csv", row.names = T)
write.csv(diffData, file ="./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_39_Diff_DESeq2_ALL.csv", row.names = T)
write.csv(diffDataUP, file ="./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_39_Diff_DESeq2_UP.csv", row.names = T)
write.csv(diffDataDown, file ="./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_39_Diff_DESeq2_DOWN.csv", row.names = T)
# 写入DEG 列表
write.table(Gene_diff_DESeq2,file="./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_39_DElnc_DESeq2.txt",sep = "\t",quote = F,row.names = F,col.names = F)

###### edgeR ######
#读取数据
exprSet <- datacount[,c(10:12,1:9)]
##分为两组
groups <- c(rep("control",3),rep("disease",9))
exprSet <- exprSet[rowMeans(exprSet) > 1,]


design = factor(groups)
design <- model.matrix(~groups)

##构建DGEList对象
exprSet <- DGEList(counts = exprSet, group = groups)
exprSet <- calcNormFactors(exprSet)
##估计离散度
exprSet <- estimateCommonDisp(exprSet)
exprSet <- estimateTagwiseDisp(exprSet)
##找差异基因，并按阈值进行筛选
list <- exactTest(exprSet,pair = c("control","disease"))
Gene_diff <- topTags(list, n=nrow(exprSet))
##筛选差异表达基因
Gene_edgeR <- Gene_diff$table
#在总表中对应的标记出来基因具体的上下调
Gene_edgeR[which(Gene_edgeR$logFC > F.change & Gene_edgeR$FDR < p.adj),'sig'] <- 'up'
Gene_edgeR[which(Gene_edgeR$logFC < -F.change & Gene_edgeR$FDR < p.adj),'sig'] <- 'down'
Gene_edgeR[which(abs(Gene_edgeR$logFC) <= F.change| Gene_edgeR$FDR >= p.adj),'sig'] <- 'none'
Gene_diff_edgeR <- subset(Gene_edgeR, sig %in% c('up', 'down'))
DEG_edgeR <-row.names(Gene_diff_edgeR)
edgeR_UP <-  Gene_diff_edgeR[Gene_diff_edgeR$logFC > F.change, ]
edgeR_DOWN <- Gene_diff_edgeR[Gene_diff_edgeR$logFC < -F.change, ]
exp_edgeR <- exprSet$counts
write.csv(exp_edgeR, file = paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_39_lnc_edgeR_exp.csv"))
write.csv(all_edgeR, file = "./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_39_lnc_edgeR.csv")
write.csv(Gene_diff_edgeR, file = "./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_39_Diff_edgeR.csv")
write.csv(edgeR_UP, file = "./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_39_Diff_edgeR_UP.csv")
write.csv(edgeR_DOWN, file = "./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_39_Diff_edgeR_DOWN.csv")

# 写入DEG 列表
write.table(DEG_edgeR,file="./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_39_DElnc_edgeR.txt",sep = "\t",quote = F,row.names = F,col.names = F)


###### limma ######
##### limma-trend #####    
exprSet <- datacount[,c(10:12,1:9)]
exprSet=exprSet[rowMeans(exprSet)>1,] 

dge <- DGEList(counts = exprSet)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3)

group_list=c(rep("control",3),rep("disease",9))
group_list=factor(group_list)
design <- model.matrix(~group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(count)

fit <- lmFit(logCPM, design)
fit <- eBayes(fit)
diff_trend <- topTable(fit, coef=2,n=Inf) ## 可以取出不显著的值
DEG_trend <- diff_trend[with(diff_trend, ((logFC > F.change| logFC< (-F.change)) & adj.P.Val < p.adj)), ] #可根据实际情况掉配阈值的参数类型，然后进行后续的多分析方法的交集分析
limma_trend <- row.names(DEG_trend)
DEG_trend_UP <- DEG_trend[with(DEG_trend, (logFC > F.change & adj.P.Val < p.adj )), ]
DEG_trend_DOWN <- DEG_trend[with(DEG_trend, (logFC < -F.change & adj.P.Val < p.adj)), ]

# save the diff and DEG of the trend analysis pathway
write.csv(logCPM,file=paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_39_lnc_trend_exp.csv"),row.names=T)
write.csv(diff_trend,file=paste0("./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_39_lnc_trend_ALL.csv"),row.names=T)
write.csv(DEG_trend, file = "./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_39_Diff_trend_ALL.csv")
write.csv(DEG_trend_UP, file = "./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_39_Diff_trend_UP.csv")
write.csv(DEG_trend_DOWN, file = "./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_39_Diff_trend_DOWN.csv")

# 写入DEG 列表
write.table(limma_trend,file="./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_39_DElnc_trend.txt",sep = "\t",quote = F,row.names = F,col.names = F)

##### limma-voom #####
dge <- DGEList(counts = exprSet)
dge <- calcNormFactors(dge)
group_list=factor(c(rep("control",3),rep("disease",9)))
design <- model.matrix(~group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(exprSet)

v <- voom(dge, design)
fit <- lmFit(v, design)
fit <- eBayes(fit)
output <- topTable(fit, coef=2,n=Inf)
diff_voom = na.omit(output)
DEG_voom <- diff_voom[with(diff_voom, ((logFC > F.change| logFC < (-F.change)) & adj.P.Val < p.adj )), ]
limma_voom <- row.names(DEG_voom)
DEG_voom_UP<- DEG_voom[with(DEG_voom, ((logFC > F.change) & adj.P.Val < p.adj )), ]
DEG_voom_DOWN<- DEG_voom[with(DEG_voom, ((logFC < -F.change) & adj.P.Val < p.adj )), ]

voom_exp=v$E

# 写入差异基因数据
write.csv(voom_exp, file = "./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_39_mRNA_voom_exp.csv")
write.csv(diff_voom, file = "./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_39_mRNA_voom.csv")
write.csv(DEG_voom, file = "./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_39_mRNA_Diff_voom_ALL.csv")
write.csv(DEG_voom_UP, file = "./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_39_mRNA_Diff_voom_UP.csv")
write.csv(DEG_voom_DOWN, file = "./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_39_mRNA_Diff_voom_DOWN.csv")

# 写入DEG 列表
write.table(limma_voom,file="./mRNA+lncRNA/Tseries/lnc/7days/bpd_vs_ctrl_7days_39_DEmRNA_voom.txt",sep = "\t",quote = F,row.names = F,col.names = F)


#### 28days ####
# 设置批量保存文件的参数 
#  设置文件保存的位置
name28 =c("control","control","control","msc","msc","msc","mnc","mnc","mnc","bpd","bpd","bpd")

d28 =list(
  name1 = c("control","control","control","bpd","bpd","bpd"),
  name2= c("bpd","bpd","bpd","msc","msc","msc"),
  name3= c("bpd","bpd","bpd","mnc","mnc","mnc")
)

G28cp<-list(compare1="F vs E",
            compare2="G vs F",
            compare3="H vs F"
)
G28name<-c("Con","BPDM","BPDN","BPD")

# 路径设置
dir.create(path="./mRNA+lncRNA/Tseries/lnc/28days",recursive = T)

#因为在实际筛选的过程中发现由于R2较小，得到的DEG总量明显小于7days 分组，因此采用p.adj/FDR和pValue两种矫正方法,可根据实际情况酌情改变阈值的设定
dir.create(path = "./mRNA+lncRNA/Gseries/lnc/28days/pValue")


##### DESeq2 #####
library(DESeq2)
datcount=Tlnc28d[,-1]
rownames(datcount)=Tlnc28d$`#id`
datcount=datcount[,c(1:12)]
datcount=datcount[,c(10:12,4:9,1:3)]
head(datcount)


for (i in 2:3){
  datacount <- cbind(datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==G28name[4]],datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==G28name[i]]) #要比较的组进行前置
  ##构建dds对象
  datacount <- round(as.matrix(datacount))

  groups<-d28[[i]]
  Data <- data.frame(condition = as.factor(groups))
  rownames(Data) <- colnames(datacount)
  dds <-DESeqDataSetFromMatrix(countData = datacount, 
                               colData = Data, design= ~condition)
  dds$condition <- relevel(dds$condition, ref = "bpd")
  
  ## 去掉所有条件都没有read的基因
  dds <- dds[rowSums(counts(dds))>1,]
  ##使用DESeq函数预估离散度
  dds <-DESeq(dds)
  res <- results(dds)
  ##设定阈值筛选差异基因
  
  res <- res[order(res$pvalue), ]
  Gene_diff <- subset(res, pvalue < P.Val & (log2FoldChange > F.change | log2FoldChange < -F.change))
   #直接保存基因条目方便后续做VENN图进行分析计算
  resData <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  resData <-resData[order(resData$pvalue),]
  diffData<-as.data.frame(Gene_diff)
  diffData<-na.omit(diffData)
  diffDataUP = diffData[(diffData$pvalue < P.Val & (diffData$log2FoldChange > F.change)),]
  diffDataDOWN = diffData[(diffData$pvalue < P.Val & (diffData$log2FoldChange < (-F.change))),]
  
  Gene_diff_DESeq2 <- row.names(diffData)
  
  exp_res_Data=resData[,c(8:length(colnames(resData)))]
  row.names(exp_res_Data)=resData$Row.names
 
   # 写入文件
  write.csv(exp_res_Data,file =paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_lnc_DESeq2_exp.csv"), row.names = T)
  write.csv(resData, file =paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_lnc_DESeq2.csv"), row.names = T)
  write.csv(diffData, file =paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_lnc_Diff_DESeq2_ALL.csv"), row.names = T)
  write.csv(diffDataUP, file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_lnc_Diff_DESeq2_UP.csv"), row.names = T)
  write.csv(diffDataDOWN, file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_lnc_Diff_DESeq2_DOWN.csv"), row.names = T)
  write.table(Gene_diff_DESeq2,file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_DElnc_DEseq2.txt"),sep = "\t",quote = F,row.names = F,col.names = F)
}

# ------ bpd vs con ------  
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
Gene_diff <- subset(res, pvalue < P.Val & (log2FoldChange > F.change| log2FoldChange < -F.change))

resData <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
resData <-resData[order(resData$pvalue),]
diffData<-as.data.frame(Gene_diff)
diffData=na.omit(diffData)
diffDataUP = diffData[(diffData$pvalue < P.Val & (diffData$log2FoldChange > F.change)),]
diffDataDOWN = diffData[(diffData$pvalue < P.Val & (diffData$log2FoldChange < (-F.change))),]


Gene_diff_DESeq2 <- row.names(diffData) #直接保存基因条目方便后续做VENN图进行分析计

exp_res_Data=resData[,c(8:length(colnames(resData)))]
row.names(exp_res_Data)=resData$Row.names

# 写入文件
write.csv(exp_res_Data, file =paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_lnc_DESeq2_exp.csv"), row.names = T)
write.csv(resData, file =paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_lnc_DESeq2.csv"), row.names = T)
write.csv(diffData, file =paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_lnc_Diff_DESeq2.csv"), row.names = T)
write.csv(diffDataUP, file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_lnc_Diff_DESeq2_UP.csv"), row.names = T)
write.csv(diffDataDOWN, file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_lnc_Diff_DESeq2_DOWN.csv"), row.names = T)
write.table(Gene_diff_DESeq2,file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_DElnc_DEseq2.txt"),sep = "\t",quote = F,row.names = F,col.names = F)



##### edgeR #####
library(edgeR)
#读取数据
exprSet <- datcount
###### pValue  #####
for (i in 2:3){ 
  exprSet <- cbind(datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==G28name[4]],datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==G28name[i]])
  
  ##分为两组
  groups <- factor(d28[[i]])
  exprSet <- exprSet[rowMeans(exprSet)>1,]
  
  design = factor(groups)
  design <- model.matrix(~groups)
  
  ##构建DGEList对象
  exprSet <- DGEList(counts = exprSet, group = groups)
  exprSet <- calcNormFactors(exprSet)
  ##估计离散度
  exprSet <- estimateCommonDisp(exprSet)
  exprSet <- estimateTagwiseDisp(exprSet)
  ##找差异基因，并按阈值进行筛选
  list <- exactTest(exprSet)
  Gene_diff <- topTags(list, n=nrow(exprSet))
  ##筛选差异表达基因
  Gene_edgeR <- Gene_diff$table
  Gene_edgeR= Gene_edgeR[is.na(Gene_edgeR$FDR)==FALSE,]
  
  Gene_diff_edgeR <- subset(Gene_edgeR, PValue < P.Val & (logFC > F.change| logFC < -(F.change)))
  DEG_edgeR <- row.names(Gene_diff_edgeR)
  edgeR_UP <- Gene_diff_edgeR[Gene_diff_edgeR$logFC > F.change & Gene_diff_edgeR$PValue < P.Val,]
  edgeR_DOWN <- Gene_diff_edgeR[Gene_diff_edgeR$logFC < -(F.change) & Gene_diff_edgeR$PValue < P.Val,]
  exp_edgeR <- exprSet$counts
  
  write.csv(exp_edgeR, file = paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_lnc_edgeR_exp.csv"), row.names = T)
  write.csv(Gene_edgeR, file =paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_lnc_edgeR.csv"), row.names = T)
  write.csv(Gene_diff_edgeR, file =paste0("./mRNA+lncRNA/Gseries/lnc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_lnc_Diff_edgeR.csv"), row.names = T)
  write.csv(edgeR_UP, file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_lnc_Diff_edgeR_UP.csv"), row.names = T)
  write.csv(edgeR_DOWN, file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_lnc_Diff_edgeR_DOWN.csv"), row.names = T)
  write.table(DEG_edgeR,file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_DElnc_edgeR.txt"),sep = "\t",quote = F,row.names = F,col.names = F)
}

# ----- bpd vs control 
i = 1
exprSet <- cbind(datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==G28name[i]],datcount[,str_sub(colnames(datcount),start = 3,end = -2L)==G28name[4]])
##分为两组
groups <- factor(d28[[i]])
exprSet <- exprSet[rowMeans(exprSet)>1,]
##构建DGEList对象
exprSet <- DGEList(counts = exprSet, group = groups)
exprSet <- calcNormFactors(exprSet)
##估计离散度
exprSet <- estimateCommonDisp(exprSet)
exprSet <- estimateTagwiseDisp(exprSet)
##找差异基因，并按阈值进行筛选
list <- exactTest(exprSet,pair = c("control","bpd"))
Gene_diff <- topTags(list, n=nrow(exprSet))
##筛选差异表达基因
Gene_edgeR <- Gene_diff$table
Gene_edgeR= Gene_edgeR[is.na(Gene_edgeR$FDR)==FALSE,]

Gene_diff_edgeR <- subset(Gene_edgeR, PValue < P.Val & (logFC > F.change| logFC < -(F.change)))
DEG_edgeR <- row.names(Gene_diff_edgeR)
edgeR_UP <- Gene_diff_edgeR[Gene_diff_edgeR$logFC > F.change & Gene_diff_edgeR$PValue < P.Val,]
edgeR_DOWN <- Gene_diff_edgeR[Gene_diff_edgeR$logFC < -(F.change) & Gene_diff_edgeR$PValue < P.Val,]
exp_edgeR <- exprSet$counts

write.csv(exp_edgeR, file = paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_lnc_edgeR_exp.csv"), row.names = T)
write.csv(Gene_edgeR, file =paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_lnc_edgeR.csv"), row.names = T)
write.csv(Gene_diff_edgeR, file =paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_lnc_Diff_edgeR.csv"), row.names = T)
write.csv(edgeR_UP, file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_lnc_Diff_edgeR_UP.csv"), row.names = T)
write.csv(edgeR_DOWN, file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_lnc_Diff_edgeR_DOWN.csv"), row.names = T)
write.table(DEG_edgeR,file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_DElnc_edgeR.txt"),sep = "\t",quote = F,row.names = F,col.names = F)


##### #limma ######
##### limma-trend #####    
###### Pvalue  ######
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
  DEG_trend <- diff_trend[with(diff_trend, ((logFC > F.change| logFC< (-F.change)) & P.Value  < P.Val )), ]
  limma_trend <- row.names(DEG_trend)
  DEG_trend_UP <- DEG_trend[with(DEG_trend, (logFC > F.change & P.Value < P.Val )), ]
  DEG_trend_DOWN <- DEG_trend[with(DEG_trend, (logFC < -F.change & P.Value < P.Val )), ]
  
  write.csv(logCPM,file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_lnc_trend_exp.csv"), row.names = T)
  write.csv(diff_trend, file =paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_lnc_trend.csv"), row.names = T)
  write.csv(DEG_trend, file =paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_lnc_Diff_trend_ALL.csv"), row.names = T)
  write.csv(DEG_trend_UP, file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_lnc_Diff_trend_UP.csv"), row.names = T)
  write.csv(DEG_trend_DOWN, file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_lnc_Diff_trend_DOWN.csv"), row.names = T)
  write.table(limma_trend,file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_DElnc_trend.txt"),sep = "\t",quote = F,row.names = F,col.names = F)
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
DEG_trend <- diff_trend[with(diff_trend, ((logFC > F.change| logFC< (- F.change)) & P.Value < P.Val )), ]
limma_trend <- row.names(DEG_trend)
DEG_trend_UP <- DEG_trend[with(DEG_trend, (logFC > F.change & P.Value < P.Val )), ]
DEG_trend_DOWN <- DEG_trend[with(DEG_trend, (logFC < - F.change & P.Value < P.Val )), ]

write.csv(logCPM,file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_lnc_trend_exp.csv"), row.names = T)
write.csv(diff_trend, file =paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_lnc_trend.csv"), row.names = T)
write.csv(DEG_trend, file =paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_lnc_Diff_trend_ALL.csv"), row.names = T)
write.csv(DEG_trend_UP, file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_lnc_Diff_trend_UP.csv"), row.names = T)
write.csv(DEG_trend_DOWN, file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_lnc_Diff_trend_DOWN.csv"), row.names = T)
write.table(limma_trend,file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_DElnc_trend.txt"),sep = "\t",quote = F,row.names = F,col.names = F)


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
  voom_exp<- v$E
  
  write.csv(voom_exp, file =paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_lnc_voom_exp.csv"), row.names = T)
  write.csv(voom_exp_1,file =paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_lnc_voom_exp_1.csv"), row.names = T)
  write.csv(diff_voom, file =paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_lnc_voom.csv"), row.names = T)
  write.csv(DEG_voom, file =paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_lnc_Diff_voom_ALL.csv"), row.names = T)
  write.csv(DEG_voom_UP, file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_lnc_Diff_voom_UP.csv"), row.names = T)
  write.csv(DEG_voom_DOWN, file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_lnc_Diff_voom_DOWN.csv"), row.names = T)
  write.table(limma_voom,file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[i],"_",G28name[4],"_28days_",G28cp[[i]],"_DElnc_voom.txt"),sep = "\t",quote = F,row.names = F,col.names = F)
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
DEG_voom_DOWN<- DEG_voom[with(DEG_voom, ((logFC < - F.change) & P.Value < P.Val)), ]   
voom_exp <- v$E

write.csv(voom_exp,file =paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_lnc_voom_exp.csv"), row.names = T)
write.csv(voom_exp_1,file =paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_lnc_voom_exp_1.csv"), row.names = T)
write.csv(diff_voom, file =paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_lnc_voom.csv"), row.names = T)
write.csv(DEG_voom, file =paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_lnc_Diff_voom_ALL.csv"), row.names = T)
write.csv(DEG_voom_UP, file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_lnc_Diff_voom_UP.csv"), row.names = T)
write.csv(DEG_voom_DOWN, file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_lnc_Diff_voom_DOWN.csv"), row.names = T)
write.table(limma_voom,file=paste0("./mRNA+lncRNA/Tseries/lnc/28days/",G28name[4],"_",G28name[i],"_28days_",G28cp[[i]],"_DElnc_voom.txt"),sep = "\t",quote = F,row.names = F,col.names = F)



