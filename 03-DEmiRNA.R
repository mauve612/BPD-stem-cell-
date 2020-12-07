rm(list=ls())
dir<- "./count/米RNA/"
file=list.files(path = dir)

library(stringr)
# 测试数据的读取过程  
test1 <- file[1]

tmp <- as.data.frame(readxl::read_xlsx(path = paste0(dir,test1),sheet = 1)
)
temp = tmp[,-1]
temp=temp[,c(1:6)]
rownames(temp)=tmp[,1]
head(tmp)

expr <- lapply(file,
               function(x){
                 tmp <- as.data.frame(readxl::read_xlsx(path = paste0(dir,x),sheet = 1)
                 )
                 temp = tmp[,-1]
                 temp=temp[,c(1:6)]
                 rownames(temp)=tmp[,1]
                 head(temp)
                 return(temp)
               })

#整理数据
df <- do.call(cbind, expr)
#是否去除NA值需要根据自己的需求决定
df <- na.omit(df)
colnames(df)<-str_sub(colnames(df),start = 7,end = -1L)
df=df[,!duplicated(colnames(df))]
df2=df[order(colnames(df2),decreasing=T)]
write.table(df2,file = "./count/miRNA_count.txt",sep = "\t",quote = F,row.names = T,col.names = T)
save(df2,file="miRNA_count.Rdata")


