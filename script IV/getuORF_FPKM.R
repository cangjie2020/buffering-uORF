library(dplyr)
library(DT)
library(SummarizedExperiment)
args <- commandArgs(trailingOnly = TRUE)
name<-args[1]                                                                               # 肿瘤名字
getuORF<-read.table(paste(name,"/",name,".getuORF3.txt",sep=""),sep="\t",header=F)      # getuORF.txt
colnames(getuORF)<-c("tran_name","gene_site","trans_site","strand","pre","off",
                      "sample","to_sample","GRCh38","id","uORF_start","submitter_id","TCGA")          #加名字

query_trans_FPKM <- GDCquery(
  project = paste("TCGA",name,sep="_"),
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - FPKM"
)
GDCdownload(query_trans_FPKM,method = "api",files.per.chunk = 6)
expdate<-GDCprepare(query = query_trans_FPKM)
count_matrix_FPKM=assay(expdate)
count_matrix<-as.data.frame(count_matrix_FPKM)                                             #FPKM的矩阵

gene<-read.table("Homo_sapiens.GRCh38.104.gtf.geneInfo",sep="\t",header=T)               #导入gene位置信息
colnames(gene)[2]<-"tran_name"
over_uORF<-merge(gene[,(1:2)],getuORF,by="tran_name")                                        #将tran—name作为合并标准
getuORF_FPKM<-merge(over_uORF,count_matrix,by="gene_id")                                    #合并表达矩阵和突变信息
AS<-getuORF_FPKM                                                                            #
col1<-ncol(getuORF_FPKM)                                                                    #矩阵多少列
row1<-nrow(getuORF_FPKM)                                                                    #矩阵多少行
as<-getuORF_FPKM
getuORF_FPKM$A01<-NA
for (m in 1:row1){                                                                          #按行数循环
  row_o<-getuORF_FPKM$sample[m]                                                             #第m行的sample
  row_n<-strsplit(row_o,'-')[[1]][1:4]                                                      #前4个字符
  row_n1<-paste(row_n,collapse = '-')                                                       #合并

  for (i in 1:col1){                                                                        #按列循环
    col_o<-colnames(getuORF_FPKM)[i]                                                        #第列的列名
    col_n<-strsplit(col_o,'\\.')[[1]][1:4]                                                  #前4个字符
    col_n1<-paste(col_n,collapse = '-')
    ifelse(row_n1==col_n1 ,getuORF_FPKM$A01[m]<-getuORF_FPKM[m,i], getuORF_FPKM$A01[m]<-getuORF_FPKM$A01[m])         #列名和行sample相等，输出表达量在A01这列上
  }
}
AS<-getuORF_FPKM[,-(1:14)]                                                                  #删前14列，基因组突变数据
AS<-AS[,-ncol(AS)]                                                                          #删除最后一列，A01的结果

col1<-ncol(AS)
row1<-nrow(AS)
as<-AS
for (i in 1:col1){                                                                          #列循环
    col_o<-colnames(AS)[i]                                                                  #第i列列名
    col_n<-strsplit(col_o,'\\.')[[1]][4]                                                    #第4个字符
    col_n1<-strsplit(col_n,'')[[1]][1:2]                                                    #数字是多少
    col_n2<-paste(col_n1,collapse = '')
    col_n2<-as.numeric(col_n2)
    ifelse( col_n2 < 10 ,AS[row1+1,i]<- -100, AS[row1+1,i]<-  -200)                         #标记出是癌症的或癌旁的
}

as <- AS[,which(AS[nrow(AS),] == -100)]                                                     #将癌症的表达量提出
as<-as[-nrow(AS),]                                                                          #删除最后标记的那行
as$A01P<-rowSums(as)                                                                        #每行求和至A01P
getuORF_FPKM$A01P<-as$A01P                                                                  #然后将结果转至总表
getuORF_FPKM$A01P<-(getuORF_FPKM$A01P - getuORF_FPKM$A01)/(ncol(as)-2)                      #计算癌症表达量除了自己这个之外的平均值

AB<-getuORF_FPKM[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,ncol(getuORF_FPKM)-1,ncol(getuORF_FPKM))]     #前面突变信息和突变基因表达量以及平均表达量

as <- AS[,which(AS[nrow(AS),] == -200)]                                                     #将正常的表达量提出
as<-as[-nrow(AS),]                                                                          #删除最后标记的那行
as$A11P<-rowSums(as)                                                                        #每行求和至A11P
AB$A11P<-as$A11P/(ncol(as)-1)                                                               #计算平均值

AB$LFC_M_N<-log2((AB$A01+0.1)/(AB$A11P+0.1))                                                #突变和正常LFC
AB$LFC_M_T<-log2((AB$A01+0.1)/(AB$A01P+0.1))                                                #突变和肿瘤LFC
AB$LFC_T_N<-log2((AB$A01P+0.1)/(AB$A11P+0.1))                                               #肿瘤和正常LFC



write.table(AB,paste(name,"/",name,".getuORF_FPKM_LFC",sep=""),sep="\t",col.names = T,row.names=F,quote = F)
