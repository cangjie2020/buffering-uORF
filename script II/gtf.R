
args <- commandArgs(trailingOnly = TRUE)                #导入文件
path_in1<-args[1]
path_in2<-args[2]
path_in3<-args[3]                                               #导入所有正常细胞RSEM之后表达量的结果

table1<-read.table(path_in1,sep="\t",header = T)
table2<-read.table(path_in2,sep="\t",header = T)
table3<-read.table(path_in3,sep="\t",header = T)                #导入
gtf <-read.table("/data/caizixin/uorf/RSEM/Homo_sapiens.GRCh38.104.gtf.geneInfo",sep="\t",header = T)    #导入gene类别
GenePred <-read.table("/data/caizixin/uorf/Homo_sapiens.GRCh38.104.GenePred",sep="\t",header = F)        #注释文件前奏

table1.1<-table1[,c(1,2,5)]                                      #只要1,2,5列的
table2.1<-table2[,c(1,2,5)]
table3.1<-table3[,c(1,2,5)]

all=merge(table1.1,table2.1,by=c("transcript_id","gene_id"))
all=merge(all,table3.1,by=c("transcript_id","gene_id"))           #将结果合在一起

all$sum<-all[,3] + all[,4] + all[,5]                              #计算每个转录本的总表达量
gtf<-gtf[,c(2,11)]                                                #提取含有gene类别的信息
colnames(gtf)[1]<-"transcript_id"                                 #一样的列名
all<-merge(all,gtf,by="transcript_id")
all<-all[which(all$transcript_biotype == "protein_coding"),]      #只要编码蛋白的转录本
all<-all[order(-all[,6]),]                                        #排序
all<-all[!duplicated(all[,2]),]                                   #去除重复
A<-all[,c(2,6,1)]
colnames(GenePred)[1]<-"transcript_id"
A<-merge(A,GenePred,by="transcript_id")
A<-A[,-c(2,3)]                                                    #符合GenePred格式

write.table(A, file="thebestrans.GenePred" ,sep = "\t",row.names = F,quote = F,col.names =F)
