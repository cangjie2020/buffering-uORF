args <- commandArgs(trailingOnly = TRUE)
path_in<-args[1]                                                                  # siteandtrans.txt
path_in2<-args[2]                                                                 # SRR${i}.13-.sam
path_out<-args[3]                                                                 # reads终止长度（最短是多少）
start <- as.numeric(args[4])                                                      #转数字
end <- as.numeric(args[5])
library(ggplot2)
sam<-read.table(path_in,sep="\t",header=F)                                          #siteandtrans.txt信息
colnames(sam)<-c("read.name","chrom","site","length","trans.name")                  #每一列是什么
CDS.site<-read.csv("/data/caizixin/uorf/RSEM/Homo_sapiens.GRCh38.104.gtf.geneInfo",sep="\t",header = T)           #导入转录本类别
CDS.site2<-CDS.site[,c(1,2,6,7,9,10,13,14)]

colnames(CDS.site2)<-c("gene_id","trans.name","CDS.long","5utr.long","gene_name","gene_biotype","chrom","strand")  #选择这几列
trans<-read.csv(path_in2,sep="\t",comment.char = "",header=F)

trans2<-trans[,c(1,3,4,6)]
colnames(trans2)<-c("read.name","trans.name","trans.site","long")                   #sam留下来的信息，主要是在转录本上的位置
A<-merge(trans2,sam,by=c("read.name","trans.name"))                                 #结合转录本位置和有多少reads是满足我们的要求的reads
B<-A[,c(1,2,3,7)]                                                                   #"read.name","trans.name"，"trans.site"，"length"

CDS.site2<-CDS.site2[which(CDS.site2$gene_biotype == "protein_coding"),]           #之前做过
CDS.site2$CDS.start<-CDS.site2$`5utr.long` + 1                                     #CDS起始位置是5utr长度+1
CDS.site2$CDS.end<-CDS.site2$`5utr.long` + CDS.site2$CDS.long                      #CDS结束位置是5utr长度+CDS长度
CDS.site3<-CDS.site2[,-c(3,4,5,6,7,8)]
C<-merge(B,CDS.site3,by="trans.name")
C$distance<-C$trans.site - C$CDS.start                                             #reads起始位点和CDS起始位点的距离


u=path_out
for  (i in seq(start,end,1)){                                                      #从start位置开始循环（最短序列长度）

  A<-C[which(C$length == i),]                                                      #先只要这个长度的
  A<-A[which(A$distance > -30),]                                                   #距离CDS在-30到60之间
  A<-A[which(A$distance < 60 ),]

  threshold <- as.factor(ifelse(A$distance %% 3 == 0,"1",ifelse(A$distance %% 3 == 1,"2",ifelse(A$distance %% 3 == 2,"3","0"))))  #三碱基周期
  NAME<-paste(u,i,"start.pdf",sep = "_")                                           #为产生的图命名
  pdf(NAME,family="Helvetica",width=10,height=9)                                   #绘制三碱基周期性图，在一个范围区间内
  p<-ggplot(A, aes(x=distance,fill=threshold))+ geom_bar()+xlab("length")+ylab("count") +
            scale_x_continuous(breaks=seq(-30, 50, 10))+ theme_bw()+
            theme(axis.text.x = element_text(size = 30,angle = 45,vjust = 1.2,hjust = 1.2 ),)+
            theme(axis.text.y = element_text(size = 30),)+
            theme(axis.title.x = element_text(size =30, vjust = 0.5, hjust = 0.5))+
            theme(axis.title.y = element_text(size = 30, vjust = 0.5, hjust = 0.5)) +
            theme(panel.grid.major = element_blank()) +
            theme(panel.grid =element_blank()) +
            theme(legend.key=element_blank()) +
            theme(panel.border = element_blank()) +
            theme(axis.line = element_line(size=1, colour = "black"))+
            theme(legend.title = element_text(size=20, face="bold")) +
            theme(legend.text = element_text( size = 20))
  print(p)
  dev.off()

P<-as.data.frame(threshold)
colnames(P)[1]<-"phase"
NAME<-paste(u,i,"p.pdf",sep = "_")
pdf(NAME ,family="Helvetica",width=9,height=10)                                      #绘制三碱基周期性图，只有三条柱子
p<-ggplot(P, aes(x=phase,fill=phase))+ geom_bar()+xlab("phase")+ylab("count") +
        theme_bw()+theme(axis.text.x = element_text(size = 30),)+
        theme(axis.text.y = element_text(size = 30),)+
        theme(axis.title.x = element_text(size = 30, vjust = 0.5, hjust = 0.5))+
        theme(axis.title.y = element_text(size = 30, vjust = 0.5, hjust = 0.5)) +
        theme(panel.grid.major = element_blank()) +  theme(panel.grid =element_blank()) +
        theme(legend.key=element_blank()) + theme(panel.border = element_blank()) +
        theme(axis.line = element_line(size=1, colour = "black"))+
        theme(legend.title = element_text(size=30, face="bold")) +
        theme(legend.text = element_text( size = 30))
print(p)
dev.off()
}
