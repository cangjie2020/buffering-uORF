args <- commandArgs(trailingOnly = TRUE)
path_in<-args[1]                             #siteandtrans.txt
path_in2<-args[2]                            #SRR${i}.13-.sam
path_OUT<-args[3]                            #CDS.count.txt
path_OUT2<-args[4]                           #uORF.gene.count.txt
path_OUT3<-args[5]                           #uORF.count.txt
path_OUT4<-args[6]                           #all.phase.pdf

start1<-as.numeric(args[7])
end1<-as.numeric(args[8])
psite1<-as.numeric(args[9])
start2<-as.numeric(args[10])
end2<-as.numeric(args[11])
psite2<-as.numeric(args[12])
start3<-as.numeric(args[13])
end3<-as.numeric(args[14])
psite3<-as.numeric(args[15])


A<-read.table("/data/caizixin/uorf/RSEM/Homo_sapiens.GRCh38.104.gtf.geneInfo",sep="\t",header=T)
A<-A[,c(1,2,5,6,7,9,14,10)]
B<-read.table(path_in,sep="\t",header=F)
colnames(B)<-c("read.name","chrom","weizhi","long","tx_name")               #B
C<-merge(A,B,by="tx_name")
D<-read.csv(path_in2,sep="\t",comment.char = "",header=F)
D<-D[c(1,3,4)]
colnames(D)<-c("read.name","tx_name","trans.site")                          #D
E<-merge(D,C,by=c("read.name","tx_name"))

E$CDS.start<-E$utr5_len + 1                                                 # 转录本CDS起始位置
E$CDS.end<-E$utr5_len + E$cds_len                                           # 转录本cds结束位置
E<-E[-c(5,7,11)]
for  (i in seq(start1,end1,1)){                                             #相同的Psite，从多长到多长
  R<-E[which(E$long == i ),]
  R$trans.site<-R$trans.site + psite1                                       #将其右移Psite碱基，找到P位点
  assign(paste("Psite",i,sep = "_"),R)
}
for  (i in seq(start2,end2,1)){
R<-E[which(E$long == i ),]
R$trans.site<-R$trans.site + psite2
assign(paste("Psite",i,sep = "_"),R)
}

for  (i in seq(start3,end3,1)){
R<-E[which(E$long == i ),]
R$trans.site<-R$trans.site + psite3
assign(paste("Psite",i,sep = "_"),R)
}
i = start1                                                                 #起始reads长度
u = as.numeric(start1)+1                                                   #第二reads长度
c = as.numeric(start1)+2                                                   #第三reads长度
gene_cds_count<-rbind(get(paste("Psite",i,sep = "_")),get(paste("Psite",u,sep = "_")))     # 行合并
for (i in seq(c,end3,1)){
  d<-paste("Psite",i,sep = "_")
  gene_cds_count<-rbind(gene_cds_count,get(d))                                              #从第三reads长度到最后
}                                                                                           #全部行合并
#########################################################################################CDS
SRRz<-gene_cds_count[which(gene_cds_count$trans.site >= gene_cds_count$CDS.start),]         #筛选出校正后的CDS内的reads
SRRz<-SRRz[which(SRRz$trans.site <= SRRz$CDS.end),]
B_tab = table(SRRz$tx_name)                                                                 #统计在CDS内的reads有多少
B_new<-as.data.frame(B_tab)

CDS.LONG<-A[c(1,2,3,4,6,8)]                                                                 #geneInfo
colnames(B_new)[1]<-"tx_name"
out_CDS<-merge(CDS.LONG,B_new,by="tx_name")                                                 #带上转录本的位置信息输出
write.table(out_CDS,file=path_OUT,sep = "\t",row.names = F,quote = F,col.names = T)
########################################################################################################uORF
uORF<-read.table("/data/caizixin/uorf/RSEM/Homo_sapiens.uORF.TRUE.tiv.txt",sep="\t",header=F)

uORF<-uORF[,c(1,2,3)]
uORF$uORF_ID<-paste(uORF$V1, uORF$V2,sep="_")
colnames(uORF)<-c("tx_name","uORF.start","uORF.end","uORF_ID")                            #输入文件
gene_uORF_count<-merge(gene_cds_count,uORF,by="tx_name")                                  #uORF情况
SRRz<-gene_uORF_count[which(gene_uORF_count$trans.site >= gene_uORF_count$uORF.start),]   #留下在uORF范围内的uORFreads
SRRz<-SRRz[which(SRRz$trans.site <= SRRz$uORF.end),]
B_tab = table(SRRz$gene_id)                                                               #在一个gene中的uORFreads数
B_new<-as.data.frame(B_tab)
uORF$long<-uORF$uORF.end - uORF$uORF.start +1                                             #统计个数和长度
write.table(B_new, file=path_OUT2,sep = "\t",row.names = F,quote = F,col.names = T)       #输出

A_tab = table(SRRz$uORF_ID)                                                               #在一个uORF中的uORFreads数
A_new<-as.data.frame(A_tab)

colnames(A_new)[1]<-"uORF_ID"
out.uORF<-merge(A_new,uORF,by="uORF_ID")                                                  #结合uorf数目和相关信息
out.uORF<-out.uORF[,c(1,2,3,6)]
out.uORF<-merge(out.uORF,CDS.LONG,by="tx_name")                                           #
write.table(out.uORF, file=path_OUT3,sep = "\t",row.names = F,quote = F,col.names = T)

library(ggplot2)
gene_cds_count$ditance<-gene_cds_count$trans.site - gene_cds_count$CDS.start              #总的在CDS内的reads的三碱基周期性图
SRRz<-gene_cds_count[which(gene_cds_count$trans.site >= gene_cds_count$CDS.start),]
SRRz<-SRRz[which(SRRz$trans.site <= SRRz$CDS.end),]

threshold <- as.factor(ifelse(SRRz$ditance %% 3 == 0,"1",ifelse(SRRz$ditance %% 3 == 1,"2",ifelse(SRRz$ditance %% 3 == 2,"3","0"))))
P<-as.data.frame(threshold)
colnames(P)[1]<-"phase"
pdf(path_OUT4 ,family="Helvetica",width=9,height=10)
p<-ggplot(P, aes(x=phase,fill=phase))+ geom_bar()+xlab("phase")+ylab("count") +   theme_bw()+theme(axis.text.x = element_text(size = 30),)+ theme(axis.text.y = element_text(size = 30),)+theme(axis.title.x = element_text(size = 30, vjust = 0.5, hjust = 0.5))+theme(axis.title.y = element_text(size = 30, vjust = 0.5, hjust = 0.5)) + theme(panel.grid.major = element_blank()) +  theme(panel.grid =element_blank()) + theme(legend.key=element_blank()) + theme(panel.border = element_blank()) + theme(axis.line = element_line(size=1, colour = "black"))+  theme(legend.title = element_text(size=30, face="bold")) +theme(legend.text = element_text( size = 30))
print(p)
dev.off()
