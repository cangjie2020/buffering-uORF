
args <- commandArgs(trailingOnly = TRUE)
path_in<-args[1]
path_out<-args[2]
long<-read.table(path_in,sep="\t")

long<-long[which(long[,4] < 36),]
library(ggplot2)
pdf(path_out,family="Helvetica",width=10,height=9)
p<-ggplot(long, aes(x=V4)) + geom_bar(fill='blue')+xlab("length")+ylab("count") +
  theme_bw()+theme(axis.text.x = element_text(size = 30 ),)+
  theme(axis.text.y = element_text(size = 30),) +
  theme(axis.title.x = element_text(size = 30, vjust = 0.5, hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 30, vjust = 0.5, hjust = 0.5)) +
  theme(panel.grid.major = element_blank()) +  theme(panel.grid =element_blank()) +
  theme(legend.key=element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(size=1, colour = "black"))+
  theme(legend.title = element_text(size=30, face="bold")) +
  theme(legend.text = element_text( size = 30))
print(p)
dev.off()

#################注释版####################
args <- commandArgs(trailingOnly = TRUE)
path_in<-args[1]                           #导入之前的siteandtrans.txt
path_out<-args[2]                          #导出名称
long<-read.table(path_in,sep="\t")

long<-long[which(long[,4] < 36),]          #只要长度小于36的结果
library(ggplot2)
pdf(path_out,family="Helvetica",width=10,height=9)
p<-ggplot(long, aes(x=V4)) + geom_bar(fill='blue')+xlab("length")+ylab("count") +
  theme_bw()+theme(axis.text.x = element_text(size = 30 ),)+
  theme(axis.text.y = element_text(size = 30),) +
  theme(axis.title.x = element_text(size = 30, vjust = 0.5, hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 30, vjust = 0.5, hjust = 0.5)) +
  theme(panel.grid.major = element_blank()) +  theme(panel.grid =element_blank()) +
  theme(legend.key=element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(size=1, colour = "black"))+
  theme(legend.title = element_text(size=30, face="bold")) +
  theme(legend.text = element_text( size = 30))
print(p)
dev.off()
