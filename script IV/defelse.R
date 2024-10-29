
args <- commandArgs(trailingOnly = TRUE)
i<-args[1]

A<-read.table(paste(i,"/",i,".getuORF2.txt",sep=""),sep="\t",header=F)
B<-read.table(paste(i,"/",i,".getuORF1.txt",sep=""),sep="\t",header=F)
C<-cbind(A,B)
C<-C[which(C[,4]=="True"),]                                                                       #只要是TRUE的行
D<-C[,(6:18)]
write.table(D,paste(i,"/",i,".getuORF3.txt",sep=""),sep="\t",col.names = F,row.names=F,quote = F)
