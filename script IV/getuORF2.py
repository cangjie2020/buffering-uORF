import argparse

parser = argparse.ArgumentParser(description='识别突变产生的uORF2')
parser.add_argument('-i', help='输入文件')
parser.add_argument('-o',  help='输出文件名以及路径')
args = parser.parse_args()
ATG={}
from Bio import SeqIO
cdna=open("Homo_sapiens.GRCh38.cdna.all.fa")                                #cDNA序列信息
input=open(args.i)                                                          #getuORF1.txt
fh = open(args.o, "w")
for line in input:
    line1 = line.split("\t")
    if line1[0] in ATG:
        ATG[line1[0]].append(line)                                          #转录本ID为键 ，值为一整行
    else:
        ATG[line1[0]] = [line]
atg = ATG.keys()
seqs = SeqIO.parse(cdna,'fasta')
for record in seqs:                                                         #循环序列
    tx_name = record.id
    tx_name = tx_name.split('.')[0]
    stop_codon = ['TAA', 'TAG', 'TGA']                                      #定义终止密码子
    if tx_name in atg:
        seq = record.seq
        for m in range(len(ATG[tx_name])):
            distance = ATG[tx_name][m].strip("\n")
            j = int(distance.split("\t")[10])+2                             # 产生uORF的起始位点加三，就是下一个三碱基的起始位点
            last_tidx = len(seq)
            stop_found = False
            while j <= last_tidx:                                           #当终止密码子还在序列内的
                triplet = seq[j:(j + 3)]                                    #三个碱基一次循环
                if triplet in stop_codon:                                   #如果产生了循环中有终止密码子
                    stop_found = True                                       #表示出来
                    break
                j = j + 3                                                   #否则继续加三
            t_end = j + 3 if stop_found else len(seq)                       #输出uORF终止密码子位置
            uORF_len = t_end - int(distance.split("\t")[10]) + 1            #uORF长度
            out = [tx_name, distance.split("\t")[10], t_end, stop_found,uORF_len]
            print('\t'.join([str(k) for k in out]),file=fh)
fh.close()
