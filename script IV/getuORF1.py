from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='转变在5UTR中gene位置为转录本位置')
parser.add_argument('-i', help='输入文件')
parser.add_argument('-o',  help='输出文件名以及路径')
args = parser.parse_args()

cdna=open("Homo_sapiens.GRCh38.cdna.all.fa")
maf=open(args.i)                                                      #trans_site2.txt
fh = open(args.o, "w")
Distance={}
ATG={}
for line in maf:                                                     #循环每个突变
    line1 = line.split("\t")
    if line1[0] in Distance:                                         #转录本ID为键，整行为值
        Distance[line1[0]].append(line)
    else:
        Distance[line1[0]] = [line]
DIS = Distance.keys()
seqs = SeqIO.parse(cdna,'fasta')                                     #提取cDNA信息
name_TCGA = args.i.split("/")
name_TCGA = name_TCGA[1]
for record in seqs:                                                  #循环每个转录本
    tx_name = record.id
    tx_name = tx_name.split('.')[0]                                  #使cDNA和文件的转录本ID相同
    if tx_name in DIS :
        seq = record.seq                                             #转录本序列
        for m in range(len(Distance[tx_name])):                      #循环一个转录本id内的所有突变
            distance = Distance[tx_name][m].strip("\n")              #
            distance_site = int(distance.split("\t")[2])
            strand = distance.split("\t")[3]
            snp = distance.split("\t")[4]
            to_snp = distance.split("\t")[5]
            sample = distance.split("\t")[6]                         #突变后样本
            to_sample = distance.split("\t")[7]                      #突变前样本
            name = distance.split("\t")[6]
            name1 = name.split("-")[0:3]
            name2 = '-'.join([str(k) for k in name1])                #突变后样本名前三号
            if distance_site < 0 :
                continue
            if distance_site > len(seq):                             # distance_site要符合要求，在范围内
                continue
            if to_snp == "A":                                        #如果突变成的A，则看后面两个碱基是不是TG
                TG= seq[distance_site :distance_site +2]             #distance_site位置要减一才是突变碱基位置
                if TG == "TG":
                    print(distance + "\t" + str(distance_site) + "\t" + name2  + "\t" + name_TCGA,file=fh)
            elif to_snp == "T":
                if distance_site < 2:                                #至少在第二个碱基发生此突变
                    continue
                AG= seq[distance_site-2] + seq[distance_site]        #前一个碱基和后一个碱基
                if AG == "AG":
                    print(distance + "\t" + str(distance_site -1) + "\t" + name2 + "\t" + name_TCGA,file=fh)
            elif to_snp == "G":
                if distance_site < 3:                                #至少在第三个碱基发生此突变
                    continue
                AT= seq[distance_site-3:distance_site-1]             #前两个碱基
                if AT == "AT":
                    print(distance + "\t" + str(distance_site - 2) + "\t" + name2 + "\t" + name_TCGA,file=fh)
fh.close()
