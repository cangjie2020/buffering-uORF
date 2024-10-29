

import argparse

parser = argparse.ArgumentParser(description='转变在5UTR中gene位置为转录本位置')
parser.add_argument('-i', help='输入文件')
args = parser.parse_args()
maf=open(args.i)
GenePred=open("./Homo_sapiens.GRCh38.104.GenePred")                           #提供外显子位置
B= {}
for line in maf:                                                              #按行循环maf文件
    if line[0] == "#":                                                        #删除注释文件
        continue
    tran = line.split("\t")[48]                                               #转录本ID
    C = line.split("\t")[6]                                                   #突变位点
    if line.split("\t")[8]=="5'UTR":                                          #位置在5'UTR
        if tran in B:                                                         #以转录本ID为键，一整行为值
            B[tran].append(line)
        else:
            B[tran]=[line]
    if line.split("\t")[8] == "Translation_Start_Site":                       #位置在Translation_Start_Site
        if tran in B:
            B[tran].append(line)
        else:
            B[tran] = [line]
key_list = B.keys()                                                           #第一个字典键为转录本ID，值为一整行
#############################################################################第二层
for line1 in GenePred:
    if line1.split("\t")[0] in key_list:                                      #将GenePred中的转录本ID和之前的字典对应
        for m in range(len(B[line1.split("\t")[0]])):                         #将在转录本ID为键下的全部突变循环
            site = B[line1.split("\t")[0]][m]                                 #第m个突变
            refer_snp = site.split("\t")[10]                                  #参考序列碱基
            snp = site.split("\t")[11]                                        #正常组织中的碱基
            to_snp = site.split("\t")[12]                                     #癌细胞突变产生的碱基
            sample = site.split("\t")[15]                                     #癌细胞的样本名
            to_sample = site.split("\t")[16]                                  #正常细胞的样本名
            snp_site = int(site.split("\t")[5])                               #突变在染色体上的位置
            exon_num = line1.split("\t")[7]                                   #这个转录本外显子个数
            if len(to_snp) > 1:                                               #只要snp
                continue
            if len(snp) > 1:
                continue
            if len(refer_snp) > 1:
                continue
            if snp == "-":
                continue
            if to_snp == "-":
                continue                                                      #将除了snp都去掉
            if refer_snp == "-":
                continue
            if snp != refer_snp:                                             # 如果参考序列和癌旁不一样就不要此情况
                continue
            if line1.split("\t")[2] == "+":                                   #分成正链和负链的分别处理
                trans_start=line1.split("\t")[3]                              #转录本起始位置
                exon_start=line1.split("\t")[8]                               #外显子起始位置
                exon_end = line1.split("\t")[9]                               #外显子终止位置
                distance = 0
                for i in range(int(exon_num)):                                #循环外显子个数
                    start = int(exon_start.split(",")[i])                     #外显子的起始位置
                    end = int(exon_end.split(",")[i])                         #外显子的终止位置
                    long=end - start                                          #外显子长度

                    if start <= snp_site <= end:                              #突变位点在哪个转录本上
                        distance = distance + (snp_site - start)              #最后一个外显子上距离
                        distance_site= line1.split("\t")[0] + "\t" + str(snp_site) + "\t" + str(distance) + "\t" + str(line1.split("\t")[2]) + "\t" + snp +"\t" + to_snp +"\t" + sample +"\t" + to_sample + "\t" + site.split("\t")[3] + "\t" +  site.split("\t")[13]
                        print(distance_site)
                    #输出转录本ID，突变位点，突变位点距离转录本起始位点距离，正负链，突变前碱基，突变后碱基，突变后样本名，突变前样本名， GRCh38，dbSNP中的名字
                        continue
                    distance = distance + long                                #外显子距离累加
            if line1.split("\t")[2] == "-":                                   #负链的情况
                trans_start=line1.split("\t")[4]
                exon_start=line1.split("\t")[9]
                exon_end = line1.split("\t")[8]
                distance = 0
                for ii in range(int(exon_num),0,-1):                          #从最后一个外显子开始数
                    i= ii - 1                                                 #保证最后是0
                    start = int(exon_start.split(",")[i])                     #
                    end = int(exon_end.split(",")[i])
                    long=start - end                                          #外显子长度

                    if end <= snp_site <= start:
                        distance = distance + (start - snp_site +1)           #最后一个外显子上距离
                        if to_snp == "A":                                     #改变负链上突变完成后的碱基，
                            to_snp = "T"
                            distance_site=  line1.split("\t")[0] + "\t" + str(snp_site) + "\t" + str(distance) + "\t" + str(line1.split("\t")[2]) + "\t" + snp + "\t" +  to_snp +"\t" + sample +"\t" + to_sample + "\t" + site.split("\t")[3] + "\t" +  site.split("\t")[13]
                        elif to_snp == "T":
                            to_snp = "A"
                            distance_site=  line1.split("\t")[0] + "\t" + str(snp_site) + "\t" + str(distance) + "\t" + str(line1.split("\t")[2]) + "\t" + snp + "\t" +  to_snp +"\t" + sample +"\t" + to_sample + "\t" + site.split("\t")[3] + "\t" +  site.split("\t")[13]
                        elif to_snp == "C":
                            to_snp = "G"
                            distance_site=  line1.split("\t")[0] + "\t" + str(snp_site) + "\t" + str(distance) + "\t" + str(line1.split("\t")[2]) + "\t" + snp + "\t" +  to_snp +"\t" + sample +"\t" + to_sample + "\t" + site.split("\t")[3] + "\t" +  site.split("\t")[13]
                        elif to_snp == "G":
                            to_snp = "C"
                            distance_site=  line1.split("\t")[0] + "\t" + str(snp_site) + "\t" + str(distance) + "\t" + str(line1.split("\t")[2]) + "\t" + snp + "\t" +  to_snp +"\t" + sample +"\t" + to_sample + "\t" + site.split("\t")[3] + "\t" +  site.split("\t")[13]
                        print(distance_site)
                        continue
                    distance = distance + long
