from collections import defaultdict
import pandas as pd
import sys
import os
from time import clock
from numpy import  *
import re
import multiprocess
import time
import subprocess
import argparse
import signal
import shutil
import operator
from functools import reduce
from Bio import SeqIO
import argparse
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='识别突变产生的uORF2')
parser.add_argument('-outdir',  help='输出文件夹')
parser.add_argument('-t',  help='线程')
args = parser.parse_args()
def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

outdir = args.outdir
chrs = list(range(1,23))
chrs.extend(['X','Y'])

def Get_inuORF(snp_site, outfile):
    cdna=open("Homo_sapiens.GRCh38.cdna.all.fa")                                                        # cDNA文件
    maf=open(snp_site)
    fh = open(outfile, 'w')                                                                              # trans_site2.txt
    uORF=open("Homo_sapiens.uORF.TRUE.tiv.txt")                                                         #uORF位置
    aa=open("AA.txt")                                                                                   #氨基酸类别
    Distance={}
    inuORF={}
    AAs={}
    for line in maf:
        line1 = line.split("\t")
        if line1[0] in Distance:
            Distance[line1[0]].append(line)                                                             #转录本ID为键，整行为值
        else:
            Distance[line1[0]] = [line]
    DIS = Distance.keys()
    for line in uORF:
        line1 = line.split("\t")
        if line1[0] in DIS:
            for m in range(len(Distance[line1[0]])):
                distance = Distance[line1[0]][m].strip("\n")
                distance_site = int(distance.split("\t")[2])                                            #突变位点
                if int(line1[1]) <= distance_site < int(line1[2]):                                      #uORF的起始到终止位点之间有突变位点
                    if line1[0] in inuORF:
                        inuORF[line1[0]].append(distance + "\t" + line1[1] + "\t" + line1[0] + "_" + line1[1])    #字典的值为此
                    else:
                        inuORF[line1[0]]=[distance + "\t" + line1[1] + "\t" + line1[0] + "_" + line1[1]]
    inuorf=inuORF.keys()

    seqs = SeqIO.parse(cdna,'fasta')
    for record in seqs:
        tx_name = record.id
        tx_name = tx_name.split('.')[0]                                                                  #相同的转录本id
        if tx_name in inuorf :
            seq = record.seq                                                                             #Cdna序列
            for m in range(len(inuORF[tx_name])):
                uorf = inuORF[tx_name][m].strip("\n")
                uorf1=uorf.split("\t")
                uorf_site=int(uorf1[6])                                                                 #uORF的起始位点
                snp_site=int(uorf1[2])                                                                   #突变位点
                phase = (snp_site - uorf_site)%3                                                         #突变的位点与uORF起始位点的距离除以3取余
                if phase == 0:                                                                           #说明是三个碱基中开头的那个碱基
                    AA = seq[snp_site - 1:snp_site + 2]                                                  #原始碱基就是snp_site - 1以及后面两个
                    to_AA = uorf1[5] + seq[snp_site:snp_site + 2]                                        #突变后碱基是突变后的碱基加原来的后两个
                    if AA in AAs:                                                                        #生成字典
                        AAs[AA].append(uorf + "\t" + AA + "\t" + to_AA )
                    else:
                        AAs[AA] = [uorf + "\t" + AA + "\t" + to_AA ]
                elif phase == 1:                                                                         #说明是三个碱基中中间的那个碱基
                    AA = seq[snp_site - 2:snp_site + 1]
                    to_AA = seq[snp_site - 2] + uorf1[5] + seq[snp_site]
                    if AA in AAs:
                        AAs[AA].append(uorf + "\t" + AA + "\t" + to_AA )
                    else:
                        AAs[AA] = [uorf + "\t" + AA + "\t" + to_AA ]
                elif phase == 2:                                                                         #说明是三个碱基中最后的那个碱基
                    AA = seq[snp_site - 3:snp_site]
                    to_AA = seq[snp_site - 3:snp_site - 1] + uorf1[5]
                    if AA in AAs:
                        AAs[AA].append(uorf + "\t" + AA + "\t" + to_AA )
                    else:
                        AAs[AA] = [uorf + "\t" + AA + "\t" + to_AA ]
    stop_codon = ['TAA', 'TAG', 'TGA']
    for line in aa:
            line1 = line.split("\t")
            l = line1[1].strip("\n")
            l = l.split(",")
            for m in range(len(l)):
                if l[m] in AAs.keys():                                                               #如果
                    uorf = AAs[l[m]]
                    for i in range(len(uorf)):
                        uorf2=uorf[i].split("\t")
                        if uorf2[8] == "ATG":
                            print(uorf[i] + "\t" + "deuORF" + "\t",file=fh)
                            continue
                        if uorf2[8] not in stop_codon  and uorf2[9] in stop_codon:
                            print(uorf[i] + "\t" + "nonsense mutation" +  "\t",file=fh)
                            continue
                        if uorf2[9] in l:
                            print(uorf[i] + "\t" + "Synonymous mutation" +  "\t" ,file=fh)
                        else:
                            print(uorf[i] + "\t" + "Missense mutation"  +  "\t",file=fh)
    fh.close()
if __name__ == '__main__':
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    results_q = []
    tmp_ReadOut = []
    t0=clock()
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
    pool = multiprocess.Pool(int(args.t),init_worker)
    t1=clock()
    print('Split chr file')
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
    for cc in chrs:
        cc2='chr' + str(cc)
        snp_site = "%s/%s.snv_gene.tmp" % (outdir,cc2)
        out_basename = "%s/%s.snv_inuORF.tmp" % (outdir,cc2)
        result = pool.apply_async(Get_inuORF,(snp_site,out_basename))
        results_q.append(result)
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
    try:
        print("Waiting 1 seconds")
        time.sleep(5)
    except KeyboardInterrupt:
        print("Caught KeyboardInterrupt, terminating workers")
        pool.terminate()
        pool.join()
    else:
        print("Quitting normally")
        pool.close()
        pool.join()
        print('Pool finished')
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
    Args_m = ['cat %s/chr*.snv_inuORF.tmp > %s_snv_inuORF.txt'% (outdir,outdir)]
    subprocess.check_call(Args_m, shell=True)
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
    print("inuORF annotation finished!!!")
