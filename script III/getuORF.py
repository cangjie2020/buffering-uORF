from Bio import SeqIO
import argparse
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

parser = argparse.ArgumentParser(description='得到产生的uORF')
parser.add_argument('-outdir',  help='输出文件夹')
parser.add_argument('-geneInfo',  help='5UTR长度')
parser.add_argument('-t',  help='线程')
args = parser.parse_args()
def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)
geneInfo = args.geneInfo
outdir = args.outdir
chrs = list(range(1,22))
chrs.extend(['X','Y'])

def Get_uORF(geneInfo,snp_site, outfile):
    cdna=open("Homo_sapiens.GRCh38.cdna.all.fa")
    maf=open(snp_site)
    geneInfo =open(geneInfo)                                               #trans_site2.txt
    fh = open(outfile, "w")
    Distance={}
    ATG={}
    snpPos2 = pd.read_csv(gtf2,sep='\t',header=None,index_col =None)
    for line in maf:                                                     #循环每个突变
        line1 = line.split("\t")
        for mm in snpPos2.index:
            if line1[0] == snpPos2.loc[mm,][1]:
                if int(line1[2]) > int(snpPos2.loc[mm,][5]):
                    continue
                if line1[0] in Distance:                                         #转录本ID为键，整行为值
                    Distance[line1[0]].append(line)
                else:
                    Distance[line1[0]] = [line]
    DIS = Distance.keys()
    seqs = SeqIO.parse(cdna,'fasta')                                     #提取cDNA信息
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
                if distance_site < 0 :
                    continue
                if distance_site > len(seq):                             # distance_site要符合要求，在范围内
                    continue
                if to_snp == "A":                                        #如果突变成的A，则看后面两个碱基是不是TG
                    TG= seq[distance_site :distance_site +2]             #distance_site位置要减一才是突变碱基位置
                    if TG == "TG":
                        print(distance + "\t" + str(distance_site) ,file=fh)
                elif to_snp == "T":
                    if distance_site < 2:                                #至少在第二个碱基发生此突变
                        continue
                    AG= seq[distance_site-2] + seq[distance_site]        #前一个碱基和后一个碱基
                    if AG == "AG":
                        print(distance + "\t" + str(distance_site -1) ,file=fh)
                elif to_snp == "G":
                    if distance_site < 3:                                #至少在第三个碱基发生此突变
                        continue
                    AT= seq[distance_site-3:distance_site-1]             #前两个碱基
                    if AT == "AT":
                        print(distance + "\t" + str(distance_site - 2) ,file=fh)
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
        out_basename = "%s/%s.snv_getuORF.tmp" % (outdir,cc2)
        gtf_m = ["awk '$13 == \"%s\"' %s > %s/%s.geneInfo.tmp" % (cc, geneInfo, outdir,cc2)]
        subprocess.check_call(gtf_m, shell=True)
        geneInfo="%s/%s.geneInfo.tmp"% (outdir,cc2)
        result = pool.apply_async(Get_uORF,(geneInfo,snp_site,out_basename))
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
    Args_m = ['cat %s/chr*.snv_getuORF.tmp > snv_getuORF.txt'% (outdir)]
    subprocess.check_call(Args_m, shell=True)
    print(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
    print("getuORF annotation finished!!!")
