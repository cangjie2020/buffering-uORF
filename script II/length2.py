import argparse
import numpy as np
parser = argparse.ArgumentParser(description='计算片段长度信息')
parser.add_argument('-gene_info_file', help='输入转录组比对的Sam文件')
args = parser.parse_args()
ff = open(args.gene_info_file)                         #之前产生的sam文件
num=[]
for line in ff :
    line1 = line.split('\t')
    length=len(line1[9])                               #碱基序列长度
    num.append(int(length))                            #吧碱基序列长度数据加入列表
length_std=np.std(num,ddof=1)                          #计算标准差
print(length_std)
