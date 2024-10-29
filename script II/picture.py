
#为画图准备的前期工作
import argparse
parser = argparse.ArgumentParser(description='计算片段长度信息')
parser.add_argument('-gene_info_file', help='输入转录组比对的Sam文件')
parser.add_argument('-out',  help='输出文件名以及路径')
args = parser.parse_args()
ff = open(args.gene_info_file)
fh = open(args.out, "w")
for line in ff :
    line1 = line.split('\t')
    length=len(line1[9])                                        #序列长度
    trans = line1[17].split(':')[2]                             #提取结合到的转录本
    print(line1[0],line1[2],line1[3],length,trans.replace("\n", ""),sep="\t",file=fh)
    #输出分别是reads_id,染色体号，染色体位置，reads长度，转录本id
fh.close()
