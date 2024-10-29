import sys
import pandas as pd
import re
import gzip
from Bio import SeqIO

fh = open("Homo_sapiens.noncanonical_uORF.tiv.txt", "w")
init_codon = ['ATG']
stop_codon = ['TAA', 'TAG', 'TGA']

gene_info_file = "Homo_sapiens.GRCh38.104.gtf.geneInfo"  # 'data/Ensembl_r96/Homo_sapiens.GRCh38.96.gtf.geneInfo'
gene_info = pd.read_csv(gene_info_file, sep='\t', header=0)            #将文件转换成表格形式
txfine = gene_info[(gene_info.utr5_len > 0) & (gene_info.transcript_biotype == 'protein_coding')]    #取5utr长度大于0的，且是protein_coding，变量矩阵中的列名前要加上变量矩阵的名字
utr5_len = {i:j for i, j in zip(txfine.tx_name, txfine.utr5_len)}          #取这两列，记得变量名的名字要写,就是每个转录本的名字和5utr的长度

cdna_path = "Homo_sapiens.GRCh38.cdna.all.fa.gz"                  #导入基因组序列信息
seqs = SeqIO.parse(gzip.open(cdna_path, 'rt'),'fasta')           #建立基因组序列信息,得到Seq是碱基序列，ID为转录本ID，还有Description,负链的碱基会反向互补

split_id = gene_info_file.find('.104.gtf') >= 0                 #检查是不是这个注释文件

for record in seqs:

    tx_name = record.id                   #转录本id
    #print(tx_name)
    if split_id:                          #如果split_id大于0，是正确的
        tx_name = tx_name.split('.')[0]   #不要.之后的号码，就是第几号转录本，因为要和utr5_len匹配
    if tx_name not in utr5_len:
        continue                          #就是只要相同的，如果满足tx_name不等于utr5_len，则往下继续，但是之后就没了，若tx_name等于utr5_len则循环下一个
    tx_utr5_len = utr5_len[tx_name]       #取tx_name在utr5_len对应的长度
    seq = record.seq                      #转录本碱基序列
    for single in init_codon:
        augs = re.finditer(single, str(seq))            #找转录本序列中的ATG
        for aug in augs:
            t_start = aug.start() + 1  # 起始密码子开始的位置
            if t_start > tx_utr5_len:  # 不要起始密码子在5utr之后的
                continue
            j = aug.start() + 3       # 起始密码子结束位置
            last_tidx = len(seq) - 3  # 碱基片段长度减3，去掉起始密码子的总长度
            stop_found = False
            while j <= last_tidx:         # 当起始密码子位置在转录本内的时候，然后分成3个一组
                triplet = seq[j:(j + 3)]  # 三个一组
                if triplet in stop_codon:  # 如果有一组是终止密码子
                    stop_found = True  # 显示找到了终止密码子
                    break  # 结束
                j = j + 3  # 再循环下3个碱基
            t_end = j + 3 if stop_found else len(seq)  # uorf结束的位置
            out = [tx_name, t_start, t_end, stop_found, tx_utr5_len, single]  # 输出的结果
            print('\t'.join([str(k) for k in out]))
