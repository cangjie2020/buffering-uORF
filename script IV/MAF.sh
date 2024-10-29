for i in $(cat file.txt)
do                                                                 #样本名，TCGA类型
mkdir ${i}                                                             #建立该样本文件夹，该类型
tar=${2}                                                               #压缩文件
tar xzvf ${tar} -C ./${i}                                              #解压下载的对应TCGA类型的MAF文件，解压到该样本文件夹下
zip=*"maf.gz"                                                          #文件夹内全部压缩文件
gunzip ./${i}/*/${zip}                                                 #解压之前解压出来的压缩文件，在文件夹下所有的压缩文件
ls ./${i}/*/*".maf" > nema                                             #将全部的文件名制成文件
rm ./${i}/${i}.trans_site.txt                                          #因为这里之后是>>生成文件，所以如果一开始有trans_site.txt不行
for line in $(cat nema)                                                #分成不同文件操作
do
        echo ${line}
        python3 gene_to_trans.py -i ${line}  >>./${i}/${i}.trans_site.txt
        # 提取MAF文件中的内容，转变在5UTR中gene位置为转录本位置，以及替换负链突变之后碱基
        # 输出转录本ID，突变位点，突变位点距离转录本起始位点距离，正负链，突变前碱基，突变后碱基，突变后样本名，突变前样本名， GRCh38，dbSNP中的名字
done

cat ./${i}/${i}.trans_site.txt | sort |uniq > ./${i}/${i}.trans_site2.txt            #去除完全相同的行
python3 getuORF1.py -i ./${i}/${i}.trans_site2.txt -o ./${i}/${i}.getuORF1.txt       #找到突变产生uORF的情况，产生的uORF相对于转录本的起始位置
python3 getuORF2.py -i ./${i}/${i}.getuORF1.txt -o ./${i}/${i}.getuORF2.txt          #找到产生uORF的终止密码子和uORF长度
Rscript defelse.R ${i}                                                               #选择TRUE的，就是有终止密码的uORF
python3 inuORF.py -i ./${i}/${i}.trans_site2.txt > ./${i}/${i}.inuORF.txt            #通过在uORF内的突变情况，分成起始密码子突变，同义突变，非同义突变
rm ${2}

done