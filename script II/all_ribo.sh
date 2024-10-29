#1预处理步骤，得到我们需要的fastq
#/data/caizixin/uorf/RSEM/kidney/SRP044937

base_file=/data/caizixin/uorf/RSEM/kidney/SRP044937

for i in `seq 25 1 35`
do
fastq-dump.2.10.8 --split-e SRR20644${i}.1;       #将SRR转成fastq
fastx_clipper  \
  -Q33 -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG \
  -l 25 -c -n –v \
  -i SRR20644${i}.1.fastq -o SRR20644${i}.fastq;
#去除接头序列，并且将太短的去除
fastq_quality_filter  \
  -Q33 -v -q 20 -p 50  \
  -i ./SRR20644${i}.fastq -o SRR20644${i}.2.fastq;
#去除质量低的序列，
bowtie2  -N 0  \
  -x  /data/caizixin/uorf/rRNA/rrna \
  -U SRR20644${i}.2.fastq \
  --un no.SRR20644${i}.fastq -S SRR;
#去除rRNA和tRNA序列
done

fastq=/data/caizixin/uorf/RSEM/kidney/SRP044937/no.SRR2064${i}.fastq
####################################################################################
#2计算转录本表达量
#/data/caizixin/uorf/RSEM/kidney/SRP044937/RSEM
mkdir  /data/caizixin/uorf/RSEM/kidney/SRP044937/RSEM
cd  /data/caizixin/uorf/RSEM/kidney/SRP044937/RSEM                              #进入操作界面

for i in `seq 24 1 25`;                                                         #使用正常细胞的RNA-seq fastq文件
do

mkdir  /data/caizixin/uorf/RSEM/kidney/SRP044937/RSEM/SRR20644${i}              #新建文件夹
mkdir  /data/caizixin/uorf/RSEM/kidney/SRP044937/RSEM/SRR20644${i}/star
mkdir  /data/caizixin/uorf/RSEM/kidney/SRP044937/RSEM/SRR20644${i}/RSEM

all_file=/data/caizixin/uorf/RSEM/kidney/SRP044937/RSEM/SRR20644${i}
star_file=/data/caizixin/uorf/RSEM/kidney/SRP044937/RSEM/SRR20644${i}/star
RSEM_file=/data/caizixin/uorf/RSEM/kidney/SRP044937/RSEM/SRR20644${i}/RSEM
fastq=/data/caizixin/uorf/RSEM/kidney/SRP044937/no.SRR20644${i}.fastq           #上一步产生的fastq文件

STAR --runThreadN 32 --genomeDir /data/caizixin/uorf/RSEM/index   \
      --readFilesIn ${fastq} --outFilterType BySJout \
      --outFilterMismatchNmax 2 --outSAMtype BAM SortedByCoordinate   \
      --quantMode TranscriptomeSAM --outFilterMultimapNmax 1 \
      --outFilterMatchNmin 16 --outFileNamePrefix    ${star_file}/SRR${i}
     # star将fastq文件比对到参考基因组上，产生转录本的表达量信息
samtools view -h  ${star_file}/SRR${i}Aligned.sortedByCoord.out.bam  > ${star_file}/SRR${i}Aligned.sortedByCoord.out.sam
#将产生的bam转成sam方便后续处理
awk 'NF == 15' ${star_file}/SRR${i}Aligned.sortedByCoord.out.sam > ${star_file}/SRR${i}Aligned15.sortedByCoord.out.sam
#将产生的bam转成sam方便后续处理
samtools sort  ${star_file}/SRR${i}Aligned.toTranscriptome.out.bam -o  ${star_file}/SRR${i}.-.bam
#将转录本的bam排序
A=`python3  /data/caizixin/uorf/RSEM/length.py   -gene_info_file ${star_file}/SRR${i}Aligned15.sortedByCoord.out.sam` ;
B=`python3  /data/caizixin/uorf/RSEM/length2.py   -gene_info_file ${star_file}/SRR${i}Aligned15.sortedByCoord.out.sam` ;
#python用于计算标准差和平均值

rsem-calculate-expression  --sort-bam-by-coordinate --alignments -p 32  \
                          --fragment-length-mean ${A} --fragment-length-sd ${B}   ${star_file}/SRR${i}.-.bam  \
                          /data/caizixin/uorf/RSEM/index/INDEX  ${RSEM_file}/SRR${i}
#RSEM --strandedness forward可以加的链特异性 线程32 单端测序的fastq
done
#循环结束接下来将数据合并
##############################################################################
#length.py
#length2.py
#####################################################

#3合并转录本表达量数据找到表达量最高的转录本，产生新的注释文件
#/data/caizixin/uorf/RSEM/kidney/SRP044937
cd ${base_file}                                                                #进入操作界面
Rscript gtf.R SRR24.isoforms.results SRR25.isoforms.results SRR26.isoforms.results
# 例子中需要的三个文件，由RSEM产生，最后产生thebestrans.GenePred文件
/data/caizixin/uorf/genePredToGtf  file  thebestrans.GenePred  new.trans.gtf
#产生新的注释文件
#######################################################
#gtf.R
#####################################################
#4 RNA-seq数据，将全部的RNA-seq数据处理
for i in `seq 424 1 429`;
do

mkdir  /data/caizixin/uorf/RSEM/kidney/SRP044937/SRR2064${i}
mkdir  /data/caizixin/uorf/RSEM/kidney/SRP044937/SRR2064${i}/star

fastq=/data/caizixin/uorf/RSEM/kidney/SRP044937/no.SRR2064${i}.fastq
GTF=/data/caizixin/uorf/RSEM/kidney/SRP044937/new.trans.gtf
all_file=/data/caizixin/uorf/RSEM/kidney/SRP044937/SRR2064${i}
star_file=/data/caizixin/uorf/RSEM/kidney/SRP044937/SRR2064${i}/star

STAR --runThreadN 32 --genomeDir /data/caizixin/uorf/RSEM/index \
      --readFilesIn ${fastq} --outFilterType BySJout \
      --outFilterMismatchNmax 2 --outSAMtype BAM SortedByCoordinate \
      --outFilterMultimapNmax 1 --outFilterMatchNmin 16 \
      --outFileNamePrefix    ${star_file}/SRR${i}
#star 比对结果
featureCounts -t exon -g gene_id -a ${GTF} -o ${star_file}/SRR${i}.gene.txt  ${star_file}/SRR${i}Aligned.sortedByCoord.out.bam
#featureCounts使用新产生的注释文件对reads注释
done

######################################
#5 Ribo_seq
for i in `seq 430 1 435`;
do

mkdir  /data/caizixin/uorf/RSEM/kidney/SRP044937/SRR2064${i}
mkdir  /data/caizixin/uorf/RSEM/kidney/SRP044937/SRR2064${i}/star
mkdir  /data/caizixin/uorf/RSEM/kidney/SRP044937/SRR2064${i}/result

all_file=/data/caizixin/uorf/RSEM/kidney/SRP044937/SRR2064${i}
star_file=/data/caizixin/uorf/RSEM/kidney/SRP044937/SRR2064${i}/star
result_file=/data/caizixin/uorf/RSEM/kidney/SRP044937/SRR2064${i}/result
read_site=/data/caizixin/uorf/RSEM/picture.py
resd_distribution=/data/caizixin/uorf/RSEM/length.R

STAR --runThreadN 32 --genomeDir /data/caizixin/uorf/RSEM/index \
      --readFilesIn ${fastq} --outFilterType BySJout --outFilterMismatchNmax 2 \
      --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM \
      --outFilterMultimapNmax 1 --outFilterMatchNmin 16 \
      --outFileNamePrefix    ${star_file}/SRR${i}
#star 跑出结果，有转录本数据
featureCounts -s 1 -R SAM -t exon -g gene_id -a ${GTF} -o ${result_file}/SRR${i}.gene.txt  ${star_file}/SRR${i}Aligned.sortedByCoord.out.bam
#将之前产生的bam文件数数，得到 每个转录本的reads数
awk 'NF == 18' ${result_file}/SRR${i}Aligned.sortedByCoord.out.bam.featureCounts.sam  >  ${result_file}/SRR${i}.18.sam
#将featurecounts产生的转录本sam文件筛选有18列的结果
python3 ${read_site} -gene_info_file ${result_file}/SRR${i}.18.sam -out ${result_file}/siteandtrans.txt
#将Sam文件的结果提取出来，对应reads和转录本
Rscript ${resd_distribution} ${result_file}/siteandtrans.txt ${result_file}/SRR${i}.length.pdf
#根据以上结果，绘制出不同长度reads分布的条形图
samtools view -h ${star_file}/SRR${i}Aligned.toTranscriptome.out.bam > ${star_file}/SRR${i}Aligned.toTranscriptome.out.sam
#将之前star产生的转录本bam文件转成sam的，有reads的位点信息，用于之后绘制图
awk 'NF == 13' ${star_file}/SRR${i}Aligned.toTranscriptome.out.sam  > ${star_file}/SRR${i}.13-.sam
#筛选我们需要的情况
done
##################
#picture.py
#length.R
##################
#6 生成三碱基周期性图，查看哪段reads积累的多

txt=/data/caizixin/uorf/RSEM/kidney/SRP044937/SRR2064${1}/result/siteandtrans.txt
sam=/data/caizixin/uorf/RSEM/kidney/SRP044937/SRR2064${1}/star/SRR${1}.13-.sam
#result=/data/caizixin/uorf/RSEM/brain/SRP031501/ SRR15625${1}/result/SRR${1}
Rscript /data/caizixin/uorf/RSEM/Psiteoffset.R  ${txt}  ${sam}  \
        /data/caizixin/uorf/RSEM/kidney/SRP044937/SRR2064${1}/result/SRR${1}  ${2}  ${3}
# siteandtrans.txt:代表着在范围内的reads，sam：代表着在哪个转录本的那个位置上的信息 ${1}是SRR号  ${2}是起始reads长度  ${3}终止reads长度

#########################
#Psiteoffset.R      大脚本
#########################
#7 导出全部uORF和CDS信息，以及整个样本的总图

file=/data/caizixin/uorf/RSEM/kidney/SRP044937/SRR2064${1}                         #哪个SRR文件下
txt=${file}/result/siteandtrans.txt
sam=${file}/star/SRR${1}.13-.sam
Rscript  /data/caizixin/uorf/RSEM/psite.R ${txt} ${sam} \
    ${file}/result/SRR${1}.CDS.count.txt \
    ${file}/result/SRR${1}.uORF.gene.count.txt \
    ${file}/result/SRR${1}.uORF.count.txt \
    ${file}/result/SRR${1}.all.phase.pdf ${2}  ${3}  ${4}  ${5}  ${6}  ${7} ${8} ${9} ${10}
#输出CDS的counts数，uORF在gene上的COUNTS，uORF在uORF上的counts，以及总的周期图，${2}是第一个起始reads长度，${3}是第一个终止reads长度，${4}是Psite偏差数
#接下来就是有三个不同的Psite可以选择范围
exit
###################################
#psite.R              大脚本
##################################33
