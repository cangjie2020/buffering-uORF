#bwa index ../../hg38.fa
##gatk IndexFeatureFile -I 1000G_phase1.snps.high_confidence.hg38.vcf           给vcf加idx
#gatk CreateSequenceDictionary -R hg38.fa  给fa加dic
# samtools faidx hg38.fa     给fa加索引
# gatk ReorderSam -I norna.925._.bam.marked.bam -O norna.925._.bam.marked2.bam -SD Homo_sapiens.GRCh38.dna.primary_assembly.fa 将bam顺序换成fa的

for i in `seq 25 2 51`;
do

mkdir  /data/caizixin/uorf/gatk/LIHC/SRR69399${i}
fastq=/data/caizixin/uorf/bam.sort/norrna/norna.9${i}
all_file=/data/caizixin/uorf/gatk/LIHC/SRR69399${i}
bwa_index=/data/caizixin/uorf/gatk/bwa.index/bwa.index

time bwa mem -t 24 -M -k 16 -R '@RG\tID:id\tSM:sample\tPL:Illumina' ${bwa_index} ${fastq} | samtools view -b >  ${all_file}/norna.9${i}.bam && echo "** bwa mapping done **"
time samtools sort ${all_file}/norna.9${i}.bam -o ${all_file}/norna.9${i}._.bam && echo "** BAM sort done"
rm -f ${all_file}/norna.9${i}.bam
time gatk MarkDuplicates -I ${all_file}/norna.9${i}._.bam -M ${all_file}/norna.9${i}.txt -O ${all_file}/norna.9${i}.marked.bam  && echo "** markdup done **"
time samtools index ${all_file}/norna.9${i}.marked.bam && echo "** index done **"
time gatk BaseRecalibrator -R /data/caizixin/uorf/gatk/hg38.fa -I ${all_file}/norna.9${i}.marked.bam  \
  --known-sites /data/caizixin/uorf/gatk/Homo_sapiens_assembly38.known_indels.vcf  \
  --known-sites /data/caizixin/uorf/gatk/1000G_phase1.snps.high_confidence.hg38.vcf  \
  --known-sites  /data/caizixin/uorf/gatk/Homo_sapiens_assembly38.dbsnp138.vcf  \
  -O ${all_file}/norna.9${i}.table && echo "** recalibration table done **"
time gatk ApplyBQSR --bqsr-recal-file ${all_file}/norna.9${i}.table  \
  -R /data/caizixin/uorf/gatk/hg38.fa -I ${all_file}/norna.9${i}.marked.bam -O ${all_file}/norna.9${i}.BQSR.bam && echo "** BQSR done **"

time gatk HaplotypeCaller -R /data/caizixin/uorf/gatk/hg38.fa -I ${all_file}/norna.9${i}.BQSR.bam -O ${all_file}/norna.9${i}.vcf
rm -f ${all_file}/norna.9${i}.marked.bam
time gatk VariantRecalibrator \
   -R /data/caizixin/uorf/gatk/hg38.fa \
   -V ${all_file}/norna.9${i}.vcf \
   -resource:hapmap,known=false,training=true,truth=true,prior=15.0 /data/caizixin/uorf/gatk/hapmap_3.3.hg38.vcf \
   -resource:omini,known=false,training=true,truth=false,prior=12.0 /data/caizixin/uorf/gatk/1000G_omni2.5.hg38.vcf \
   -resource:1000G,known=false,training=true,truth=false,prior=10.0 /data/caizixin/uorf/gatk/1000G_phase1.snps.high_confidence.hg38.vcf \
   -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /data/caizixin/uorf/gatk/Homo_sapiens_assembly38.dbsnp138.vcf  \
   -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
   -mode SNP \
   -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
   --tranches-file ${all_file}/norna.9${i}.snps.tranches \
   -O ${all_file}/norna.9${i}.snps.recal && echo "** SNPs VQSR1 done **"

time gatk ApplyVQSR \
    -R /data/caizixin/uorf/gatk/hg38.fa \
    -V ${all_file}/norna.9${i}.vcf \
    --truth-sensitivity-filter-level 99.0 \
    --tranches-file ${all_file}/norna.9${i}.snps.tranches \
    --recal-file ${all_file}/norna.9${i}.snps.recal \
    -mode SNP \
    -O ${all_file}/norna.9${i}.VQSR.vcf && echo "** SNPs VQSR2 done **"

  done
