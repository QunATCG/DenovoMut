#!/bin/bash
# Exome-seq data processing pipeline
# 2020-11-12
# LiQun
# version 1.0
# mode: step by step

#################### before #####################
# ref = hg19									#
# sequence layout = paire						#
# conda install -c bioconda bwa, fastp, gatk4...#
#################################################

#################### input ######################
#### default ####
LAYOUT=paired
DIRNOW=`pwd`
VCFOUT=${DIRNOW}/vcfout
TRIM_DIR=${DIRNOW}/trim
TRIM_REPORT_DIR=${DIRNOW}/trim_report
THREADS=8
GENOME=/data/public/references/hg19_BWAIndex/genome.fa
BWAINDEX=/data/public/references/hg19_BWAIndex/
####  other  ####
STEPS=$1
SAMPLE_DIR=$2
#### help ####
if [ -z $1 ] || [ -z $2 ] ;then
		echo -e "\033[32m{Useage: $0 STEPS SAMPLE_DIR}\033[0m"
		echo -e "\033[32m  STEPS: example:1,2,3,4,5,6,7\033[0m"
		#echo -e "\033[32m  STEP0: BWA index\033[0m"
		echo -e "\033[32m  STEP1: get filenames\033[0m"
		echo -e "\033[32m  STEP2: RawData quality control(QC)\033[0m"
		echo -e "\033[32m  STEP3: Start Mapping & Convert to bam and Sort\033[0m"
		#echo -e "\033[32m  STEP4: Convert to bam and Sort\033[0m"
		echo -e "\033[32m  STEP4: Mark Duplicate\033[0m"
		echo -e "\033[32m  STEP5: RBQS\033[0m"
		echo -e "\033[32m  STEP6: Varient calling\033[0m"
		echo -e "\033[32m  STEP7: Mutation\033[0m"
		echo -e "\033[32m  SAMPLE_DIR: raw samples(paired & *fastq.gz) path \033[0m"
		#echo -e "\033[32m  GENOME_TAG: hg19/hg38 \033[0m"
		exit
fi

#### KONWN ANNOTATION INFORMATION ####
DBSNP=/data/public/references/gatkhg19/dbsnp_138.hg19.vcf
DBSNP1000G=/data/public/references/gatkhg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf
DBINDEL1000G=/data/public/references/gatkhg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
OMNI=/data/public/references/gatkhg19/1000G_omni2.5.hg19.sites.vcf
HAPMAP=/data/public/references/gatkhg19/hapmap_3.3.hg19.sites.vcf
#### START ####
echo -e "\033[32m======Your Analysis steps: $STEPS ======\033[0m"
echo -e "\033[32m======Sequence layout is $LAYOUT ======\033[0m"
echo -e "\033[32m======Sample DIR is $SAMPLE_DIR ======\033[0m"
echo -e "\033[32m======Trim out DIR is $TRIM_DIR ======\033[0m"
echo -e "\033[32m======WES analysis is running ======\033[0m"

#### STEP1 Get filenames
if [[ $STEPS =~ '1' ]];then
	FILENAMES=$(ls ${SAMPLE_DIR} > filenames.txt)
	R1_FILENAMES=$(cat filenames.txt | grep -E "R1" > R1.txt)
	#R2_FILENAMES=$(cat filenames.txt | grep -E "R2" > R2.txt)
fi

#### STEP2 RawData quality control(QC) ####
if [[ $STEPS =~ '2' ]];then
	echo -e "\033[32m======Start RawData quality control(QC)======\033[0m"
	if [ ! -d $TRIM_DIR ];then
		mkdir -p $TRIM_DIR
	fi
	if [ ! -d $TRIM_REPORT_DIR ];then
		mkdir -p $TRIM_REPORT_DIR
	fi
	while read line
	do
		TAG_R=${line##*R1}
		TAG_L=${line%%R1*}
		R1_FASTQ=$TAG_L"R1"$TAG_R
		R2_FASTQ=$TAG_L"R2"$TAG_R
		R1_FASTQNAME=${R1_FASTQ%%.*}
		R2_FASTQNAME=${R2_FASTQ%%.*}
		#echo $TAG_L
		#echo $TAG_R
		echo -e "\033[32m======QC......$line======\033[0m"
		fastp -i $SAMPLE_DIR/$R1_FASTQ -o $TRIM_DIR/${R1_FASTQNAME}_clean.fastq.gz -I $SAMPLE_DIR/$R2_FASTQ -O $TRIM_DIR/${R2_FASTQNAME}_clean.fastq.gz -h $TRIM_REPORT_DIR/${TAG_L}fastp.html -j $TRIM_REPORT_DIR/${TAG_L}fastp.json
	done < R1.txt
	echo -e "\033[32m======QC analysis DONE ======\033[0m"
fi

#### STEP3 Mapping and Sort ####
if [[ $STEPS =~ '3' ]];then
	echo -e "\033[32m======Start Mapping ======\033[0m"
	while read line
	do
		TAG_R=${line##*R1}
		TAG_L=${line%%R1*}
		R1_FASTQ=$TAG_L"R1"$TAG_R
		R2_FASTQ=$TAG_L"R2"$TAG_R
		R1_FASTQNAME=${R1_FASTQ%%.*}
		R2_FASTQNAME=${R2_FASTQ%%.*}
		SAMPLE=${line%%_R1*}
		echo -e "\033[32m======BWA......$SAMPLE ======\033[0m"
		bwa mem -t $THREADS -R "@RG\tID:${SAMPLE}\tPL:illumina\tLB:WES\tSM:${SAMPLE}" $GENOME  $TRIM_DIR/${R1_FASTQNAME}_clean.fastq.gz $TRIM_DIR/${R2_FASTQNAME}_clean.fastq.gz > $TRIM_DIR/${SAMPLE}.sam
		echo -e "\033[32m======SORTSAM ======\033[0m"
		gatk --java-options "-Xmx20G" SortSam -I $TRIM_DIR/${SAMPLE}.sam -O $TRIM_DIR/${SAMPLE}.sort.bam -SO coordinate --CREATE_INDEX true
		#rm $TRIM_DIR/${SAMPLE}.sam
		#samtools flagstat -@ $THREADS $TRIM_DIR/${SAMPLE}.sort.bam
		rm $TRIM_DIR/${SAMPLE}.sam
	done < R1.txt
	echo -e "\033[32m======Mapping DONE ======\033[0m"
	#rm $TRIM_DIR/*.sam
fi
# STEP4  Mark Duplicate
if [[ $STEPS =~ '4' ]];then
		echo -e "\033[32m======Mark Duplicate======\033[0m"
		while read line
		do
			TAG_R=${line##*R1}
			TAG_L=${line%%R1*}
			R1_FASTQ=$TAG_L"R1"$TAG_R
			R2_FASTQ=$TAG_L"R2"$TAG_R
			R1_FASTQNAME=${R1_FASTQ%%.*}
			R2_FASTQNAME=${R2_FASTQ%%.*}
			SAMPLE=${line%%_R1*}
			echo -e "\033[32m======Mark Duplicate......$SAMPLE ======\033[0m"
			gatk --java-options "-Xmx20G" MarkDuplicates -I $TRIM_DIR/${SAMPLE}.sort.bam -O $TRIM_DIR/${SAMPLE}.sorted.marked.bam -M $TRIM_DIR/${SAMPLE}.metrics
		done < R1.txt
		echo -e "\033[32m======Mark Duplicate DONE ======\033[0m"
fi

# STEP5 RBSQ
if [[ $STEPS =~ '5' ]];then
		echo -e "\033[32m======RBQS....======\033[0m"
		while read line
		do
			TAG_R=${line##*R1}
			TAG_L=${line%%R1*}
			R1_FASTQ=$TAG_L"R1"$TAG_R
			R2_FASTQ=$TAG_L"R2"$TAG_R
			R1_FASTQNAME=${R1_FASTQ%%.*}
			R2_FASTQNAME=${R2_FASTQ%%.*}
			SAMPLE=${line%%_R1*}
			echo -e "\033[32m======RBQS......$SAMPLE ======\033[0m"
			gatk --java-options "-Xmx20G" BaseRecalibrator -R $GENOME -I $TRIM_DIR/${SAMPLE}.sorted.marked.bam -O $TRIM_DIR/${SAMPLE}.recal_data.table --known-sites $DBSNP --known-sites $DBSNP1000G --known-sites $DBINDEL1000G
			gatk --java-options "-Xmx20G" ApplyBQSR -R $GENOME -I $TRIM_DIR/${SAMPLE}.sorted.marked.bam --bqsr-recal-file ${TRIM_DIR}/${SAMPLE}.recal_data.table -O $TRIM_DIR/${SAMPLE}.sorted.marked.BQSR.bam
		done < R1.txt
		echo -e "\033[32m======RBSQ DONE ======\033[0m"
fi

# STEP6 Variant calling
if [[ $STEPS =~ '6' ]];then
		echo -e "\033[32m======Varient calling======\033[0m"
	if [ ! -d $VCFOUT ];then
		mkdir -p $VCFOUT
	fi
		while read line
		do
			TAG_R=${line##*R1}
			TAG_L=${line%%R1*}
			R1_FASTQ=$TAG_L"R1"$TAG_R
			R2_FASTQ=$TAG_L"R2"$TAG_R
			R1_FASTQNAME=${R1_FASTQ%%.*}
			R2_FASTQNAME=${R2_FASTQ%%.*}
			SAMPLE=${line%%_R1*}
			echo -e "\033[32m======Variant calling......$SAMPLE ======\033[0m"
			gatk --java-options "-Xmx20G" HaplotypeCaller -R $GENOME -I $TRIM_DIR/${SAMPLE}.sorted.marked.BQSR.bam -ERC GVCF -O $VCFOUT/${SAMPLE}.g.vcf
		done < R1.txt
		echo -e "\033[32m======Variant DONE ======\033[0m"
fi

# STEP7 mutation
if [[ $STEPS =~ '7' ]];then
		echo -e "\033[32m======Denovo calling======\033[0m"
	#获取g.vcf文件
		ls $VCFOUT/ | grep -E "*.g.vcf$" | while read line; do echo "vcfout/"${line}; done > input.list
		while read line
		do
			TAG_R=${line##*R1}
			TAG_L=${line%%R1*}
			R1_FASTQ=$TAG_L"R1"$TAG_R
			R2_FASTQ=$TAG_L"R2"$TAG_R
			R1_FASTQNAME=${R1_FASTQ%%.*}
			R2_FASTQNAME=${R2_FASTQ%%.*}
			SAMPLE=${line%%_R1*}
			echo -e "\033[32m======Varient calling......$SAMPLE ======\033[0m"
			gatk --java-options "-Xmx20G" CombineGVCFs -R $GENOME --variant input.list -O combined.g.vcf
			gatk --java-options "-Xmx20G" GenotypeGVCFs -R $GENOME -V combined.g.vcf  -O genetype.output.vcf
			gatk VariantRecalibrator -R $GENOME -V genetype.output.vcf --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP --resource:omni,known=false,training=true,truth=false,prior=12.0 $OMNI --resource:1000G,known=false,training=true,truth=false,prior=10.0 $DBSNP1000G --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP  -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP  -O output.snp.recal --tranches-file output.snp.tranches --rscript-file output.snp.plots.R
			gatk ApplyVQSR -R $GENOME  -V genetype.output.vcf -O output.snp.vcf  --tranches-file output.snp.tranches --recal-file output.snp.recal -mode SNP
			gatk VariantRecalibrator -R $GENOME -V genetype.output.vcf --max-gaussians 4 --resource:mills,known=false,training=true,truth=true,prior=12.0 $DBINDEL1000G --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -mode INDEL  -O output.indel.recal --tranches-file output.indel.tranches --rscript-file output.indel.plots.R
			gatk ApplyVQSR -R $GENOME  -V genetype.output.vcf -O output.indel.vcf --tranches-file output.indel.tranches --recal-file output.indel.recal -mode INDEL
		done < R1.txt
		echo -e "\033[32m======Denovo DONE ======\033[0m"
fi
