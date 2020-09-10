#!/bin/bash

################################################################################
############################ Estimate nucleotide diversity #####################
################################## JU775 MY16 N2 ###############################
################################################################################

###load modules
#module load samtools bcftools easybuild  icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132 VCFtools/0.1.15-Perl-5.24.1 bedtools


###program versions:
	#samtools Version: 1.5 (using htslib 1.5)
	#bcftools Version: 1.5 (using htslib 1.5)
	#vcftools VCFtools (0.1.15)
	#bedtools v2.25.0


################################################################################
###download the CeNDR data
wget --no-check-certificate https://elegansvariation.org.s3.amazonaws.com/bam/JU775.bam
wget --no-check-certificate https://elegansvariation.org.s3.amazonaws.com/bam/JU775.bam.bai
wget --no-check-certificate https://elegansvariation.org.s3.amazonaws.com/bam/MY16.bam
wget --no-check-certificate https://elegansvariation.org.s3.amazonaws.com/bam/MY16.bam.bai
wget --no-check-certificate https://elegansvariation.org.s3.amazonaws.com/bam/N2.bam
wget --no-check-certificate https://elegansvariation.org.s3.amazonaws.com/bam/N2.bam.bai


################################################################################
###variant calling

ref="/projects/phillipslab/ateterina/CeNDR/ref_245/c_elegans.PRJNA13758.WS245.genomic.fa"

mkdir -p N2_MY16_JU775
cd N2_MY16_JU775


#cat bams_CITP.file
#N2.bam
#MY16.bam
#JU775.bam


###call variants
samtools mpileup -q 30 -Q 15 -uf $ref -s $(cat bams_CITP.file) | bcftools call -mv > N2_MY16_JU775.raw.vcf

###some filtering
vcftools --vcf N2_MY16_JU775.raw.vcf --out N2_MY16_JU775.filt.Q20 --minQ 20  --recode --recode-INFO-all ;

################################################################################
###collect all sites that have too low/high coverage al least in 1 strain
for file in N2.bam MY16.bam JU775.bam; do

	samtools depth -a $file	| awk '{ if ($3 < 10 || $3 > 500) { print $1 "\t" $2-1 "\t" $2 }}' - >> N2_MY16_JU775.DP10-500.mask.bed

done

###create a mask

samtools faidx $ref
bedtools sort -faidx ${ref}.fai -i N2_MY16_JU775.DP10-500.mask.bed | bedtools merge -i - > N2_MY16_JU775.DP10-500.mask.s.bed


################################################################################

###remove indels and non biallelic sites
awk '{ if($1 ~ /#/ || $5 ~/^[ATGC]$/){print;} }' N2_MY16_JU775.filt.Q20.recode.vcf > N2_MY16_JU775.filt.Q20.ok.recode.vcf

###add a header
grep "#" N2_MY16_JU775.filt.Q20.ok.recode.vcf > N2_MY16_JU775.filt.Q20.DP10-500.ok2.recode.vcf

###remove some sites
bedtools subtract -a N2_MY16_JU775.filt.Q20.ok.recode.vcf -b N2_MY16_JU775.DP10-500.mask.s.bed >>N2_MY16_JU775.filt.Q20.DP10-500.ok2.recode.vcf


################################################################################

###estimate nucleotide diversity for each polymorphic site
vcftools --vcf N2_MY16_JU775.filt.Q20.DP10-500.ok2.recode.vcf --site-pi --out N2_MY16_JU775.filt.Q20.DP10-500.pi


################################################################################

###the genome size
genome=$(awk '{sum+=$2} END {print sum}' ${ref}.fai)
#100286401

###the sun of pi
totalpi=$(cat N2_MY16_JU775.filt.Q20.DP10-500.pi.sites.pi | awk '{sum+=$3} END {print sum}' -)
#176206

### the number of position with too low/high coverage in N2, MY16, JU775
masked=$(cat N2_MY16_JU775.DP10-500.mask.s.bed | awk '{sum+=$3-$2} END {print sum}')
#5336625

###the number of removed variants
removed=$(($(wc -l < N2_MY16_JU775.filt.Q20.recode.vcf)- $(wc -l < N2_MY16_JU775.filt.Q20.DP10-500.ok2.recode.vcf)))
#114888

echo $totalpi $genome $masked $removed
#176206 100286401 5336625 114888

###nucleotide diversity
echo $(($totalpi/($genome-$masked-$removed)))
#it's 0.001858029294029
