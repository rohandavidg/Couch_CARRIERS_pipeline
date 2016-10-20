#!/bin/bash


#This script subsets gVCF to each
#region in the target bed file
#


JAVA=/usr/local/biotools/java/jdk1.7.0_67/bin/java
GATK=/projects/bsi/bictools/apps/alignment/GenomeAnalysisTK/3.1-1/GenomeAnalysisTK.jar
REFERENCE=/data2/bsi/reference/sequence/human/ncbi/hg19/indexed/allchr.fa
QSUB=/home/oge/ge2011.11/bin/linux-x64/qsub


usage() {
cat <<EOF
usage: $0 options
This script 
OPTIONS:
    -h   Show this message
    -g   gVCF directory path
    -b   BED FILE 
EOF
}

while getopts "hg:b:" OPTION
do
    case $OPTION in
        h) usage ; exit 1 ;;
	g) GVCF=$OPTARG ;;
	b) BED_FILE=$OPTARG ;;
        ?) usage ; exit 1 ;;
    esac
done

if [ $# -eq 0 ];then
    echo "No arguments supplied"
    usage
    exit 1
elif [ $# -lt 2 ];then
    usage
else
    true
fi

if [ -z $GVCF ];then
    echo "path to gVCF directory missing"
    exit 1
else
    true
fi

if [ -f $BED_FILE ];then
    echo "using $BED_FILE"
else
    true
fi

function gvcf_string {
    gVCF_file=`ls $GVCF/*g.vcf.gz | grep -v "tbi"`
    for gvcf in $gVCF_file;
    do
	echo "--variant" $gvcf
    done
}

function parse_bed_file {
    while read line;do
	chrom=`echo "$line" | cut -f1`
	start=`echo "$line" | cut -f2`
	stop=`echo "$line" | cut -f3`
	echo "$chrom:$start-$stop"
    done < $BED_FILE
}

    