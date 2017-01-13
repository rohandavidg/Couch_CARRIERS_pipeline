#! /bin/bash

#######
#concats the induvidual vcfs
######

usage() {
cat <<EOF
usage: $0 options
This script concats the induvidual vcf to
create 1 variants.vcf.gz file
OPTIONS:
    -h   Show this message
    -v   path to variants dir
EOF
}

while getopts "hv:" OPTION
do
    case $OPTION in
        h) usage ; exit 1 ;;
        v) VARIANTS_DIR=$OPTARG ;;
        ?) usage ; exit 1 ;;
    esac
done


if [ $# -eq 0 ];then
    echo "No arguments supplied"
    usage
    exit 1
else 
    true
fi


OUT=$VARIANTS_DIR
SUB=$VARIANTS_DIR/subsets2
IN=$VARIANTS_DIR/subsets
TABIX=/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bin/tabix
BGZIP=/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bin/bgzip
VCFSORT=/data5/bsi/epibreast/m087494.couch/Scripts/Progams/vcftools_0.1.12a/perl/vcf-sort
echo $OUT
echo $SUB
echo $IN
if  [ -d $SUB ];then
    true
else
    mkdir $SUB
fi
### GET THE FIRST HEADER
FIRST_VCF=`ls $IN/*.bed.vcf.gz | head -1`
echo $FIRST_VCF
zcat $FIRST_VCF | grep '^#' > $SUB/SubMerge.vcf

for file in `ls $IN/target_*.bed.vcf.gz`
do
    echo $file
    zcat $file | grep -v '^#' >> $SUB/SubMerge.vcf
#    printf ". "
done
echo ""
cat $SUB/SubMerge.vcf | $VCFSORT | $BGZIP -c > $OUT/variants.vcf.gz
$TABIX -p vcf $OUT/variants.vcf.gz

