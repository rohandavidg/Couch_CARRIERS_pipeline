#/bin/bash

#####################################################
#this script is used to estimate the
#performance of the primer for each lane
#####################################################

usage() {
cat <<EOF
usage: $0 options
This script is used to estimate the performance of
the carriers primers for each lane in the run
OPTIONS:
    -h   Show this message
    -b   list of bams for analysis
    -s   sample info file
    -o   output directory
EOF
}

MAYO_EXCEL=/data5/bsi/epibreast/m087494.couch/Couch/Huge_Breast_VCF/Mayo_panel_primer.xlsx
TEST_MAYO_EXCEL=/data2/bsi/staff_analysis/m149947/couch_test/primer_coverage/test_input/TEST_MAYO_PANEL_PRIMER.xlsx
CARRIERS_TARGET_BED=/data5/bsi/epibreast/m087494.couch/Couch/Huge_Breast_VCF/CARRIERS_PANC.targets.annotated.bed
CARRIERS_PRIMER_CODE=/data2/bsi/staff_analysis/m149947/couch_code/Couch_CARRIERS_pipeline/bin/CARRIERS_primer_performance.py
PYTHON=/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bcbio/anaconda/bin/python 
QSUB=/home/oge/ge2011.11/bin/linux-x64/qsub

while getopts "hb:s:o:" OPTION
do
    case $OPTION in
        h) usage ; exit 1 ;;
        b) BAM_LIST=$OPTARG ;;
	s) SAMPLE_INFO=$OPTARG ;;
	o) OUTDIR=$OPTARG ;;
        ?) usage ; exit 1 ;;
    esac
done

LOGS_DIR=$OUTDIR/primer_analysis/logs
PRIMER_DIR=$OUTDIR/primer_analysis

function create_directory {
    if [ -d $OUTDIR ];then
	if [ -d $OUTDIR/primer_analysis/logs ];then
	    true
	else
	    mkdir -p $OUTDIR/primer_analysis/logs
	fi
	if [ -d  $OUTDIR/primer_analysis ];then
	    true
	else
	    mkdir $OUTDIR/primer_analysis
	fi
	if [ -d $OUTDIR/primer_analysis/working ];then
	    true
        else
            mkdir $OUTDIR/primer_analysis/working
        fi
    else
	echo "output directory not given"
	exit
    fi
}

function run_primer_analysis {
    for bam in `cat $BAM_LIST`;
    do
	name=`basename $bam | cut -d"." -f1`
	sleep .5
	$QSUB -N prime_$name -V -M gnanaolivu.rohandavid@mayo.edu -wd $OUTDIR/primer_analysis/working -l h_stack=20M -l h_vmem=5G -pe threaded 2 -e $LOGS_DIR -o $PRIMER_DIR -b y $PYTHON $CARRIERS_PRIMER_CODE -p $TEST_MAYO_EXCEL -b $CARRIERS_TARGET_BED -l $bam -s $SAMPLE_INFO
#	$QSUB -N prime_$name -V -M gnanaolivu.rohandavid@mayo.edu -wd $OUTDIR/primer_analysis/working -l h_stack=20M -l h_vmem=5G -pe threaded 2 -e $LOGS_DIR -o $LOGS_DIR -b y $PYTHON $CARRIERS_PRIMER_CODE -p $MAYO_EXCEL -b $CARRIERS_TARGET_BED -l $bam -s $SAMPLE_INFO
#	echo $QSUB -hold_jid -N prime_$name
    done
}


if [ $# -eq 0 ];then
    echo "No arguments supplied"
    usage
    exit 1
elif [ $# -eq 1 ];then
    usage
    exit 1
else
    create_directory
    run_primer_analysis
fi


	    