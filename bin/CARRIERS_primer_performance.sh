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
    -c   output filename
EOF
}

MAYO_EXCEL=/data5/bsi/epibreast/m087494.couch/Couch/Huge_Breast_VCF/Mayo_panel_primer.xlsx
TEST_MAYO_EXCEL=/data2/bsi/staff_analysis/m149947/couch_test/primer_coverage/test_input/TEST_MAYO_PANEL_PRIMER.xlsx
CARRIERS_TARGET_BED=/data5/bsi/epibreast/m087494.couch/Couch/Huge_Breast_VCF/CARRIERS_PANC.targets.annotated.bed
CARRIERS_PRIMER_CODE=/data2/bsi/staff_analysis/m149947/couch_code/Couch_CARRIERS_pipeline/bin/CARRIERS_primer_performance.py
PYTHON=/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bcbio/anaconda/bin/python 
QSUB=/home/oge/ge2011.11/bin/linux-x64/qsub
CARRIERS_CONCAT_PRIMER=/data2/bsi/staff_analysis/m149947/couch_code/Couch_CARRIERS_pipeline/bin/CARRIERS_concat_paste_primerstsv.sh

while getopts "hb:s:o:c:" OPTION
do
    case $OPTION in
        h) usage ; exit 1 ;;
        b) BAM_LIST=$OPTARG ;;
	s) SAMPLE_INFO=$OPTARG ;;
	o) OUTDIR=$OPTARG ;;
	c) OUTPUT_FILENAME=$OPTARG ;;
        ?) usage ; exit 1 ;;
    esac
done

LOGS_DIR=$OUTDIR/primer_analysis/logs
PRIMER_DIR=$OUTDIR/primer_analysis
echo $SAMPLE_INFO

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

a=0
JobIds=""
function run_primer_analysis {
    for bam in `cat $BAM_LIST`;
    do
	a=$((a+1))
#	echo $a
        if [ $a -eq 1 ]; then
            JobIds="primer_run_1"
        else
            JobIds=$JobIds','primer_run_$a
        fi
	name=`basename $bam | cut -d"." -f1`
#	sleep .5
#	$QSUB -N primer_run_$a -V -M gnanaolivu.rohandavid@mayo.edu -wd $OUTDIR/primer_analysis/working -l h_stack=20M -l h_vmem=5G -pe threaded 2 -e $LOGS_DIR -o $PRIMER_DIR -b y $PYTHON $CARRIERS_PRIMER_CODE -p $TEST_MAYO_EXCEL -b $CARRIERS_TARGET_BED -l $bam -s $SAMPLE_INFO
	$QSUB -N primer_run_$a -V -M gnanaolivu.rohandavid@mayo.edu -wd $OUTDIR/primer_analysis/working -l h_stack=20M -l h_vmem=5G -pe threaded 2 -e $LOGS_DIR -o $LOGS_DIR -b y $PYTHON $CARRIERS_PRIMER_CODE -p $MAYO_EXCEL -b $CARRIERS_TARGET_BED -l $bam -s $SAMPLE_INFO
    done
}

if [ $# -eq 0 ];then
    echo "No arguments supplied"
    usage
    exit 1
elif [ $# -lt 5 ];then
    usage
    exit 1
else
    create_directory
    run_primer_analysis
    STUNT=`echo $JobIds | rev | cut -d"," -f1-140 | rev`
    sleep 10m
    $QSUB -N Primer_merge$RANDOM  -hold_jid $STUNT -V -M gnanaolivu.rohandavid@mayo.edu -wd $OUTDIR/primer_analysis/working -l h_stack=20M -l h_vmem=5G -pe threaded 2 -e $LOGS_DIR -o $PRIMER_DIR -b y /bin/bash $CARRIERS_CONCAT_PRIMER $OUTDIR/primer_analysis/working $OUTPUT_FILENAME
fi
