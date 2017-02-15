#! /bin/bash

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
    -g   gVCF list
    -t   path to TOOL info   
    -o   output directory
EOF
}

while getopts "hg:t:o:" OPTION
do
    case $OPTION in
        h) usage ; exit 1 ;;
        g) gVCF_LIST=$OPTARG ;;
        t) TOOL_INFO=$OPTARG ;;
        o) OUT_DIR=$OPTARG ;;
        ?) usage ; exit 1 ;;
    esac
done

MEM="-l h_vmem=32G"
QUE="-q ngs-sec -l medp=TRUE"
#QUE="-q lg-mem"

output_dir=$OUT_DIR/variants/subsets
if [ -d $output_dir ];then
    true
else
    echo $output_dir
    mkdir $output_dir
fi

BEDTOOLS=/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bin/bedtools
export JAVA7="/usr/local/biotools/java/jdk1.7.0_03/bin"
export REF_GENOME="/data2/bsi/reference/sequence/human/ncbi/hg19/allchr.fa"
export GenotypeGVCFs_JVM="-XX:CompileThreshold=1000 -XX:MaxHeapFreeRatio=70 -XX:ReservedCodeCacheSize=256m -Xmx58g -Xms8g"
export GATK="/data5/bsi/bictools/alignment/gatk/3.4-46/"
export THREADS="1"  # Better memory effiency
export gatk_param="-R $REF_GENOME -et NO_ET -K $GATK/Hossain.Asif_mayo.edu.key"		
export output=$output_dir
export gvcfList=$gVCF_LIST
export TOOLINFO=$TOOL_INFO

function set_intervals {
    cd $OUT_DIR
    if [ ! -d intervals ]; 
    then
	TARGET_BED=`cat $TOOL_INFO | grep "ONTARGET" | cut -d"=" -f2`
	echo $TARGET_BED
	echo $TARGET_BED | xargs cat | awk '{print $1"\t"$2-200"\t"$3+200}' | $BEDTOOLS merge > target.wide.bed
	echo "Make Dir & Split Intervals"
	mkdir intervals
	cd intervals
	split -a 4 -d -l 30 ../target.wide.bed target_
	for file in *; do mv $file $file.bed; done
	cd ../
    fi
}

function run_gatk {
    cd $OUT_DIR
    for file in intervals/*; 
    do 
	f=`basename $file`;
	outputvcf="$f.vcf.gz"; 
	range="-L $OUT_DIR/$file"
	#echo "
	qsub -e $OUT_DIR/logs -o $OUT_DIR/logs $QUE -N GVCF_$f -m a -M gnanaolivu.rohandavid@mayo.edu -l h_stack=20M -pe threaded 2 -V -cwd $MEM -b y $JAVA7/java $GenotypeGVCFs_JVM -Djava.io.tmpdir=$output/temp/ -jar $GATK/GenomeAnalysisTK.jar -T GenotypeGVCFs -V $gvcfList -nt $THREADS -o $output/$outputvcf $range $gatk_param
	/projects/bsi/bictools/scripts/dnaseq/GENOME_GPS/tags/4.0.1/check_qstat.sh $TOOLINFO 8000
    done
}

if [ $# -eq 0 ];then
    echo "No arguments supplied"
    usage
    exit 1
elif [ $# -lt 3 ];then
    usage
    exit 1
else
    if [ -d $OUT_DIR ];then
	cd $OUT_DIR
	if [ ! -d logs ];then
	    mkdir logs
	    set_intervals
	    run_gatk	    
	else
	    rm -rf logs
	    mkdir logs
	    set_intervals
	    run_gatk	    	    
	fi
    else
	echo "$OUT_DIR doesnt exist"
    fi
fi

