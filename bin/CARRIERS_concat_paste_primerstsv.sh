#/bin/bash

OUTPUT_DIR=$1
OUTFILE_NAME=$2


function create_combined_primer {
    CURRENT_DIR=$(pwd)
    OUTPUT=$1
    OUT_FILE=$2
    cd $OUTPUT
    FIRST_FILE=`ls CARRIERS*.tsv | head -1`
    NEW_OUT=`cat $FIRST_FILE | cut -f1-8 > $CURRENT_DIR/$OUT_FILE`
    echo $NEW_OUT
    TEMP_FILE=$CURRENT_DIR/test.tmp
    for file in `ls CARRIERS*.tsv`;
    do 
	echo $file
	FIRST=`cat $file | cut -f8 | paste $CURRENT_DIR/$OUT_FILE - > $TEMP_FILE`
	echo $FIRST
	mv -f $TEMP_FILE $OUT_FILE
    done
}

echo $OUTPUT_DIR
echo $OUTFILE_NAME
if [ $# -eq 0 ];then
    echo "No arguments supplied"
    exit 1
else
    create_combined_primer $OUTPUT_DIR $OUTFILE_NAME
    rm -rf outfile*
    rm -rf temp*
    rm -rf qiagen_primer*
fi
