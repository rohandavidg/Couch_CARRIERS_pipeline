#!/data5/bsi/epibreast/m087494.couch/Scripts/Progams/bcbio/anaconda/bin/python


"""
This script annotates the vcf with 
exon numbers
"""


import pybedtools
import os
import argparse
import sys
import csv
import gzip
import multiprocessing
import shlex

transcript_mapping='/data5/bsi/epibreast/m087494.couch/Reference_Data/TranscriptMapping/NormalizedTranscripts.tsv'


def main():
    args = parse_args()
    run(args.gtf_file, args.CNV_txt_list, args.output_dir)


def run(gtf_file, CNV_txt_list, output_dir):
    subset_gtf_file = output_dir + '/' + 'gtf_exon_number.bed'
    gene_desired_transcript_dict = map_gene_transcript(transcript_mapping)
    if os.path.isfile(subset_gtf_file):
        pass
    else:
        subset_gtf = parse_gtf(gtf_file, gene_desired_transcript_dict, subset_gtf_file)
    create_cnv_seg_intersected = parse_cnv_txt(CNV_txt_list, subset_gtf_file, output_dir)
#    print create_cnv_seg_intersected


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-c', dest='CNV_txt_list', type=argparse.FileType('r'),
                        help="path to the ggps bam directory", required=True)
    parser.add_argument('-g', dest='gtf_file', 
                        default='/data2/bsi/staff_analysis/m112876/MAPRSeq/MAPRSeq_Alt_Species/Human/Homo_sapiens.GRCh37.75.gtf.gz',
                        help='path to gtf to use')
    parser.add_argument('-o', dest='output_dir', required=True,
                        help='path to create tmp and final output')
    args = parser.parse_args()
    return args


def map_gene_transcript(TRANSCRIPT_MAPPING):
    gene_transcript_dict = {}
    with open(TRANSCRIPT_MAPPING) as fin:
        for raw_line in fin:
            value = raw_line.strip().split('\t')
            gene = value[0]
            transcript = value[2]
            gene_transcript_dict[gene] = transcript
    return gene_transcript_dict


def parse_gtf(gtf_file, gene_transcript_dict, outfile):
    with gzip.open(gtf_file) as csvfile, open(outfile, 'w') as fout:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            try:
                if row[2] == 'CDS':
                    chrom = 'chr' + row[0]
                    start = row[3]
                    stop = row[4]
                    transcript_id = row[8].split(';')[1]
                    exon_number = row[8].split(';')[2]
                    exon = shlex.split(exon_number.split(" ")[2])[0]
                    gene_name = row[8].split(';')[3]
                    gene = shlex.split(gene_name.split(" ")[2])[0]
                    transcript = shlex.split(transcript_id.split(" ")[2])[0]
                    try:
                        if transcript == gene_transcript_dict[gene]:
                            out = (chrom, start, stop, gene, transcript, exon)
                            fout.write('\t'.join(str(i) for i in out) + '\n')
                        else:
                            pass
                    except KeyError:
                        pass
                else:
                    pass
            except IndexError:
                pass
                

def parse_cnv_txt(CNV_txt_list, gtf_file, output_dir):
    for cnv_file in CNV_txt_list:
        get_seg_bed, sample_name = get_seg_bed_file(cnv_file)
        CNV_txt_bed = output_dir + '/tmp_out.bed'
        out_headers = convert_txt_to_bed(cnv_file, CNV_txt_bed)
        intersect_seg_CNV_bed = output_dir + '/CNV_seg_intersected.bed'
        intersected_CNV_txt = intersect_txt_seg(CNV_txt_bed, get_seg_bed, intersect_seg_CNV_bed)
        intersected_CNV_seq_gtf = intersect_cnv_txt_with_gtf(intersect_seg_CNV_bed, gtf_file, sample_name, out_headers, output_dir) 
        write_out_cnv_txt_exon = get_req_columns(intersected_CNV_seq_gtf, sample_name, output_dir)
        remove_tmp = cleanup_dir(CNV_txt_bed)
        remove_tmp = cleanup_dir(intersect_seg_CNV_bed)
        remove_tm = cleanup_dir(intersected_CNV_seq_gtf)

def get_seg_bed_file(CNV_file):
    sample_file = os.path.basename(CNV_file)
    base_dir = os.path.dirname(CNV_file)
    sample = sample_file.split('(')[0]
    for root, dirs, filename in os.walk(base_dir):
        for name in filename:            
            if name.startswith(sample) and name.endswith('seg_.bed'):
                path = root + "/" + name
                return path, sample


def convert_txt_to_bed(CNV_file, out_file):
    header_list = []
    with open(CNV_file.strip()) as fin, open(out_file, 'w') as fout:
        header = fin.readline().strip().split('\t')
        header_list += header
        for raw_line in fin:
            fout.write(raw_line.strip() + '\n')
    return header_list


def intersect_txt_seg(CNV_txt_bed, seg_bed, out_file):
    with open(out_file, 'w') as fout:
        a = pybedtools.BedTool(CNV_txt_bed)
        b = pybedtools.BedTool(seg_bed)
        CNV_intersected = a.intersect(b)
        fout.write(''.join(str(i) for i in CNV_intersected) + '\n')


def intersect_cnv_txt_with_gtf(cnv_bed, subset_gtf_bed, sample, header, out_dir):
    new_header = header + ['', '', '', '', 'Transcript', 'Exon_number']
    outfile = out_dir + "/" + sample + '_germline.transcript.exon.txt'
    with open(outfile, 'w') as fout:
        fout.write('\t'.join(str(i) for i in new_header) + '\n')
        a = pybedtools.BedTool(cnv_bed)
        b = pybedtools.BedTool(subset_gtf_bed)
        germline_intersected = a.intersect(b, wb=True)
        get_snps = a.intersect(b, v=True)
        string_germline_intersected = ''.join(str(i) for i in germline_intersected)
        fout.write(string_germline_intersected + '\n')
        fout.write(''.join(str(i) for i in get_snps) + '\n')
    return outfile


def get_req_columns(outfile, sample, output_dir):
    new_outfile = output_dir + '/' + sample+ '_CNV_germline.transcript.exon.txt'
    with open(outfile, "r") as fin, open(new_outfile, 'w') as fout:
        data = fin.readlines()
        for raw_lines in data:
            line = raw_lines.strip()
            value = line.split('\t')
            value2 = filter(None, value)
            if value2:
                result = value2[:9] + value2[-2:]
                fout.write('\t'.join(str(i) for i in result) + '\n')


def cleanup_dir(filename):
    assert os.path.isfile(filename)
    os.remove(filename)

if __name__ == '__main__':
    main()
