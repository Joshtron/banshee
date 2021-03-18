from .sniffles_func import sniffles_parser
from .pbsv_func import pbsv_parser
from .nstd152_func import nstd152_parser
from .manta_func import manta_parser
from .delly_func import delly_parser
from .cnmops_func import cnmops_parser
import os

def prepare_files(path, ts, mode):

    header_list = []
    variant_list = []

    with open(path, 'r') as file:
        for line in file:
            if not line.lstrip().startswith('#'):
                variant_list.append(line.strip().split('\t'))
            else:
                header_list.append(line.strip().split('\t'))
    file.close()

    sv_method = ''

    if len(variant_list[0]) == 10:
        if 'Sniffles' in variant_list[0][7]:
            sv_method = 'sniffles'
            bed_file = sniffles_parser(path)
        elif 'pbsv' in variant_list[0][2]:
            sv_method = 'pbsv'
            bed_file = pbsv_parser(path)
        elif 'Manta' in variant_list[0][2]:
            sv_method = 'manta'
            bed_file = manta_parser(path)
        elif 'DELLY' in variant_list[0][7]:
            sv_method = 'delly'
            bed_file = delly_parser(path)
        #elif 'source=svaba' in header_list[2][0]:
        #    sv_method = 'svaba'
        #    bed_file = svaba_parser(path)
    else:
        if 'cnmops' in header_list[0][0]:
            sv_method = 'cnmops'
            bed_file = cnmops_parser(path)
        elif 'nssv' in variant_list[0][2]:
            sv_method = 'nstd152'
            bed_file = nstd152_parser(path, ts)

    if mode == 'query':

        del_file_query = open('temp_bed_files/del_query.bed', 'a+')
        ins_file_query = open('temp_bed_files/ins_query.bed', 'a+')
        dup_file_query = open('temp_bed_files/dup_query.bed', 'a+')
        inv_file_query = open('temp_bed_files/inv_query.bed', 'a+')

        for entry in bed_file:
            if 'DEL' in entry[3]:
                del_file_query.write('\t'.join(entry) + '\n')
            if 'INS' in entry[3]:
                ins_file_query.write('\t'.join(entry) + '\n')
            if 'DUP' in entry[3]:
                dup_file_query.write('\t'.join(entry) + '\n')
            if 'INV' in entry[3]:
                inv_file_query.write('\t'.join(entry) + '\n')

    #os.system('bedtools merge -d 50 -c 4,4,5 -o count,collapse,collapse -i temp_bed_files/del_query.bed > temp_bed_files/merged_del_query.bed')

    if mode == 'truth':

        del_file_truth = open('temp_bed_files/del_truth.bed', 'a+')
        ins_file_truth = open('temp_bed_files/ins_truth.bed', 'a+')
        dup_file_truth = open('temp_bed_files/dup_truth.bed', 'a+')
        inv_file_truth = open('temp_bed_files/inv_truth.bed', 'a+')

        for entry in bed_file:
            if 'DEL' in entry[3]:
                del_file_truth.write('\t'.join(entry) + '\n')
            if 'INS' in entry[3]:
                ins_file_truth.write('\t'.join(entry) + '\n')
            if 'DUP' in entry[3]:
                dup_file_truth.write('\t'.join(entry) + '\n')
            if 'INV' in entry[3]:
                inv_file_truth.write('\t'.join(entry) + '\n')


    return(bed_file)