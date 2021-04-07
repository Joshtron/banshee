# Imports
from .sniffles_func import sniffles_parser
from .pbsv_func import pbsv_parser
from .nstd152_func import nstd152_parser
from .manta_func import manta_parser
from .delly_func import delly_parser
from .cnmops_func import cnmops_parser
import os

# This function recognizes the SV caller and uses the appropriate parse function that will extract
# information for bed file creation
def prepare_multiple_files(path, ts, mode, query_count):

    header_list = []
    variant_list = []
    sv_types = []

    # This loop stores header and variants in different lists
    with open(path, 'r') as file:
        for line in file:
            if not line.lstrip().startswith('#'):
                current_line = line.strip().split('\t')
                sv_type = ''.join(current_line).split('SVTYPE=')[1].split(';')[0].replace('/','')
                if not sv_type in sv_types:
                    sv_types.append(sv_type)
                variant_list.append(current_line)
            else:
                header_list.append(line.strip().split('\t'))
    file.close()

    sv_method = ''

    # The caller gets recognized by unique format features
    # VCF files get parsed by the appropriate parsing function
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

    os.mkdir('temp_bed_files/query_'+str(query_count))

    # This creates temporary bed files for the query
    if mode == 'query':

        for sv_type in sv_types:
            current_sv_file = open('temp_bed_files/query_'+str(query_count) + '/' + str(sv_type + '_query_' + str(query_count) + '.bed'), 'a+')
            for entry in bed_file:
                if sv_type == entry[3].replace('/',''):
                    current_sv_file.write('\t'.join(entry) + '\n')


    return(bed_file)