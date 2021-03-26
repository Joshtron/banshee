import argparse

# The module above enables easy command line usability
from typing import List

parser = argparse.ArgumentParser(description='Simplify your .vcf file  --> python vcf_parser.py [tool] [path to file]')

# This adds new arguments that have to be entered to make the script work. I just thought about adding different tools
# over time and store each tool in a function. Also the path to the vcf file has to be specified
parser.add_argument("--stats",
                    help="shows stats for SV-types", action="store_true")
parser.add_argument('--format', default='bed',
                    help='please specify the output format (.bed, .csv)')
parser.add_argument("--truthset", default='None',
                    help='please specify truthset')
parser.add_argument('--svlen',
                    help='please specify your structural variant length')
parser.add_argument('path', type=str,
                    help='please specify the path to the file')
args = parser.parse_args()

import string
import svaba_script


# This is necessary to filter out all numbers from a string, e.g. the length of a variant coming from an info field
class Del:
    def __init__(self, keep=string.digits):
        self.comp = dict((ord(c), c) for c in keep)

    def __getitem__(self, k):
        return self.comp.get(k)


DD = Del()

software_list = ['sniffles', 'pbsv', 'manta', 'delly', 'svaba', 'cnmops', 'nstd152 (Yorubian Truthset)']

complete_list = []

number_del = 0
number_ins = 0
number_dup = 0
number_inv = 0
number_nes = 0

header_list = []
variant_list = []

# This opens the file and stores every line that does not strart with '#' in a list. With that I can skip the header
with open(str(args.path), 'r') as file:
    for line in file:
        if not line.lstrip().startswith('#'):
            variant_list.append(line.strip().split('\t'))
        else:
            header_list.append(line.strip().split('\t'))
file.close()

# This gets the calling method from the info field

sv_method = ''

if len(variant_list[0]) == 10:
    if 'Sniffles' in variant_list[0][7]:
        sv_method = 'sniffles'
    elif 'pbsv' in variant_list[0][2]:
        sv_method = 'pbsv'
    elif 'Manta' in variant_list[0][2]:
        sv_method = 'manta'
    elif 'DELLY' in variant_list[0][7]:
        sv_method = 'delly'
    elif 'source=svaba' in header_list[2][0]:
        sv_method = 'svaba'
else:
    if 'cnmops' in header_list[0][0]:
        sv_method = 'cnmops'
    elif 'nssv' in variant_list[0][2]:
        sv_method = 'nstd152'

if sv_method == 'sniffles':

    variant_list = []

    # This opens the file and stores every line that does not strart with '#' in a list. With that I can skip the header
    with open(str(args.path), 'r') as file:
        for line in file:
            if not line.lstrip().startswith('#'):
                variant_list.append(line.strip().split('\t'))
    file.close()


    # This filters for the end position of a variant as the position of that in the .vcf file is 'always' the same when
    # coming from sniffles
    # If the SV type is 'BND' the info field SVLENGTH gets filtered for the number and added to the start of the variant
    def get_end_position(variant_list_entry):
        if 'SVTYPE=BND' in variant_list_entry[7]:
            info_field = variant_list_entry[7].split(';')
            end_of_variant = str(int(variant_list_entry[1]) + int(info_field[8].translate(DD)))
        else:
            end_position = variant_list_entry[7].split(';')[3]
            end_of_variant = end_position[4:]
        return end_of_variant


    # Same for the sv-type.
    #
    def get_sv_type(variant_list_entry):
        if not 'SVTYPE=BND' in variant_list_entry[7]:
            end_position = variant_list_entry[7].split(';')[8]
            sv_type = end_position[7:]
            return sv_type


    number_of_imprecise_EBV = 0

    for entry in variant_list:
        info_field = entry[7].split(';')
        # If-command filters out the decoys as I tought they might cause problems. Delete if they are needed
        # Filters also the strange variants at the end that have a altered INFO field. Functions above do not work or
        # need to be altered. We should first talk what they are exactly as I have never heard that sv type
        # Also filters out EBV variants
        if not 'IMPRECISE' in info_field[0] and not 'SVTYPE=BND' in entry[7] and not '_' in entry[0] and not 'EBV' in \
                                                                                                             entry[0]:
            chr = entry[0]
            start = entry[1]
            end = get_end_position(entry)
            sv_type = get_sv_type(entry)
            if sv_type == 'INS' or sv_type == 'DUP':
                sv_len = int(end) - int(start)
                if sv_len < 51:
                    to_add = int((51 - sv_len) / 2)
                    sublist = []
                    if (sv_len % 2) == 0:
                        sublist = [chr, str(int(start) - to_add), str(int(end) + to_add + 1), 'INS']
                    else:
                        sublist = [chr, str(int(start) - to_add), str(int(end) + to_add), 'INS']
                    sublist.append('_'.join(sublist))
                    complete_list.append(sublist)
                else:
                    sublist = [chr, start, end, 'INS']
                    sublist.append('_'.join(sublist))
                    complete_list.append(sublist)
                number_ins += 1
            elif sv_type == 'INV':
                sublist_start = [chr, str(int(start) - 25), str(int(start) + 25), sv_type]
                sublist_end = [chr, str(int(end) - 25), str(int(end) + 25), sv_type]
                sublist_start.append('_'.join(sublist_start))
                sublist_end.append('_'.join(sublist_end))
                complete_list.extend((sublist_start, sublist_end))
                number_inv += 1
            elif sv_type == 'DEL':
                sublist = [chr, start, end, sv_type]
                sublist.append('_'.join(sublist))
                complete_list.append(sublist)
                number_del += 1
            else:
                sublist = [chr, start, end, sv_type]
                sublist.append('_'.join(sublist))
                complete_list.append(sublist)
                number_nes += 1
        else:
            number_of_imprecise_EBV += 1









elif sv_method == 'nstd152':

    ts = str(args.truthset)
    if ts == 'None':
        print('\nPlease specify a valid truthset\n',
              'Yorubian Trio: NA19238, NA19239, NA19240\n',
              'Han chinese Trio: HG00512, HG00513, HG00514\n',
              'Ashkenazim Trio: HG00731, HG00732, HG00733')
    else:
        variant_list = []

        # This opens the file and stores every line that does not strart with '#' in a list. With that I can skip the header
        with open(str(args.path), 'r') as file:
            for line in file:
                if not line.lstrip().startswith('#'):
                    variant_list.append(line.strip().split('\t'))
        file.close()


        def get_end_position(variant_list_entry):
            info_field = variant_list_entry[7].split(';')
            for item in info_field:
                if 'END' in item and not 'CIEND' in item:
                    end_of_variant = item.translate(DD)
            return end_of_variant


        def get_sv_type(variant_list_entry):
            info_field = variant_list_entry[7].split(';')
            for item in info_field:
                if 'SVTYPE' in item:
                    sv_type = item[7:]
            return sv_type


        complete_list = []

        for entry in variant_list:
            info_field = entry[7].split(';')
            info_field_str = ''.join(info_field)
            chr = 'chr' + entry[0]
            start = entry[1]
            end = get_end_position(entry)
            sv_type = get_sv_type(entry)
            if ts in info_field_str and sv_type == 'INV':
                sublist_start = [chr, str(int(start) - 25), str(int(start) + 25), sv_type]
                sublist_end = [chr, str(int(end) - 25), str(int(end) + 25), sv_type]
                sublist_start.append('_'.join(sublist_start))
                sublist_end.append('_'.join(sublist_end))
                complete_list.extend((sublist_start, sublist_end))
                number_inv += 1
            elif 'NA19240' in info_field_str and not 'IMPRECISE' in info_field_str:
                if sv_type == 'INS' or sv_type == 'DUP':
                    sublist = [chr, str(int(start) - 25), str(int(end) + 25), 'INS']
                    sublist.append('_'.join(sublist))
                    complete_list.append(sublist)
                    number_ins += 1
                elif sv_type == 'DUP':
                    sublist = [chr, start, end, sv_type]
                    sublist.append('_'.join(sublist))
                    complete_list.append(sublist)
                    number_dup += 1
                elif sv_type == 'DEL':
                    sublist = [chr, start, end, sv_type]
                    sublist.append('_'.join(sublist))
                    complete_list.append(sublist)
                    number_del += 1
                else:
                    sublist = [chr, start, end, sv_type]
                    sublist.append('_'.join(sublist))
                    print(sublist)
                    complete_list.append(sublist)
                    number_nes += 1





















elif sv_method == 'manta':

    variant_list = []

    with open(str(args.path), 'r') as file:
        for line in file:
            if not line.lstrip().startswith('#'):
                variant_list.append(line.strip().split('\t'))
    file.close()


    def get_end_position(variant_list_entry):
        info_field = variant_list_entry[7].split(';')
        end_of_variant = info_field[0].translate(DD)
        return end_of_variant


    def get_sv_type(variant_list_entry):
        info_field = variant_list_entry[7].split(';')
        for item in info_field:
            if 'SVTYPE' in item:
                sv_type = item[7:]
        return sv_type


    complete_list = []

    for entry in variant_list:
        info_field = entry[7].split(';')
        info_field_str = ''.join(info_field)
        if not 'IMPRECISE' in info_field_str and not 'BND' in info_field[0] and not '_' in entry[0] and not 'EBV' in \
                                                                                                            entry[0]:
            chr = entry[0]
            start = entry[1]
            end = get_end_position(entry)
            sv_type = get_sv_type(entry)
            if sv_type == 'INS' or sv_type == 'DUP':
                if (int(end) - int(start)) < 51:
                    sv_len = int(end) - int(start)
                    to_add = int((50 - sv_len) / 2)
                    sublist = []
                    if (sv_len % 2) == 0:
                        sublist = [chr, str(int(start) - to_add), str(int(end) + to_add + 1), 'INS']
                        sublist.append('_'.join(sublist))  #
                    else:
                        sublist = [chr, str(int(start) - to_add), str(int(end) + to_add), 'INS']
                        sublist.append('_'.join(sublist))
                    complete_list.append(sublist)
                else:
                    sublist = [chr, start, end, 'INS']
                    sublist.append('_'.join(sublist))
                    complete_list.append(sublist)
                number_ins += 1
            elif sv_type == 'INV':
                sublist_start = [chr, str(int(start) - 25), str(int(start) + 25), sv_type]
                sublist_end = [chr, str(int(end) - 25), str(int(end) + 25), sv_type]
                sublist_start.append('_'.join(sublist_start))
                sublist_end.append('_'.join(sublist_end))
                complete_list.extend((sublist_start, sublist_end))
                number_inv += 1
            elif sv_type == 'DEL':
                sublist = [chr, start, end, sv_type]
                sublist.append('_'.join(sublist))
                complete_list.append(sublist)
                number_del += 1

















elif sv_method == 'pbsv':

    variant_list = []

    with open(str(args.path), 'r') as file:
        for line in file:
            if not line.lstrip().startswith('#'):
                variant_list.append(line.strip().split('\t'))
    file.close()


    def get_end_position(variant_list_entry):
        info_field = variant_list_entry[7].split(';')
        end_of_variant = info_field[1].translate(DD)
        return end_of_variant


    def get_sv_type(variant_list_entry):
        info_field = variant_list_entry[7].split(';')
        for item in info_field:
            if 'SVTYPE' in item:
                sv_type = item[7:]
        return sv_type


    complete_list = []

    for entry in variant_list:
        info_field = entry[7].split(';')
        info_field_str = ''.join(info_field)
        if not 'IMPRECISE' in info_field_str and not 'BND' in info_field[0] and not '_' in entry[0] and not 'EBV' in \
                                                                                                            entry[0]:
            chr = entry[0]
            start = entry[1]
            end = get_end_position(entry)
            sv_type = get_sv_type(entry)
            if sv_type == 'INS':
                if (int(end) - int(start)) < 51:
                    sv_len = int(get_end_position(entry)) - int(entry[1])
                    to_add = int((50 - sv_len) / 2)
                    sublist = []
                    if (sv_len % 2) == 0:
                        sublist = [chr, str(int(start) - to_add), str(int(end) + to_add + 1), sv_type]
                        sublist.append('_'.join(sublist))
                    else:
                        sublist = [chr, str(int(start) - to_add), str(int(end) + to_add), sv_type]
                        sublist.append('_'.join(sublist))
                    complete_list.append(sublist)
                else:
                    sublist = [chr, start, end, sv_type]
                    sublist.append('_'.join(sublist))
                    complete_list.append(sublist)
                number_ins += 1
            elif sv_type == 'INV':
                sublist_start = [chr, str(int(start) - 25), str(int(start) + 25), sv_type]
                sublist_end = [chr, str(int(end) - 25), str(int(end) + 25), sv_type]
                sublist_start.append('_'.join(sublist_start))
                sublist_end.append('_'.join(sublist_end))
                complete_list.extend((sublist_start, sublist_end))
                number_inv += 1
            elif sv_type == 'INS' and (int(end) - int(start)) < 10:
                sublist = [chr, str(int(start) - 25), str(int(end) + 25), sv_type]
                sublist.append('_'.join(sublist))
                complete_list.append(sublist)
                number_ins += 1
            elif sv_type == 'DUP':
                sublist = [chr, start, end, sv_type]
                sublist.append('_'.join(sublist))
                complete_list.append(sublist)
                number_dup += 1
            elif sv_type == 'DEL':
                sublist = [chr, start, end, sv_type]
                sublist.append('_'.join(sublist))
                complete_list.append(sublist)
                number_del += 1
            else:
                sublist = [chr, start, end, sv_type]
                sublist.append('_'.join(sublist))
                complete_list.append(sublist)
                number_nes += 1













elif sv_method == 'delly':

    variant_list = []

    with open(str(args.path), 'r') as file:
        for line in file:
            if not line.lstrip().startswith('#'):
                variant_list.append(line.strip().split('\t'))
    file.close()


    def get_end_position(variant_list_entry):
        info_field = variant_list_entry[7].split(';')
        end_of_variant = info_field[4].translate(DD)
        return end_of_variant


    def get_sv_type(variant_list_entry):
        info_field = variant_list_entry[7].split(';')
        for item in info_field:
            if 'SVTYPE' in item:
                sv_type = item[7:]
        return sv_type


    complete_list = []

    for entry in variant_list:
        info_field = entry[7].split(';')
        info_field_str = ''.join(info_field)
        # No IMPRECISE filter because too much would get filtered out
        if not 'BND' in info_field[1] and not '_' in entry[0] and not 'EBV' in entry[0]:
            chr = entry[0]
            start = entry[1]
            end = get_end_position(entry)
            sv_type = get_sv_type(entry)
            if sv_type == 'INS':
                if (int(get_end_position(entry)) - int(entry[1])) < 51:
                    sv_len = int(end) - int(start)
                    to_add = int((51 - sv_len) / 2)
                    sublist = []
                    if (sv_len % 2) == 0:
                        sublist = [chr, str(int(start) - to_add), str(int(end) + to_add + 1), sv_type]
                        new_sv_len = int(sublist[2]) - int(sublist[1])
                        sublist.append('_'.join(sublist))
                    else:
                        sublist = [chr, str(int(start) - to_add), str(int(end) + to_add), sv_type]
                        new_sv_len = int(sublist[2]) - int(sublist[1])
                        sublist.append('_'.join(sublist))
                    complete_list.append(sublist)
                else:
                    sublist = [chr, start, end, sv_type]
                    sublist.append('_'.join(sublist))
                    complete_list.append(sublist)
                number_ins += 1
            elif sv_type == 'INV':
                sublist_start = [chr, str(int(start) - 25), str(int(start) + 25), sv_type]
                sublist_end = [chr, str(int(end) - 25), str(int(end) + 25), sv_type]
                sublist_start.append('_'.join(sublist_start))
                sublist_end.append('_'.join(sublist_end))
                complete_list.extend((sublist_start, sublist_end))
                number_inv += 1
            elif sv_type == 'DUP':
                sublist = [chr, start, end, sv_type]
                sublist.append('_'.join(sublist))
                complete_list.append(sublist)
                number_dup += 1
            elif sv_type == 'DEL':
                sublist = [chr, start, end, sv_type]
                sublist.append('_'.join(sublist))
                complete_list.append(sublist)
                number_del += 1
            else:
                sublist = [chr, start, end, sv_type]
                sublist.append('_'.join(sublist))
                complete_list.append(sublist)
                number_nes += 1













elif sv_method == 'cnmops':
    count = 0
    for entry in variant_list:
        if 100 <= int(entry[4]) <= 900:
            if int(entry[3].translate(DD)) > 2 and int(entry[2]) - int(entry[1]) >= 50:
                sublist = [entry[0], entry[1], entry[2], 'INS']
                sublist.append('_'.join(sublist))
                number_del += 1
            elif int(entry[3].translate(DD)) < 2:
                sublist = [entry[0], entry[1], entry[2], 'DEL']
                sublist.append('_'.join(sublist))
                number_del += 1
            complete_list.append(sublist)


        else:
            count += 1

    # print(count)
    # print(len(variant_list))















elif sv_method == 'svaba':

    import re
    import os

    header_list = []
    variant_list = []

    if __name__ == "__main__":
        file = str(args.path)
        if not os.path.exists(file):
            raise IOError(file)
        alt_index = -1
        # generate mate:mate dictionary
        # load file into ram
        vcf_file = []
        with open(file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                vcf_file.append(line)
        matesDict = svaba_script.makeMateDict(vcf_file)
        with open(file, "r") as f:
            for line in f:
                # print comments
                if line.startswith("##"):
                    header_list.append(line)
                    continue
                # header contains indexes
                if line.startswith('#'):
                    split = line.split("\t")
                    for index, val in enumerate(split):
                        if val == "ALT":
                            alt_index = index
                            break
                    header_list.append(line)
                    continue
                if alt_index == -1:
                    print("ERROR: NO ALT INDEX FOUND")
                    exit(1)
                newType = svaba_script.classify(line, alt_index, matesDict)
                if newType != "NONE":
                    newLine = re.sub(r'SVTYPE=BND', "SVTYPE=" + newType, line)
                    variant_list.append(newLine.split('\t'))

    for entry in variant_list:
        sv_type = entry[7].split(';')[-1].replace('SVTYPE=', '')
        if not 'BND' in sv_type and not 'EBV' in entry[0] and not '_' in entry[0]:
            chr = entry[0]
            start = entry[1]
            end = (entry[4].split(':')[1].translate(DD))
            if sv_type == 'DEL':
                sublist = [chr, start, end, sv_type]
                sublist.append('_'.join(sublist))
                number_del += 1
                complete_list.append(sublist)
            elif sv_type == 'INS' or sv_type == 'DUP':
                sublist = [chr, start, end, sv_type]
                sublist.append(chr + '_' + start + '_' + end + '_' + sv_type)
                number_ins += 1
                complete_list.append(sublist)
            elif sv_type == 'INV':
                sublist_start = [chr, str(int(start) - 25), str(int(start) + 25), sv_type]
                sublist_start.append('_'.join(sublist_start))
                sublist_end = [chr, str(int(end) - 25), str(int(end) + 25), sv_type]
                sublist_end.append('_'.join(sublist_end))
                complete_list.extend((sublist_start, sublist_end))






















































else:
    print('This is not a tool in the list. Please use one of the following tools:')
    for software in software_list:
        print(software)

if args.stats == True:

    print('\n')
    print('Total number of variants:', len(complete_list))
    print('Number of DEL:', number_del)
    print('Number of INS:', number_ins)
    print('Number of DUP:', number_dup)
    print('Number of INV:', number_inv)
    print('Number of NES:', number_nes)

    '''
    del_array = []

    for i in range(0,len(complete_list)):
        e = complete_list[i]
        for j in range(i+1,len(complete_list)):
            ec = complete_list[j]
            if e[0] == ec[0] and int(e[1]) <= int(ec[1]) <= int(e[2]) and int(e[1]) <= int(ec[2]) <= int(e[2]) and e[3] == ec[3]:
                print(e[1], ec[1], e[3], 'and', e[2], ec[2], ec[3])
            elif e[0] == ec[0] and int(ec[1]) <= int(e[1]) <= int(ec[2]) and int(ec[1]) <= int(e[2]) <= int(ec[2]) and e[3] == ec[3]:
                print(e[1], ec[1], e[3], 'and', e[2], ec[2], ec[3])


    for entry in del_array:
        del complete_list[entry]

    '''
else:

    if str(args.format) == 'bed':
        for entry in complete_list:
            if not entry[3] == 'CNV' and not 'HLA' in entry[0]:
                print('\t'.join(entry))

    elif str(args.format) == 'csv':
        for entry in complete_list:
            if not entry[3] == 'CNV' and not 'HLA' in entry[0]:
                print(','.join(entry))
    else:
        print('This is not a useable format')





