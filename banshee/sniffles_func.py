import string

class Del:
    def __init__(self, keep=string.digits):
        self.comp = dict((ord(c), c) for c in keep)

    def __getitem__(self, k):
        return self.comp.get(k)

DD = Del()

def sniffles_parser(path):

    complete_list = []

    number_del = 0
    number_ins = 0
    number_dup = 0
    number_inv = 0
    number_nes = 0

    header_list = []

    variant_list = []

    #This opens the file and stores every line that does not strart with '#' in a list. With that I can skip the header
    with open(path, 'r') as file:
        for line in file:
            if not line.lstrip().startswith('#'):
                variant_list.append(line.strip().split('\t'))
    file.close()

    #This filters for the end position of a variant as the position of that in the .vcf file is 'always' the same when
    #coming from sniffles
    #If the SV type is 'BND' the info field SVLENGTH gets filtered for the number and added to the start of the variant
    def get_end_position(variant_list_entry):
        if 'SVTYPE=BND' in variant_list_entry[7]:
            info_field = variant_list_entry[7].split(';')
            end_of_variant = str(int(variant_list_entry[1]) + int(info_field[8].translate(DD)))
        else:
            end_position = variant_list_entry[7].split(';')[3]
            end_of_variant = end_position[4:]
        return end_of_variant

    #Same for the sv-type.
    #
    def get_sv_type(variant_list_entry):
        if not 'SVTYPE=BND' in variant_list_entry[7]:
            end_position = variant_list_entry[7].split(';')[8]
            sv_type = end_position[7:]
            return sv_type


    number_of_imprecise_EBV = 0

    for entry in variant_list:
        info_field = entry[7].split(';')
        #If-command filters out the decoys as I tought they might cause problems. Delete if they are needed
        #Filters also the strange variants at the end that have a altered INFO field. Functions above do not work or
        #need to be altered. We should first talk what they are exactly as I have never heard that sv type
        #Also filters out EBV variants
        if not 'IMPRECISE' in info_field[0] and not 'SVTYPE=BND' in entry[7] and not '_' in entry[0] and not 'EBV' in entry[0]:
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
                        sublist = [chr, str(int(start) - to_add), str(int(end) + to_add +1), 'INS']
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

    return(complete_list)
