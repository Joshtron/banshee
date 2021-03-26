import string

class Del:
    def __init__(self, keep=string.digits):
        self.comp = dict((ord(c), c) for c in keep)

    def __getitem__(self, k):
        return self.comp.get(k)

DD = Del()

def manta_parser(path):

    complete_list = []

    number_del = 0
    number_ins = 0
    number_dup = 0
    number_inv = 0
    number_nes = 0

    header_list = []

    variant_list = []

    with open(str(path), 'r') as file:
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

    return(complete_list)