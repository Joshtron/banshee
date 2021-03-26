import string

class Del:
    def __init__(self, keep=string.digits):
        self.comp = dict((ord(c), c) for c in keep)

    def __getitem__(self, k):
        return self.comp.get(k)

DD = Del()

def nstd152_parser(path, ts):

    complete_list = []

    number_del = 0
    number_ins = 0
    number_dup = 0
    number_inv = 0
    number_nes = 0

    header_list = []

    variant_list = []

    if ts == '':
        print('\nPlease specify a valid truthset\n',
              'Yorubian Trio: NA19238, NA19239, NA19240\n',
              'Han chinese Trio: HG00512, HG00513, HG00514\n',
              'Ashkenazim Trio: HG00731, HG00732, HG00733\n')
    else:
        variant_list = []

        # This opens the file and stores every line that does not strart with '#' in a list. With that I can skip the header
        with open(path, 'r') as file:
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
            elif ts in info_field_str and not 'IMPRECISE' in info_field_str:
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

    return(complete_list)