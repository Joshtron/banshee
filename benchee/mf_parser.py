import string

class Del:
    def __init__(self, keep=string.digits):
        self.comp = dict((ord(c), c) for c in keep)

    def __getitem__(self, k):
        return self.comp.get(k)

DD = Del()

def mf_parser(path):

    complete_list = []

    number_del = 0
    number_ins = 0
    number_dup = 0
    number_inv = 0
    number_nes = 0
    number_of_imprecise_EBV = 0

    variant_list = []

    # This opens the file and stores every line that does not start with '#' in a list. With that I can skip the header
    with open(path, 'r') as file:
        for line in file:
            if not line.lstrip().startswith('#'):
                variant_list.append(line.strip().split('\t'))
    file.close()

    # This function "parser_infoflied"  filters for the end position, type and rep_type of a structural variant, by searching
    # the list info_field for the strings "END", "SVYTPE" and "REPYTPE". If there is no "REPYTPE" in the list, it will
    # be substituted by a ".".
    # BND --> noch den code anpassen, da SVLEN nicht immer der selbe Index 
    # If the SV type is 'BND' the info field SVLENGTH gets filtered for the number and added to the start of the variant

    def parser_info_field(variant_list_entry):
        info_field = variant_list_entry[7].split(';')
        rep_type_position = variant_list_entry[7]
        if 'SVTYPE=BND' in variant_list_entry[7]:
            end_of_variant = str(int(variant_list_entry[1]) + int(info_field[34].translate(DD)))
        else:
            for entry in info_field:
                if entry.startswith('END'):
                    end_of_variant =entry.split('=')[1]
                if 'SVTYPE' in entry:
                    sv_type =entry.split('=')[1]
                if 'REPTYPE' in rep_type_position:
                    rep_type_index = rep_type_position.index('REPTYPE')
                    rep_type = rep_type_position[rep_type_index:].split(';')[0].split('=')[1]
                else:
                    rep_type = '.'
        return end_of_variant,sv_type,rep_type



    for entry in variant_list:
        # If-command filters out the decoys as I tought they might cause problems. Delete if they are needed
        # Filters also the strange variants at the end that have a altered INFO field. Functions above do not work or
        # need to be altered. We should first talk what they are exactly as I have never heard that sv type
        # Also filters out EBV variants
        if not 'SVTYPE=BND' in entry[7] and not '_' in entry[0] and not 'EBV' in entry[0] and not 'IMPRECISE' in entry[0]:  # gelöscht: "not 'IMPRECISE' in info_field[0] and" nach dem ersten if
            chr = entry[0]
            start = entry[1]
            end,sv_type,rep_type= parser_info_field(entry)
            if sv_type == 'INS' or sv_type == 'DUP':
                if rep_type == 'DUP' or sv_type == 'DUP':
                    sv_type = 'DUP'
                    number_dup += 1
                else:
                    number_ins += 1
                sv_len = int(end) - int(start)
                if sv_len < 51:
                    to_add = int((51 - sv_len) / 2)
                    if (sv_len % 2) == 0:
                        sublist = [chr, str(int(start) - to_add), str(int(end) + to_add + 1), sv_type, rep_type]  # rep_type am 28.05.21
                    else:
                        sublist = [chr, str(int(start) - to_add), str(int(end) + to_add), sv_type, rep_type]  # rep_type am 28.05.21
                    sublist.append('_'.join(sublist))
                    complete_list.append(sublist)
                else:
                    sublist = [chr, start, end, sv_type, rep_type] # rep_type am 28.05.21
                    tag_list =[ chr, start, end, sv_type] #hier rep_type aus dem Tag aber nicht aus der Liste wir können es auch vollständig für alle rausnehmen.Hier nur versuchsweise.
                    sublist.append('_'.join(tag_list))
                    complete_list.append(sublist)
            elif sv_type == 'INV':
                sublist_start = [chr, str(int(start) - 25), str(int(start) + 25), sv_type,rep_type]  # rep_type am 28.05.21
                sublist_end = [chr, str(int(end) - 25), str(int(end) + 25), sv_type, rep_type]  # rep_type am 28.05.21
                sublist_start.append('_'.join(sublist_start))
                sublist_end.append('_'.join(sublist_end))
                complete_list.extend((sublist_start, sublist_end))
                number_inv += 1
            elif sv_type == 'DEL':
                sublist = [chr, start, end, sv_type, rep_type]  # rep_type am 28.05.21
                sublist.append('_'.join(sublist))
                complete_list.append(sublist)
                number_del += 1
            else:
                sublist = [chr, start, end, sv_type, rep_type]  # rep_type am 28.05.21
                sublist.append('_'.join(sublist))
                complete_list.append(sublist)
                number_nes += 1
        else:
            number_of_imprecise_EBV += 1

    return (complete_list)
