import string

class Del:
    def __init__(self, keep=string.digits):
        self.comp = dict((ord(c), c) for c in keep)

    def __getitem__(self, k):
        return self.comp.get(k)

DD = Del()

def cnmops_parser(path):

    variant_list = []

    with open(path, 'r') as file:
        for line in file:
            if not line.lstrip().startswith('#'):
                variant_list.append(line.strip().split('\t'))
    file.close()

    complete_list = []

    count = 0
    for entry in variant_list:
        if 100 <= int(entry[4]) <= 900:
            if int(entry[3].translate(DD)) > 2 and int(entry[2])-int(entry[1]) >= 50:
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

    return(complete_list)