# To Do
# Are there other seperators in the info field than ;


# Imports
import click
from .identify_caller import prepare_files
from .identify_multiple_caller import prepare_multiple_files
import pybedtools
import glob, os

# Command line options
@click.group()
@click.option('--query', help='Path to query vcf')
@click.option('--mq', help='Path to multiple query vcfs', multiple=True)
@click.option('--truth', help='Path to truth vcf')
@click.option('--ts', help='Specify truthset', required=False, default='')
@click.pass_context

# Main function that creates a ctx object which will come in handy later and takes the command line options as input
def main(ctx, query: str, truth:str, ts: str, mq: str):

    if mq:

        os.system('mkdir temp_bed_files')

        query_count = 1

        for query in mq:
            prepare_multiple_files(query, query_count, 'query', query_count)
            query_count += 1

        prepare_files(truth, ts, 'truth')

        ctx.obj['mq'] = mq


    else:

        # This folder will contain the temporary bed files that come from the VCF file
        os.system('mkdir temp_bed_files')

        # This creates bed files for the query
        query_bed_file = prepare_files(query, ts, 'query')
        ctx.obj['query'] = query_bed_file

        # This creates bed files for the truth
        truth_bed_file = prepare_files(truth, ts, 'truth')
        ctx.obj['truth'] = truth_bed_file

@main.command()
@click.pass_context
def benchmark(ctx):

    query_dict = {}
    truth_dict = {}

    os.chdir('temp_bed_files')

    # This goes through all query files and merges similar variants inside each file
    # The result becomes a value in a dictionary with the key being the SV type
    for file in glob.glob('*query.bed'):
        if not os.stat(file).st_size == 0:
            temp_file = pybedtools.BedTool(file)
            merged_temp_file = temp_file.sort().merge(d=50, c=[4, 4, 5], o=['count', 'collapse', 'collapse'])
            query_dict[file.split('_')[0]] = merged_temp_file


    # This goes through all truth files and merges similar variants inside each file
    # The result becomes a value in a dictionary with the key being the SV type
    for file in glob.glob('*truth.bed*'):
        if not os.stat(file).st_size == 0:
            temp_file = pybedtools.BedTool(file)
            merged_temp_file = temp_file.sort().merge(d=50, c=[4, 4, 5], o=['count', 'collapse', 'collapse'])
            truth_dict[file.split('_')[0]] = merged_temp_file

    query_keys = list(query_dict.keys())
    truth_keys = list(truth_dict.keys())
    intersect_keys = list(set(truth_keys) & set(query_keys))

    print()

    tp_all = 0
    fp_all = 0
    fn_all = 0

    for sv_type in intersect_keys:
        intersect = query_dict[sv_type].intersect(truth_dict[sv_type], r=True, f=0.51, wa=True, wb=True)
        print('Number of ' + sv_type + 's in query: ', len(query_dict[sv_type]))
        print('Number of ' + sv_type + 's in truthset: ', len(truth_dict[sv_type]))
        print('Number of ' + sv_type + 's in intersect: ', len(intersect))
        print()
        tp = len(intersect)
        tp_all += tp
        fp = len(truth_dict[sv_type]) - len(intersect)
        fp_all += fp
        fn = len(query_dict[sv_type]) - len(intersect)
        fn_all += fn
        print('TP (' + sv_type + '): ', tp)
        print('FP (' + sv_type + '): ', fp)
        print('FN (' + sv_type + '): ', fn)
        print()
        print('Precision (' + sv_type + '): ', round(tp / (tp + fp), 3))
        print('Recall (' + sv_type + '): ', round(tp / (tp + fn), 3))
        print('F1 (' + sv_type + '): ', round((2 * (tp / (tp + fp)) * (tp / (tp + fn))) / (
                    (tp / (tp + fp)) + (tp / (tp + fn))), 3))
        print()

    print('TP (ALL): ', tp_all)
    print('FP (ALL): ',fp_all)
    print('FN (ALL): ',fn_all)
    print()
    print('Precision (ALL): ', round(tp_all/(tp_all+fp_all),3))
    print('Recall (ALL): ', round(tp_all/(tp_all+fn_all),3))
    print('F1 (ALL): ', round((2*(tp_all/(tp_all+fp_all))*(tp_all/(tp_all+fn_all)))/((tp_all/(tp_all+fp_all))+(tp_all/(tp_all+fn_all))), 3))
    print()

    os.chdir('..')
    os.system('rm -r temp_bed_files')

@main.command()
@click.pass_context
def multibenchmark(ctx):

    query_list = []
    truth_dict = {}
    naming_order = []

    os.chdir('temp_bed_files')

    # This goes through all query files and merges similar variants inside each file
    # The result becomes a value in a dictionary with the key being the SV type
    for subfolder in glob.glob('*/'):
        query_dict = {}
        for file in glob.glob(subfolder + '/*query*'):
            if not os.stat(file).st_size == 0:
                temp_file = pybedtools.BedTool(file)
                merged_temp_file = temp_file.sort().merge(d=50, c=[4, 4, 5], o=['count', 'collapse', 'collapse'])
                query_dict[file.split('/')[1].split('_')[0]] = merged_temp_file
                naming_order.append(int(file.split('/')[0].split('_')[1]))
                #query_dict[file.split('/')[1].split('_')[0] + '_' + file.split('/')[0].split('_')[1]] = merged_temp_file
        query_list.append(query_dict)

    # This goes through all truth files and merges similar variants inside each file
    # The result becomes a value in a dictionary with the key being the SV type
    for file in glob.glob('*truth.bed*'):
        if not os.stat(file).st_size == 0:
            temp_file = pybedtools.BedTool(file)
            merged_temp_file = temp_file.sort().merge(d=50, c=[4, 4, 5], o=['count', 'collapse', 'collapse'])
            truth_dict[file.split('_')[0]] = merged_temp_file

    mq_count = 0
    ordered_naming_order = list(dict.fromkeys(naming_order))

    for current_dict in query_list:
        query_keys = list(current_dict.keys())
        truth_keys = list(truth_dict.keys())
        intersect_keys = list(set(truth_keys) & set(query_keys))

        print('Query file used: ', ctx.obj['mq'][ordered_naming_order[mq_count]-1])
        print()

        tp_all = 0
        fp_all = 0
        fn_all = 0

        for sv_type in intersect_keys:
            intersect = current_dict[sv_type].intersect(truth_dict[sv_type], r=True, f=0.51, wa=True, wb=True)
            print('Number of ' + sv_type + 's in query: ', len(current_dict[sv_type]))
            print('Number of ' + sv_type + 's in truthset: ', len(truth_dict[sv_type]))
            print('Number of ' + sv_type + 's in intersect: ', len(intersect))
            print()
            tp = len(intersect)
            tp_all += tp
            fp = len(truth_dict[sv_type]) - len(intersect)
            fp_all += fp
            fn = len(current_dict[sv_type]) - len(intersect)
            fn_all += fn
            print('TP (' + sv_type + '): ', tp)
            print('FP (' + sv_type + '): ', fp)
            print('FN (' + sv_type + '): ', fn)
            print()
            print('Precision (' + sv_type + '): ', round(tp / (tp + fp), 3))
            print('Recall (' + sv_type + '): ', round(tp / (tp + fn), 3))
            print('F1 (' + sv_type + '): ', round((2 * (tp / (tp + fp)) * (tp / (tp + fn))) / (
                    (tp / (tp + fp)) + (tp / (tp + fn))), 3))
            print()

        print('TP (ALL): ', tp_all)
        print('FP (ALL): ', fp_all)
        print('FN (ALL): ', fn_all)
        print()
        print('Precision (ALL): ', round(tp_all / (tp_all + fp_all), 3))
        print('Recall (ALL): ', round(tp_all / (tp_all + fn_all), 3))
        print('F1 (ALL): ', round((2 * (tp_all / (tp_all + fp_all)) * (tp_all / (tp_all + fn_all))) / (
                    (tp_all / (tp_all + fp_all)) + (tp_all / (tp_all + fn_all))), 3))

        print()
        print('########################################')
        print()

        mq_count += 1

    os.chdir('..')
    os.system('rm -r temp_bed_files')

def start():
    main(obj={})

if __name__ == '__main__':
    start()
