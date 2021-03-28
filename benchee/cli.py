# Imports
import click
from .identify_caller import prepare_files
import pybedtools
import glob, os

# Command line options
@click.group()
@click.option('--query', help='Path to query vcf')
@click.option('--truth', help='Path to truth vcf')
@click.option('--ts', help='Specify truthset', required=False, default='')
@click.pass_context

# Main function that creates a ctx object which will come in handy later and takes the command line options as input
def main(ctx, query: str, truth:str, ts: str):

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
        temp_file = pybedtools.BedTool(file)
        merged_temp_file = temp_file.sort().merge(d=50, c=[4, 4, 5], o=['count', 'collapse', 'collapse'])
        query_dict[file.split('_')[0]]=merged_temp_file

    # This goes through all truth files and merges similar variants inside each file
    # The result becomes a value in a dictionary with the key being the SV type
    for file in glob.glob('*truth.bed*'):
        temp_file = pybedtools.BedTool(file)
        merged_temp_file = temp_file.sort().merge(d=50, c=[4, 4, 5], o=['count', 'collapse', 'collapse'])
        truth_dict[file.split('_')[0]]=merged_temp_file

    # Here the actual intersect is performed
    del_intersect = query_dict['del'].intersect(truth_dict['del'], r=True, f=0.51, wa=True, wb=True)
    ins_intersect = query_dict['ins'].intersect(truth_dict['ins'], r=True, f=0.51, wa=True, wb=True)
    dup_intersect = query_dict['dup'].intersect(truth_dict['dup'], r=True, f=0.51, wa=True, wb=True)
    inv_intersect = query_dict['inv'].intersect(truth_dict['inv'], r=True, f=0.51, wa=True, wb=True)

    # Cool math
    print('')
    print('Number of DELs in intersect: ', len(del_intersect))
    print('Number of INSs in intersect: ', len(ins_intersect))
    print('Number of DUPs in intersect: ', len(dup_intersect))
    print('Number of INVs in intersect: ', len(inv_intersect))
    print()
    print('Number of DELs in truthset: ', len(truth_dict['del']))
    print('Number of INSs in truthset: ', len(truth_dict['ins']))
    print('Number of DUPs in truthset: ', len(truth_dict['dup']))
    print('Number of INVs in truthset: ', len(truth_dict['inv']))
    print()
    print('Number of DELs in query: ', len(query_dict['del']))
    print('Number of INSs in query: ', len(query_dict['ins']))
    print('Number of DUPs in query: ', len(query_dict['dup']))
    print('Number of INVs in query: ', len(query_dict['inv']))
    print()

    tp_del = len(del_intersect)
    fp_del = len(truth_dict['del'])-len(del_intersect)
    fn_del = len(query_dict['del'])-len(del_intersect)

    print('TP (DEL): ', tp_del)
    print('FP (DEL): ',fp_del)
    print('FN (DEL): ',fn_del)
    print()
    print('Precision (DEL): ', round(tp_del/(tp_del+fp_del),3))
    print('Recall (DEL): ', round(tp_del/(tp_del+fn_del),3))
    print('F1 (DEL): ', round((2*(tp_del/(tp_del+fp_del))*(tp_del/(tp_del+fn_del)))/((tp_del/(tp_del+fp_del))+(tp_del/(tp_del+fn_del))), 3))
    print()

    tp_ins = len(ins_intersect)
    fp_ins = len(truth_dict['ins'])-len(ins_intersect)
    fn_ins = len(query_dict['ins'])-len(ins_intersect)

    print('TP (INS): ', tp_ins)
    print('FP (INS): ',fp_ins)
    print('FN (INS): ',fn_ins)
    print()
    print('Precision (INS): ', round(tp_ins/(tp_ins+fp_ins),3))
    print('Recall (INS): ', round(tp_ins/(tp_ins+fn_ins),3))
    print('F1 (INS): ', round((2*(tp_ins/(tp_ins+fp_ins))*(tp_ins/(tp_ins+fn_ins)))/((tp_ins/(tp_ins+fp_ins))+(tp_ins/(tp_ins+fn_ins))), 3))
    print()

    tp_inv = len(inv_intersect)
    fp_inv = len(truth_dict['inv'])-len(inv_intersect)
    fn_inv = len(query_dict['inv'])-len(inv_intersect)

    print('TP (INV): ', tp_inv)
    print('FP (INV): ',fp_inv)
    print('FN (INV): ',fn_inv)
    print()
    print('Precision (INV): ', round(tp_inv/(tp_inv+fp_inv),3))
    print('Recall (INV): ', round(tp_inv/(tp_inv+fn_inv),3))
    print('F1 (INV): ', round((2*(tp_inv/(tp_inv+fp_inv))*(tp_inv/(tp_inv+fn_inv)))/((tp_inv/(tp_inv+fp_inv))+(tp_inv/(tp_inv+fn_inv))), 3))
    print()
    
    tp_all = tp_del + tp_ins + tp_inv
    fp_all = (len(truth_dict['del'])+len(truth_dict['ins'])+len(truth_dict['inv']))-tp_all
    fn_all = (len(query_dict['del'])+len(query_dict['ins'])+len(query_dict['inv']))-tp_all

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

def start():
    main(obj={})

if __name__ == '__main__':
    start()
