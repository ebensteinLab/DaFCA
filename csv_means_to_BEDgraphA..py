from parsing.plotting import save_image, export_csv
from parsing_scripts import get_df_from_file
import argparse
import pandas as pd
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-fp','--file_path', required=True,
                        help='path to input file', default = "")
    parser.add_argument('-bp', '--BEDgraph_path', required=False,
                        help='path to directory of output', default='.')
    parser.add_argument('-chr','--chromosome_number', required=True,
                        help='number of chromosome with prefix of chr, meaning for chromosome number 1 insert chr1', default='chr1')
    parser.add_argument('-mean1', '--mean_filter_lower', required=False,
                        help='show only pixels that have mean above this value, included', default=-100000)
    parser.add_argument('-mean2', '--mean_filter_upper', required=False,
                        help='show only pixels that have mean below this value, included', default=100000)
    parser.add_argument('-std1','--std_filter_lower', required=False,
                        help='show only pixels that have std above this value, included', default=-100000)
    parser.add_argument('-std2','--std_filter_upper', required=False,
                        help='show only pixels that have std below this value, included', default=100000)
    parser.add_argument('-occ1', '--occ_filter_lower', required=False,
                        help='show only pixels that have occurrences above this value, included', default=-100000)
    parser.add_argument('-occ2', '--occ_filter_upper', required=False,
                        help='show only pixels that have occurrences below this value, included', default=100000)
    parser.add_argument('-show', '--show_std_and_occ', required=False,
                        help='if this is True then the resulted BEDgraph will have columns of the std and occurences', default=False)
    parser.add_argument('-name', '--file_name', required=False,
                        help='write name of outputfile', default='filtered_results')	
    args = parser.parse_args()
    csv_path = args.file_path
    bed_path = args.BEDgraph_path
    chr = args.chromosome_number
    mean_up = int(args.mean_filter_upper)
    mean_low = int(args.mean_filter_lower)
    std_up = int(args.std_filter_upper)
    std_low = int(args.std_filter_lower)
    occ_up = int(args.occ_filter_upper)
    occ_low = int(args.occ_filter_lower)
    name= args.file_name
    show = args.show_std_and_occ
    df = pd.read_csv(csv_path)
    df_filtered = df[(df['mean'] >= mean_low) & (df['mean'] <= mean_up) & (df['std'] >= std_low) &(df['std'] <= std_up) &
                     (df['X_OCCURRENCES'] >= occ_low) & (df['X_OCCURRENCES'] <= occ_up)]
    df_filtered['X'] = df_filtered['X'].astype(int)
    df_filtered.insert(0, 'Chromosome', chr)
    df_filtered.insert(2, 'End Locus', df_filtered['X'] + 250)
    df_filtered['X'] = df_filtered['X'] - 250
    df_filtered = df_filtered.sort_values('X')
    df_filtered.rename(columns={'X': 'Start Locus', 'mean':'Intensity mean'}, inplace=True)
    if show is False:
        df_filtered = df_filtered.drop(columns=['std', 'X_OCCURRENCES'])
    results_name = (chr + name + '.BEDgraph')
    df_filtered.set_index('Chromosome', inplace=True)
    df_filtered.to_csv(bed_path+'/'+results_name+'', sep='\t', header=None)


