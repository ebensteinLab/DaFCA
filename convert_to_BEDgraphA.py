from parsing.plotting import save_image, export_csv
from parsing_scripts import get_df_from_file
import argparse
import os


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='convert to BEDgraph')
    parser.add_argument('-fp','--file_path', required=True,
                        help='path to input file, the raw data from checkProfileAlignment', default = '')
    parser.add_argument('-bp', '--BEDgraph_path', required=False,
                        help='path for results file', default ='C:/Users/irys/Downloads/lab/input_examples')
    parser.add_argument('-chr','--chromosome_number', required=True,
                        help='number of chromosome with prefix of chr, meaning for chromosome number 1 insert chr1', default='chr1')
    args = parser.parse_args()
    chr = args.chromosome_number
    bedgraph_path = args.BEDgraph_path
    df = get_df_from_file(args.file_path)
    df = df.fillna(value=0)
    df['Intensity'] = df.groupby(['molID', 'Start Locus'])['Intensity'].transform('mean') # only avergae for the same locus and same molecule
    df['Num of stitching'] = df.groupby(['molID', 'Start Locus'])['Start Locus'].transform('size')
    df = df.drop_duplicates(['Start Locus', 'molID'])
    df['Start Locus'] = df['Start Locus']
    df['Start Locus'] = df['Start Locus'] # needed to do this way in order to avoid floating point inaccuracies
    df['Start Locus'] = df['Start Locus'].astype(int)
    df.insert(0, 'Chromosome', chr)
    df.insert(2, 'End Locus', df['Start Locus'] + 500)
    df['Start Locus'] = df['Start Locus'] - 500
    df = df.sort_values('Start Locus')
    results_name = chr+' results.BEDgraph'
    df.set_index('Chromosome', inplace=True)
    df.to_csv(bedgraph_path+'/'+results_name, sep='\t')
