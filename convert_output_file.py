from parsing.plotting import save_image, export_csv
from parsing_scripts import get_df_from_file
import argparse
import os

FILE_PATH = "input_examples/example1.txt"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--file_path', required=True, help='path to output file')
    parser.add_argument('--csv_path', required=False,
                        help='path to csv if not set csv wont be generated')
    parser.add_argument('--csv_output_prefix', default="result", required=False,
                        help='the prefix name of the csv generated')
    parser.add_argument('--visualize', default=False, action='store_true',
                        help='set to true if you want to see the graphs')
    args = parser.parse_args()
    if not args.csv_path:
        csv_path = os.getcwd()
        print(f"CSV path is not specified will store results in {csv_path}\{args.csv_output_prefix}")
    else:
        csv_path = args.csv_path
    csv_output_prefix = args.csv_output_prefix
    df = get_df_from_file(args.file_path)
    df['mean'] = df.groupby('X')['Y'].transform('mean')
    df['std'] = df.groupby('X')['Y'].transform('std')
    df = df.fillna(value=0)
    df['X_OCCURRENCES'] = df.groupby(['X'])['X'].transform('size')
    df_for_res = df.drop_duplicates(['X'])
    del df_for_res['Y']
    df_for_res = df_for_res.sort_values('X')
    export_csv(df_for_res, csv_path, args.visualize, file_name=f"{csv_output_prefix}.csv")
    df = df.sort_values('X')

    export_csv(df, csv_path, args.visualize, file_name=f"{csv_output_prefix}_raw.csv")
    save_image(df_for_res, csv_path)
