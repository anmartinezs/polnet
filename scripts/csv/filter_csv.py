"""
Filter a CSV file with particles according to a column (Type, Code, ...) value.
For example, this script allow to extract all particles of specific type, or all macromolecules with a specific code
"""

__author__ = 'Antonio Martínez-Sánchez'


import sys, getopt, os, time
import pandas as pd

COLUMN_HEADS = ('Density Micrographs', 'PolyData', 'Tomo3D', 'Type', 'Label', 'Code', 'Polymer', 'X', 'Y', 'Z', 'Q1', 'Q2', 'Q3', 'Q4')

def main(argv):

    # Input parsing
    in_csv, out_csv = None, None
    in_column, in_value = None, None
    try:
        opts, args = getopt.getopt(argv, 'hi:o:c:v:', ['help', 'icsv', 'ocsv', 'column', 'value'])
    except getopt.GetoptError:
        print('python split_csv_by_tomos.py -i <in_csv> -o <out_dir> -c <column> -v <value>')
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print('python split_csv_by_tomos.py -i <in_csv> -o <out_dir> -c <column> -v <value>')
            print('\t-i (--icsv) <in_csv>: input CSV file')
            print('\t-o (--ocsv) <in_csv>: output CSV file')
            print('\t-c (--column) <column>: column header string')
            print('\t-v (--value) <value>: filtering value')
            sys.exit()
        elif opt in ("-i", "--icsv"):
            in_csv = arg
            if not (os.path.splitext(in_csv)[1] == '.csv'):
                print('The input file must have a .csv extension!')
                sys.exit()
        elif opt in ("-o", "--ocsv"):
            out_csv = arg
            if not (os.path.splitext(out_csv)[1] == '.csv'):
                print('The output file must have a .csv extension!')
                sys.exit()
        elif opt in ("-c", "--column"):
            in_column = arg
            if not(in_column in COLUMN_HEADS):
                print('The intput', in_column,' is not a valid column!')
                sys.exit()
        elif opt in ("-v", "--value"):
            in_value = arg
    if (in_csv is None) or (out_csv is None) or (in_column is None) or (in_value is None):
        print('python split_csv_by_tomos.py -i <in_csv> -o <out_dir> -c <column> -v <value>')
        print('\t-i (--icsv) <in_csv>: input CSV file')
        print('\t-o (--ocsv) <in_csv>: output CSV file')
        print('\t-c (--column) <column>: column header string')
        print('\t-v (--value) <value>: filtering value')
        sys.exit()

    # Load the csv as a Pandas DataFrame
    df = pd.read_csv(in_csv, delimiter='\t', names=COLUMN_HEADS, header=0)

    # Filtering
    filt_df = df.query(in_column + " == '" + in_value + "'")

    # Storing the results
    filt_df.to_csv(out_csv, sep='\t', columns=COLUMN_HEADS)

    print('Successfully terminated. (' + time.strftime("%c") + ')')


if __name__ == "__main__":
    main(sys.argv[1:])