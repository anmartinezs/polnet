"""
Split a single CSV file into several, each containing particles of a single density tomogram
"""

__author__ = 'Antonio Martínez-Sánchez'


import sys, getopt, os, time
from polnet.lio import load_csv_into_tomo_tables, write_table


def main(argv):

    # Input parsing
    in_csv, out_dir = None, None
    try:
        opts, args = getopt.getopt(argv, 'hi:o:', ['help', 'icsv', 'odir'])
    except getopt.GetoptError:
        print('python split_csv_by_tomos.py -i <in_csv> -o <out_dir>')
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print('python split_csv_by_tomos.py -i <in_csv> -o <out_dir>')
            print('\t-i (--icsv) <in_csv>: input CSV file')
            print('\t-o (--odir) <in_dir>: output directory')
            sys.exit()
        elif opt in ("-i", "--icsv"):
            in_csv = arg
            if not (os.path.splitext(in_csv)[1] == '.csv'):
                print('The input file must have a .csv extension!')
                sys.exit()
        elif opt in ("-o", "--odir"):
            out_dir = arg
            if not os.path.exists(out_dir):
                print('The output directory does not exist!')
                sys.exit()
    if (in_csv is None) or (out_dir is None):
        print('python split_csv_by_tomos.py -i <in_csv> -o <out_dir>')
        print('\t-i (--icsv) <in_csv>: input CSV file')
        print('\t-o (--odir) <in_dir>: output directory')
        sys.exit()

    print('\t-Loading the input CSV file:', in_csv)
    in_file = os.path.splitext(os.path.split(in_csv)[1])[0]
    tomo_tables = load_csv_into_tomo_tables(in_csv)

    print('\t-Writting tables in different CSV files:')
    for key in tomo_tables:
        sub_key = os.path.splitext(os.path.split(key)[1])[0]
        out_path = out_dir + '/' + in_file + '_' + sub_key + '.csv'
        print('\t\t+Table for density', key, 'in', out_path)
        write_table(tomo_tables[key], out_path)

    print('Successfully terminated. (' + time.strftime("%c") + ')')


if __name__ == "__main__":
    main(sys.argv[1:])