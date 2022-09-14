"""
Split a single CSV file into several, each containing particles of a single density tomogram
"""

__author__ = 'Antonio Martínez-Sánchez'


import sys, getopt, os, time, csv
import pandas as pd


def load_csv_into_tomo_tables(in_csv_file):
    """
    Load a CSV file as a dictionary of tables, one for each density
    :param in_csv_file: input CSV file path
    :return: a dictionary where each density path is an entry for a table, each table contains all particles of single
             density
    """
    tables_df = pd.read_csv(in_csv_file, delimiter='\t', names=['Density Micrographs', 'PolyData', 'Tomo3D', 'Type',
                                                                 'Label', 'Code', 'Polymer', 'X', 'Y', 'Z',
                                                                 'Q1', 'Q2', 'Q3', 'Q4'], header=0)
    den_tomos = set(tables_df['Tomo3D'].tolist())
    tables_dic = dict().fromkeys(den_tomos)
    for key in tables_dic:
        tables_dic[key] = dict().fromkeys(tables_df.columns.tolist())
        for kkey in tables_dic[key]:
            tables_dic[key][kkey] = list()
    for row in tables_df.iterrows():
        key = row[1]['Tomo3D']
        for item, value in row[1].items():
            tables_dic[key][item].append(value)
    return tables_dic


def write_table(table, out_file):
    """
    Store a table in a CSV file
    :param table: input table dictionary
    :param out_file: path for the output file
    :return:
    """
    with open(out_file, 'w', newline='') as csv_file:
        fieldnames = list(table.keys())
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for row in range(len(table[fieldnames[0]])):
            dic_row = dict().fromkeys(fieldnames)
            for key in fieldnames:
                dic_row[key] = table[key][row]
            writer.writerow(dic_row)


def main(argv):

    # Input parsing
    in_csv, out_dir = None, None
    try:
        opts, args = getopt.getopt(argv, 'hi:o:', ['help', 'icsv', 'odir'])
    except getopt.GetoptError:
        print('python split_csv.py -i <in_csv> -o <out_dir>')
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print('python split_csv.py -i <in_csv> -o <out_dir>')
            print('\t-i (--itomo) <in_csv>: input CSV file')
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
        print('python split_csv.py -i <in_csv> -o <out_dir>')
        print('\t-i (--itomo) <in_csv>: input CSV file')
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