"""
Classes with the info for a synthetic tomograms
"""

__author__ = 'Antonio Martinez-Sanchez'


import csv
from polnet.network import Network


class MmerFile:
    """
    For handling protein (or monomer) files
    """

    def __init__(self, in_file):
        """
        Constructor
        :param in_file: path to the input file with extension .pns
        """
        self.__mmer_id = None
        self.__mmer_svol = None
        self.__iso = None
        self.__pmer_l = None
        self.__pmer_occ = None
        self.__pmer_l_max = None
        self.__pmer_over_tol = 0
        if in_file is not None:
            self.load_mmer_file(in_file)

    def get_mmer_id(self):
        return self.__mmer_id

    def get_mmer_svol(self):
        return self.__mmer_svol

    def get_iso(self):
        return self.__iso

    def get_pmer_l(self):
        return self.__pmer_l

    def get_pmer_occ(self):
        return self.__pmer_occ

    def get_pmer_l_max(self):
        return self.__pmer_l_max

    def get_pmer_over_tol(self):
        return self.__pmer_over_tol

    def load_mmer_file(self, in_file):
        """
        Load protein parameters from an input file
        :param in_file: path to the input file with extension .pns
        """

        assert isinstance(in_file, str) and in_file.endswith('.pns')

        # Reading input file
        with open(in_file) as file:
            for linea in file:
                if len(linea.strip()) > 0:

                    # Parsing an file entry
                    var, value = linea.split('=')
                    var = var.replace(' ', '')
                    var = var.replace('\n', '')
                    value = value.replace(' ', '')
                    value = value.replace('\n', '')
                    if var == 'MMER_ID':
                        self.__mmer_id = value
                    elif var == 'MMER_SVOL':
                        self.__mmer_svol = value
                    elif var == 'MMER_ISO':
                        self.__iso = float(value)
                    elif var == 'PMER_L':
                        self.__pmer_l = float(value)
                    elif var == 'PMER_OCC':
                        self.__pmer_occ = float(value)
                    elif var == 'PMER_L_MAX':
                        self.__pmer_l_max = float(value)
                    elif var == 'PMER_OVER_TOL':
                        self.__pmer_over_tol = float(value)
                    else:
                        print('ERROR: (MmerFile - load_protein_file) input entry not recognized:', value)


class SynthTomo:
    """
    This class model a Synthetic tomogram, tomogram's info are stored in disk to handle large datasets:
        - Ground truth: integer type tomogram associating each pixel with the id (background is always 0)
    """

    def __init__(self):
        self.__den = None
        self.__tomo = None
        self.__mic = None
        self.__poly = None
        self.__motifs = list()

    def get_den(self):
        return self.__den

    def get_mics(self):
        return self.__mic

    def get_poly(self):
        return self.__poly

    def get_tomo(self):
        return self.__tomo

    def set_den(self, den):
        assert isinstance(den, str) and den.endswith('.mrc')
        self.__den = den

    def set_tomo(self, tomo):
        assert isinstance(tomo, str) and tomo.endswith('.mrc')
        self.__tomo = tomo

    def set_mics(self, mics):
        assert isinstance(mics, str) and mics.endswith('.mrc')
        self.__mics = mics

    def set_poly(self, poly):
        assert isinstance(poly, str) and poly.endswith('.vtp')
        self.__poly = poly

    def get_motif_list(self):
        return self.__motifs

    def add_network(self, net, m_type, lbl=None, code=None):
        """
        Add al motifs within the input network to the synthetic tomogram
        :param net: a network object instance
        :param m_type: string with the type of monomers contained in the network
        :param lbl: integer label, if None (default) it is taken from monomer information
        :param code: string code for the network monomers, if None (default) it is taken from monomer information
        """
        assert issubclass(type(net), Network)
        assert isinstance(m_type, str)
        if (lbl is not None): assert isinstance(lbl, int)
        if (code is not None): assert isinstance(code, str)
        for pmer_id, pmer in enumerate(net.get_pmers_list()):
            for mmer_id in range(pmer.get_num_monomers()):
                if lbl is None:
                    hold_lbl = pmer.get_mmer_id(mmer_id)
                else:
                    hold_lbl = pmer.get_mmer_code(mmer_id)
                if code is None:
                    hold_code = pmer.get_mmer_code(mmer_id)
                else:
                    hold_code = pmer.get_mmer_code(mmer_id)
                self.__motifs.append(list((m_type, hold_lbl, hold_code, pmer_id,
                                           pmer.get_mmer_center(mmer_id), pmer.get_mmer_rotation(mmer_id))))


class SetTomos:
    """
    Class for storing information for a set of simulated tomograms
    """

    def __init__(self):
        self.__tomos = list()

    def add_tomos(self, tomo):
        """
        Add a SynthTomo to the list
        :param tomo: input SyntTomo
        """
        assert isinstance(tomo, SynthTomo)
        self.__tomos.append(tomo)

    def save_csv(self, out_file):
        """
        Saves the motifs list contained in the set of tomograms
        :param out_file: output file path in .csv format
        """
        assert isinstance(out_file, str) and out_file.endswith('.csv')

        # Writing output CSV file
        with open(out_file, 'w') as csv_file:
            writer = csv.DictWriter(csv_file, dialect=csv.excel_tab, fieldnames=['Density', 'Micrographs', 'PolyData',
                           'Tomo3D', 'Type', 'Label', 'Code', 'Polymer', 'X', 'Y', 'Z', 'Q1', 'Q2', 'Q3', 'Q4'])
            writer.writeheader()

            # Tomos loop
            for tomo in self.__tomos:
                den_path = tomo.get_den()
                mics_path = tomo.get_mics()
                poly_path = tomo.get_poly()
                rec_path = tomo.get_tomo()

                # Motifs loop
                for motif in tomo.get_motif_list():
                    m_type = motif[0]
                    lbl_id = motif[1]
                    text_id = motif[2]
                    pmer_id = motif[3]
                    center = motif[4]
                    rotation = motif[5]

                    # Writing entry
                    writer.writerow({'Density':den_path, 'Micrographs':mics_path, 'PolyData':poly_path,
                                     'Tomo3D':rec_path, 'Type':m_type, 'Label':lbl_id, 'Code':text_id,
                                     'Polymer':pmer_id, 'X':center[0], 'Y':center[1], 'Z':center[2],
                                     'Q1':rotation[0], 'Q2':rotation[1], 'Q3':rotation[2], 'Q4':rotation[3]})