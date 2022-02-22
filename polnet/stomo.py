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
        if in_file is None:
            self.__mmer_id = None
            self.__mmer_svol = None
            self.__iso = None
            self.__pmer_l = None
            self.__pmer_occ = None
            self.__pmer_l_max = None
            self.__pmer_over_tol = 0
        else:
            self.load_protein_file(in_file)

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
                    var.strip(), value.strip()
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

    def add_networks(self, net, m_type, lbl, code):
        """
        Add al motifs within the input network to the synthetic tomogram
        :param net: a network object instance
        :param m_type: string with the type of monomers contained in the network
        :param lbl: integer label
        :param code: string code for the network monomers
        """
        assert issubclass(net, Network)
        assert isinstance(m_type, str) and isinstance(lbl, int) and isinstance(code, str)
        for pmer_id, pmer in enumerate(net.get_pmers_list()):
            for mmer_id in range(pmer.get_num_monomers()):
                self.__motifs.append(list((m_type, lbl, code, pmer_id,
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
        with open(out_file) as csv_file:
            csv.DictWriter(csv_file, dialect=csv.excel_tab, fieldnames=['Density', 'Micrographs', 'PolyData',
                           'Tomo3D', 'Type', 'Label', 'Code', 'Polymer', 'X', 'Y', 'Z', 'Q1', 'Q2', 'Q3', 'Q4'])

            # Tomos loop
            for tomo in self.__tomos:
                den_path = tomo.get_den()
                mics_path = tomo.get_mics()
                poly_path = tomo.get_poly()
                rec_path = tomo.get_tomo()

                # Motifs loop
                for motif in self.get_motif_list():
                    m_type = motif[0]
                    lbl_id = motif[1]
                    text_id = motif[2]
                    pmer_id = motif[3]
                    center = motif[4]
                    rotation = motif[5]

                    # Writing entry
                    file_entry = [den_path, mics_path, poly_path, rec_path, m_type, lbl_id, text_id, pmer_id, center[0],
                                  center[1], center[2], rotation[0], rotation[1], rotation[2], rotation[3]]
                    csv.DictWriter(file_entry, dialect=csv.excel_tab, fieldnames=['Density', 'Micrographs', 'PolyData',
                                'Tomo3D', 'Type', 'Label', 'Code', 'Polymer', 'X', 'Y', 'Z', 'Q1', 'Q2', 'Q3', 'Q4'])