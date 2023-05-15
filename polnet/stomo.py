"""
Classes with the info for a synthetic tomograms
"""

__author__ = 'Antonio Martinez-Sanchez'


import csv
from polnet.network import Network
from polnet.membrane import SetMembranes
from polnet import poly as pp


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
        self.__pmer_np = None
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

    def get_pmer_np(self):
        return self.__pmer_np

    def get_pmer_l_max(self):
        return self.__pmer_l_max

    def get_pmer_over_tol(self):
        return self.__pmer_over_tol

    def load_mmer_file(self, in_file):
        """
        Load protein parameters from an input file
        :param in_file: path to the input file with extension .pns
        """

        assert isinstance(in_file, str) and in_file.endswith('.pms') or in_file.endswith('.pns')

        # Reading input file
        with open(in_file) as file:
            for linea in file:
                if len(linea.strip()) > 0:

                    # Remove coments
                    linea = linea.split('#', 1)[0]

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
                        try:
                            self.__pmer_occ = float(value)
                        except ValueError:
                            value_0 = value[value.index('(') + 1:value.index(',')]
                            value_1 = value[value.index(',') + 1:value.index(')')]
                            self.__pmer_occ = (float(value_0), float(value_1))
                    elif var == 'PMER_NP':
                        self.__pmer_np = int(value)
                    elif var == 'PMER_L_MAX':
                        self.__pmer_l_max = float(value)
                    elif var == 'PMER_OVER_TOL':
                        self.__pmer_over_tol = float(value)
                    # else:
                    #     print('ERROR: (MmerFile - load_protein_file) input entry not recognized:', value)


class MmerMbFile(MmerFile):
    """
    For handling protein (or monomer) files
    """

    def __init__(self, in_file):
        """
        Constructor
        :param in_file: path to the input file with extension .pms
        """
        super().__init__(in_file)
        self.__mmer_center = None
        self.__mb_z_height = None
        self.__pmer_reverse_normals = None
        if in_file is not None:
            self.load_mmer_mb_file(in_file)

    def get_mmer_center(self):
        return self.__mmer_center

    def get_mb_z_height(self):
        return self.__mb_z_height

    def get_pmer_reverse_normals(self):
        return self.__pmer_reverse_normals

    def load_mmer_mb_file(self, in_file):
        """
        Load protein parameters from an input file
        :param in_file: path to the input file with extension .pms
        """

        super().load_mmer_file(in_file)

        # Reading input file
        with open(in_file) as file:
            for linea in file:
                if len(linea.strip()) > 0:

                    # Remove coments
                    linea = linea.split('#', 1)[0]

                    # Parsing an file entry
                    var, value = linea.split('=')
                    var = var.replace(' ', '')
                    var = var.replace('\n', '')
                    value = value.replace(' ', '')
                    value = value.replace('\n', '')
                    if var == 'MMER_CENTER':
                        hold_str = value.replace('[', '')
                        hold_str = hold_str.replace(']', '')
                        hold_str = hold_str.split(',')
                        self.__mmer_center = [float(hold_str[0]), float(hold_str[1]), float(hold_str[2])]
                    elif var == 'MB_Z_HEIGHT':
                        self.__mb_z_height = int(value)
                    elif var == 'PMER_REVERSE_NORMALS':
                        self.__pmer_reverse_normals = bool(value)




class MbFile:
    """
    For handling membrane configuration files
    """

    def __init__(self):
        self.__type = None
        self.__occ = None
        self.__thick_rg = None
        self.__layer_s_rg = None
        self.__max_ecc = None
        self.__over_tol = None
        self.__min_rad = None
        self.__den_cf_rg = None

    def get_type(self):
        return self.__type

    def get_occ(self):
        return self.__occ

    def get_thick_rg(self):
        return self.__thick_rg

    def get_layer_s_rg(self):
        return self.__layer_s_rg

    def get_max_ecc(self):
        return self.__max_ecc

    def get_over_tol(self):
        return self.__over_tol

    def get_min_rad(self):
        return self.__min_rad

    def get_den_cf_rg(self):
        return self.__den_cf_rg

    def load_mb_file(self, in_file):
        """
        Load protein parameters from an input file
        :param in_file: path to the input file with extension .mbs
        """

        assert isinstance(in_file, str) and in_file.endswith('.mbs')

        # Reading input file
        with open(in_file) as file:
            for linea in file:
                if len(linea.strip()) > 0:

                    # Remove coments
                    linea = linea.split('#', 1)[0]

                    # Parsing an file entry
                    var, value = linea.split('=')
                    var = var.replace(' ', '')
                    var = var.replace('\n', '')
                    value = value.replace(' ', '')
                    value = value.replace('\n', '')
                    if var == 'MB_TYPE':
                        self.__type = value
                    elif var == 'MB_OCC':
                        try:
                            self.__occ = float(value)
                        except ValueError:
                            value_0 = value[value.index('(') + 1:value.index(',')]
                            value_1 = value[value.index(',') + 1:value.index(')')]
                            self.__occ = (float(value_0), float(value_1))
                    elif var == 'MB_THICK_RG':
                        value_0 = value[value.index('(')+1:value.index(',')]
                        value_1 = value[value.index(',')+1:value.index(')')]
                        self.__thick_rg = (float(value_0), float(value_1))
                    elif var == 'MB_LAYER_S_RG':
                        value_0 = value[value.index('(') + 1:value.index(',')]
                        value_1 = value[value.index(',') + 1:value.index(')')]
                        self.__layer_s_rg = (float(value_0), float(value_1))
                    elif var == 'MB_MAX_ECC':
                        self.__max_ecc = float(value)
                    elif var == 'MB_OVER_TOL':
                        self.__over_tol = float(value)
                    elif var == 'MB_MIN_RAD':
                        self.__min_rad = float(value)
                    elif var == 'MB_DEN_CF_RG':
                        value_0 = value[value.index('(') + 1:value.index(',')]
                        value_1 = value[value.index(',') + 1:value.index(')')]
                        self.__den_cf_rg = (float(value_0), float(value_1))
                    else:
                        print('ERROR: (MmerFile - load_protein_file) input entry not recognized:', value)


class HelixFile:
    """
    For handling helicoidal structures configuration files
    """

    def __init__(self):
        self.__type = None
        self.__l = None
        self.__occ = None
        self.__min_p_len = None
        self.__mmer_rad = None
        self.__hp_len = None
        self.__mz_len = None
        self.__mz_len_f = None
        self.__thick_rg = None
        self.__over_tol = None
        self.__den_cf_rg = None
        self.__min_nmmer = None

    def get_type(self):
        return self.__type

    def get_l(self):
        return self.__l

    def get_occ(self):
        return self.__occ

    def get_thick_rg(self):
        return self.__thick_rg

    def get_min_p_len(self):
        return self.__min_p_len

    def get_mmer_rad(self):
        return self.__mmer_rad

    def get_hp_len(self):
        return self.__hp_len

    def get_mz_len(self):
        return self.__mz_len

    def get_mz_len_f(self):
        return self.__mz_len_f

    def get_over_tol(self):
        return self.__over_tol

    def get_den_cf_rg(self):
        return self.__den_cf_rg

    def get_min_nmmer(self):
        return self.__min_nmmer

    def load_hx_file(self, in_file):
        """
        Load protein parameters from an input file
        :param in_file: path to the input file with extension .hns
        """

        assert isinstance(in_file, str) and in_file.endswith('.hns')

        # Reading input file
        with open(in_file) as file:
            for linea in file:
                if len(linea.strip()) > 0:

                    # Remove coments
                    linea = linea.split('#', 1)[0]

                    # Parsing an file entry
                    var, value = linea.split('=')
                    var = var.replace(' ', '')
                    var = var.replace('\n', '')
                    value = value.replace(' ', '')
                    value = value.replace('\n', '')
                    if var == 'HLIX_TYPE':
                        self.__type = value
                    elif var == 'HLIX_PMER_L':
                        self.__l = float(value)
                    elif var == 'HLIX_PMER_OCC':
                        try:
                            self.__occ = float(value)
                        except ValueError:
                            value_0 = value[value.index('(') + 1:value.index(',')]
                            value_1 = value[value.index(',') + 1:value.index(')')]
                            self.__occ = (float(value_0), float(value_1))
                    elif var == 'HLIX_MIN_P_LEN':
                        self.__min_p_len = float(value)
                    elif var == 'HLIX_HP_LEN':
                        self.__hp_len = float(value)
                    elif var == 'HLIX_MZ_LEN':
                        self.__mz_len = float(value)
                    elif var == 'HLIX_MZ_LEN_F':
                        self.__mz_len_f = float(value)
                    elif var == 'HLIX_MMER_RAD':
                        self.__mmer_rad = float(value)
                    elif var == 'HLIX_OVER_TOL':
                        self.__over_tol = float(value)
                    elif var == 'HLIX_DEN_CF_RG':
                        value_0 = value[value.index('(') + 1:value.index(',')]
                        value_1 = value[value.index(',') + 1:value.index(')')]
                        self.__den_cf_rg = (float(value_0), float(value_1))
                    elif var == 'HLIX_MIN_NMMER':
                        self.__min_nmmer = int(value)


class MTFile(HelixFile):
    """
    For handling microtubular structures configuration files (inherits from Helix)
    """

    def __init__(self):
        super().__init__()
        self.__rad = None
        self.__nunits = None

    def get_rad(self):
        return self.__rad

    def get_nunits(self):
        return self.__nunits

    def load_mt_file(self, in_file):
        """
        Load protein parameters from an input file
        :param in_file: path to the input file with extension .mbs
        """

        super().load_hx_file(in_file)

        # Reading input file
        with open(in_file) as file:
            for linea in file:
                if len(linea.strip()) > 0:

                    # Remove coments
                    linea = linea.split('#', 1)[0]

                    # Parsing an file entry
                    var, value = linea.split('=')
                    var = var.replace(' ', '')
                    var = var.replace('\n', '')
                    value = value.replace(' ', '')
                    value = value.replace('\n', '')
                    if var == 'MT_RAD':
                        self.__rad = float(value)
                    elif var == 'MT_NUNITS':
                        self.__nunits = float(value)


class ActinFile(HelixFile):
    """
    For handling actin-like structures configuration files (inherits from Helix)
    """

    def __init__(self):
        super().__init__()
        self.__bprop = None
        self.__p_branch = None

    def get_bprop(self):
        return self.__bprop

    def get_p_branch(self):
        return self.__p_branch

    def load_ac_file(self, in_file):
        """
        Load protein parameters from an input file
        :param in_file: path to the input file with extension .mbs
        """

        super().load_hx_file(in_file)

        # Reading input file
        with open(in_file) as file:
            for linea in file:
                if len(linea.strip()) > 0:

                    # Remove coments
                    linea = linea.split('#', 1)[0]

                    # Parsing an file entry
                    var, value = linea.split('=')
                    var = var.replace(' ', '')
                    var = var.replace('\n', '')
                    value = value.replace(' ', '')
                    value = value.replace('\n', '')
                    if var == 'A_BPROP':
                        self.__bprop = float(value)
                    elif var == 'A_MAX_P_BRANCH':
                        self.__p_branch = float(value)



class SynthTomo:
    """
    This class model a Synthetic tomogram, tomogram's info are stored in disk to handle large datasets:
        - Ground truth: integer type tomogram associating each pixel with the id (background is always 0)
    """

    def __init__(self):
        self.__den = None
        self.__tomo = None
        self.__mics = None
        self.__poly = None
        self.__motifs = list()

    def get_den(self):
        return self.__den

    def get_mics(self):
        return self.__mics

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

    def add_network(self, net, m_type, lbl, code=None):
        """
        Add all motifs within the input network to the synthetic tomogram
        :param net: a network object instance
        :param m_type: string with the type of monomers contained in the network
        :param lbl: integer label
        :param code: string code for the network monomers, if None (default) it is taken from monomer information
        """
        assert issubclass(type(net), Network)
        assert isinstance(m_type, str)
        if (lbl is not None): assert isinstance(lbl, int)
        if (code is not None): assert isinstance(code, str)
        for pmer_id, pmer in enumerate(net.get_pmers_list()):
            for mmer_id in range(pmer.get_num_monomers()):
                if code is None:
                    hold_code = pmer.get_mmer_code(mmer_id)
                else:
                    hold_code = code
                self.__motifs.append(list((m_type, lbl, hold_code, pmer_id,
                                           pmer.get_mmer_center(mmer_id), pmer.get_mmer_rotation(mmer_id))))

    def add_set_mbs(self, set_mbs, m_type, lbl, code, dec=None):
        """
        Membrane surface point coordinates are added to the tomogram motif list
        In rotations the normal vector to each point is stored as: X->Q0, Y->Q1 , Z->Q2 and 0->Q3
        :param set_mbs: a membrane set object instance
        :param m_type: string with the type of motif contained in the network
        :param lbl: integer label
        :param code: string code for membrane
        :param dec: if not None (default) the membrane points are decimated according this factor
        """
        assert issubclass(type(set_mbs), SetMembranes)
        assert isinstance(m_type, str)
        assert isinstance(lbl, int)
        assert isinstance(code, str)

        poly_vtp = set_mbs.get_vtp()
        if dec is not None:
            poly_vtp = pp.poly_decimate(poly_vtp, dec)

        n_points = poly_vtp.GetNumberOfPoints()
        normals = poly_vtp.GetPointData().GetNormals()
        for i in range(n_points):
            x, y, z = poly_vtp.GetPoint(i)
            q0, q1, q2 = normals.GetTuple(i)
            self.__motifs.append(list((m_type, lbl, code, i, [x, y, z], [q0, q1, q2, 0])))


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