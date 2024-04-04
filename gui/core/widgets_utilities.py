import ipywidgets as widgets
from IPython.display import display, HTML
from ipyfilechooser import FileChooser
import os


def widgets_visualize_pdb():
    """
    Widgets to download PDB from protein data bank
    return: tuple with all widgets created
    """
    text_protein_widget = widgets.Text(description="Enter protein name:", value="1bxn",placeholder = "1bxn")
    text_protein_widget.style.description_width = 'initial'
    text_protein_widget.layout.width = '450px'

    save_button = widgets.Button(description="Save")
    file_destination_widget = FileChooser(os.getcwd(), title = "Select the path and tha name for the file")
    
    display(widgets.HBox([text_protein_widget,save_button]), file_destination_widget)
    return text_protein_widget,save_button,file_destination_widget
    
    
def widgets_convert_to_mrc():
    """
    Widgets to exec the PDB to MRC converter
    return: tuple with all widgets created
    """

    custom_width = '550px'

    title = widgets.Label(value="Introduce a the files you want to create the density map")
    files_button = widgets.Button(description="Select files you want to convert")
    files_button.layout.width= '200px' 
    
    dropdown = widgets.Dropdown(value=None, layout={'width': custom_width})
    
    apix_widget = widgets.BoundedFloatText(value=10, min=0, step=0.001, description="VSIZE:")
    apix_widget.style.description_width = 'initial'

    res_widget = widgets.BoundedFloatText(value=30, min=0, step=1, description="RESOLUTION:")
    res_widget.style.description_width = 'initial'

    offset_widget = widgets.BoundedFloatText(value=20, min=0, step=1, description="OFFSET:")
    offset_widget.style.description_width = 'initial'
    
    het_widget = widgets.Checkbox(value=False, description="Include HETATM atoms")
    het_widget.style.description_width = 'initial' 
    
    chains_widget = widgets.Text(value="", description="Protein Chains:")
    chains_widget.style.description_width = 'initial' 
    
    use_model_checkbox = widgets.Checkbox(value=False, description="Use Model Number:")
    use_model_checkbox.style.description_width = 'initial' 
    
    model_widget = widgets.BoundedIntText(value=None, min=1, step=1, description="Model Number:", disabled=True)
    model_widget.style.description_width = 'initial' 
    
    widgets.jslink((use_model_checkbox, 'value'), (model_widget, 'disabled'))
    
    file_destination_widget = FileChooser(os.getcwd(), show_only_dirs=True, title="Select where you want to save the files:")
    
    exec_button = widgets.Button(description="Convert")
    
    files_layaout = widgets.HBox([files_button, dropdown])
    full_layout = widgets.VBox([title, files_layaout, apix_widget, res_widget, offset_widget, het_widget, chains_widget, use_model_checkbox, model_widget, file_destination_widget, exec_button])
    
    display(full_layout)
    
    return apix_widget, res_widget, offset_widget, het_widget, chains_widget, use_model_checkbox, model_widget, files_button, dropdown, file_destination_widget, exec_button


def widgets_mmolecules_params():
    """
    Widgets to select mmolecules params
    :return: tuple with all widgets created
    """

    custom_width = '550px'
    title = widgets.Label(value="Introduce the protein you want to create the macromolecule")
    file_destination_widget = FileChooser(os.getcwd(), title  = "Select where you want to save the files")
    
    protein_file_button = widgets.Button(description="Select protein file")
    dropdown_1 = widgets.Dropdown(value = None, layout={'width': custom_width})
    
    
    save_button = widgets.Button(description="Save")
    
    # Crear widgets
    mmer_id_widget = widgets.Text(description='MMER_ID:', value='pdb_1bxn', placeholder="pdb_1bxn")
    mmer_id_widget.style.description_width = 'initial' 
    
    pmer_l_widget = widgets.BoundedFloatText(value=1.2, min=0, step=0.0001, description="PMER_L:")
    pmer_l_widget.style.description_width = 'initial' 
    
    pmer_l_mex_widget = widgets.BoundedFloatText(value=3000, min=0, step=0.0001, max=100000000, description="PMER_LMEX:")
    pmer_l_mex_widget.style.description_width = 'initial' 
    
    pmer_occ_widget = widgets.BoundedFloatText(value=0.5, min=0, step=0.0001, max=1, description="PMER_OCC:")
    pmer_occ_widget.style.description_width = 'initial' 
    
    pmer_over_tol_widget = widgets.BoundedFloatText(value=0.01, min=0, step=0.0001, max=1, description="PMER_OVER_TOL:")
    pmer_over_tol_widget.style.description_width = 'initial' 

    pmer_reverse_normals = widgets.Checkbox(value=False, disabled=True, indent=False)
    pmer_label = widgets.Label("Marks if you want to reverse membrane protein normals")
    membrane_option = widgets.Checkbox(value=False, disabled=False, indent=False)
    membrane_label = widgets.Label("Mark if it is a membrane protein")

    def update_reverse_normals_checkbox(change):
        pmer_reverse_normals.disabled = not change.new

    membrane_option.observe(update_reverse_normals_checkbox, names='value')

    widgets_top = widgets.HBox([protein_file_button, dropdown_1])
    widgets_row = widgets.HBox([membrane_label, membrane_option])
    widgets_row2 = widgets.HBox([pmer_label, pmer_reverse_normals])
    
    labels = {
        'mmer_id_widget': "String with macromolecules identifier",
        'pmer_l_widget': "Positve real value whit the variable distance between two consecutive macromolecules in a cluster",
        'pmer_l_mex_widget': "Positive real value in Ángstrom whit cluster maximum length",
        'pmer_occ_widget': "Percentage occupancy",
        'pmer_over_tol_widget': "Percentage overlapping tolerance",
    }
    
    widgets_and_labels = [
        widgets.HBox([mmer_id_widget, widgets.Label(value=labels['mmer_id_widget'])]),
        widgets.HBox([pmer_l_widget, widgets.Label(value=labels['pmer_l_widget'])]),
        widgets.HBox([pmer_l_mex_widget, widgets.Label(value=labels['pmer_l_mex_widget'])]),
        widgets.HBox([pmer_occ_widget, widgets.Label(value=labels['pmer_occ_widget'])]),
        widgets.HBox([pmer_over_tol_widget, widgets.Label(value=labels['pmer_over_tol_widget'])]),
    ]
    display(title, widgets_top, widgets_row, widgets.VBox(widgets_and_labels), widgets_row2, file_destination_widget, save_button)
    return protein_file_button, dropdown_1, save_button, mmer_id_widget, pmer_l_widget, pmer_l_mex_widget, pmer_occ_widget, pmer_over_tol_widget, pmer_reverse_normals, membrane_option, file_destination_widget


def widgets_create_axis():
    """
    Widgets to select mmolecules params
    :return: tuple with all widgets created
    """
    custom_width = '550px'
    
    title = widgets.Label(value="Introduce a the files you want to create an axis")
    select_files_button = widgets.Button(description="Select Files", layout=widgets.Layout(margin='0 5px 0 0'))  
    dropdown_1 = widgets.Dropdown(value = None, layout={'width': custom_width})

    file_destination_widget = FileChooser(os.getcwd(), show_only_dirs=True, title="Select where you want to save the files")
    
    mb_size_widget = widgets.BoundedFloatText(value=50, min=0, step=0.001, max=100000000, description="MB size:")
    mb_size_widget.style.description_width = 'initial' 
    
    zaxis_rad_widget = widgets.BoundedFloatText(value=5, min=0, step=0.001, description="Zaxis rad:")
    zaxis_rad_widget.style.description_width = 'initial' 
    
    vsize_widget = widgets.BoundedFloatText(value=10, min=0, step=0.001, description="Voxel size:")
    vsize_widget.style.description_width = 'initial' 
    
    save_button = widgets.Button(description="Save", layout=widgets.Layout(margin='8px 0 0 0'))

    scale_factor_widget = widgets.BoundedFloatText(value=1.5, min=1, step=0.1, max=3, description="Scale factor:")
    scale_factor_widget.style.description_width = 'initial' 
    
    labels = {
        'mb_size_widget': "Positive integer value in Angstroms ",
        'zaxis_rad_widget': "Positive integer value",
        'vsize_widget': "Positive integer value",
        'scale_factor_widget': "Scale axis size"
    }
    
    widgets_and_labels = [
        widgets.HBox([mb_size_widget, widgets.Label(value=labels['mb_size_widget'])]),
        widgets.HBox([zaxis_rad_widget, widgets.Label(value=labels['zaxis_rad_widget'])]),
        widgets.HBox([vsize_widget, widgets.Label(value=labels['vsize_widget'])]),
        widgets.HBox([scale_factor_widget, widgets.Label(value=labels['scale_factor_widget'])]),
    ]

    display(title, widgets.HBox([select_files_button,dropdown_1]),  widgets.VBox(widgets_and_labels), file_destination_widget, save_button)
    return mb_size_widget, zaxis_rad_widget, vsize_widget, scale_factor_widget, select_files_button, dropdown_1, save_button, file_destination_widget


def widgets_align():
    """
    Widgets to select molecules params
    :return: tuple with all widgets created
    """

    custom_width = '550px'
   
    title = widgets.Label(value="Introduce a reference membrane membrane and a protein:")
    upload_protein_button = widgets.Button(description="Select protein file")
    dropdown_1 = widgets.Dropdown(value = None, layout={'width': custom_width})
    upload_axis_button = widgets.Button(description="Select membrane file")
    dropdown_2 = widgets.Dropdown(value = None, layout={'width': custom_width})
    exec_button = widgets.Button(description="Exec")
    vsize_widget = widgets.BoundedFloatText(value=10, min=0, step=0.001, description="Voxel size:")
    vsize_widget.style.description_width = 'initial' 
    
    file_destination_widget = FileChooser(os.getcwd(), title  = "Select where you want to save the files")
    
    display(title, widgets.HBox([upload_protein_button, dropdown_1]),widgets.HBox([upload_axis_button, dropdown_2]), file_destination_widget, vsize_widget, exec_button)

    return upload_protein_button, upload_axis_button, exec_button, dropdown_1, dropdown_2, vsize_widget,file_destination_widget


def widgets_membrane_params():
    """
    Create membrane file
    :return: tuple with all widgets created
    """
    # Types
    options_mb_type = ['sphere', 'ellipse', 'torus']
    
    # Definir los widgets para cada variable
    mb_type_widget = widgets.Dropdown(options=options_mb_type, value='sphere', description='MB_TYPE:')
    mb_type_widget.style.description_width = 'initial' 
    
    mb_occ_widget = widgets.FloatRangeSlider(value=[0.001, 0.003], min=0 ,max=100, step=0.0001, readout_format='.4f',description='MB_OCC:')
    mb_occ_widget.style.description_width = 'initial'
    mb_occ_widget.layout.width = '350px'
    
    mb_thick_rg_widget = widgets.FloatRangeSlider(value=[25, 45], min=20 ,max=50, step=1, readout_format='.3f', description='MB_THICK_RG:')
    mb_thick_rg_widget.style.description_width = 'initial'
    mb_thick_rg_widget.layout.width = '350px'
    
    mb_layer_s_rg_widget = widgets.FloatRangeSlider(value=[0.5, 2], min=0.1 ,max=3,step=0.1, readout_format='.3f', description='MB_LAYER_S_RG:')
    mb_layer_s_rg_widget.style.description_width = 'initial'
    mb_layer_s_rg_widget.layout.width = '350px'
    
    mb_max_ecc_widget = widgets.BoundedFloatText(value=0.75 ,min=0,step=0.0001, description='MB_MAX_ECC:')
    mb_max_ecc_widget.style.description_width = 'initial'
    
    mb_over_tol_widget = widgets.BoundedFloatText(value=1e-9, min = 0, step = 0.00000000001, max = 1, description='MB_OVER_TOL:')
    mb_over_tol_widget.style.description_width = 'initial'
    
    mb_min_rad_widget = widgets.BoundedFloatText(value=75 ,min=0,step=0.0001,  max= 100000000, description='MB_MIN_RAD:')
    mb_min_rad_widget.style.description_width = 'initial'
    
    mb_den_cf_rg_widget = widgets.FloatRangeSlider(value=[0.3,0.5], min=0, max=1, step=0.0001,description='MB_DEN_CF_RG:')
    mb_den_cf_rg_widget.style.description_width = 'initial'

    file_destination_widget = FileChooser(os.getcwd(), title="Select destination path")

    save_button = widgets.Button(description="Save")
    
    labels = {
    'mb_type_widget': 'Membrane geometry: sphere, ellipse or torus',
    'mb_occ_widget': '% occupancy',
    'mb_thick_rg_widget': 'Membrane thickness or distance between both layers in Angstroms',
    'mb_layer_s_rg_widget': 'Membrane layer thickness in Angstroms',
    'mb_max_ecc_widget': 'Maximum ellipsoid eccentricity',
    'mb_over_tol_widget': '% overlapping tolerance',
    'mb_min_rad_widget': 'Minimum spherical membrane radius',
    'mb_den_cf_rg_widget': '% density factor',
    }
    
    widgets_and_labels = [
    widgets.HBox([mb_type_widget,widgets.Label(value=labels['mb_type_widget'])]),
    widgets.HBox([mb_occ_widget, widgets.Label(value=labels['mb_occ_widget'])]),
    widgets.HBox([mb_thick_rg_widget, widgets.Label(value=labels['mb_thick_rg_widget'])]),
    widgets.HBox([mb_layer_s_rg_widget, widgets.Label(value=labels['mb_layer_s_rg_widget'])]),
    widgets.HBox([mb_max_ecc_widget, widgets.Label(value=labels['mb_max_ecc_widget'])]),
    widgets.HBox([mb_over_tol_widget, widgets.Label(value=labels['mb_over_tol_widget'])]),
    widgets.HBox([mb_min_rad_widget, widgets.Label(value=labels['mb_min_rad_widget'])]),
    widgets.HBox([mb_den_cf_rg_widget, widgets.Label(value=labels['mb_den_cf_rg_widget'])]),
    ]
    
    display(widgets.VBox(widgets_and_labels), file_destination_widget, save_button)
    return mb_type_widget, mb_occ_widget, mb_thick_rg_widget, mb_layer_s_rg_widget, mb_max_ecc_widget, mb_over_tol_widget, mb_min_rad_widget, mb_den_cf_rg_widget, file_destination_widget, save_button


def widgets_helix_params():
    """
    Create membrane file
    :return: tuple with all widgets created
    """
    #Style
    options_hlix_type = ['actin', 'mt']
    
    hlix_type_widget = widgets.Dropdown(options=options_hlix_type, value='actin', description='HLIX_TYPE:')
    hlix_type_widget.style.description_width = 'initial'
    
    hlix_mmer_rad_widget = widgets.BoundedFloatText(value=25, min=0, step=0.0001,max= 10000, description='HLIX _MMER_RAD:')
    hlix_mmer_rad_widget.style.description_width = 'initial'
    
    hlix_pmer_l_widget = widgets.BoundedFloatText(value=1.2, min=0, step=0.0001, description='HLIX_PMER_L:')
    hlix_pmer_l_widget.style.description_width = 'initial'
    
    hlix_pmer_occ_widget =  widgets.FloatRangeSlider(value=[0.001, 0.003], min=0 ,max=1, step=0.0001, readout_format='.4f',description='MB_OCC:')
    hlix_pmer_occ_widget.style.description_width = 'initial'
    
    hlix_min_p_len_widget = widgets.BoundedFloatText(value=177000, min=0, step=0.0001, max= 1000000000000, description='HLIX_MIN_P_LEN:')
    hlix_min_p_len_widget.style.description_width = 'initial'
    
    hlix_hp_len_widget = widgets.BoundedFloatText(value=720, min=0, step=0.0001, max= 100000000000, description='HLIX_HP_LEN:')
    hlix_hp_len_widget.style.description_width = 'initial'
    
    hlix_mz_len_widget = widgets.BoundedFloatText(value=50, min=0, step=0.0001, description='HLIX_MZ_LEN:')
    hlix_mz_len_widget.style.description_width = 'initial'
    
    hlix_mz_len_f_widget = widgets.BoundedFloatText(value=0.2, min=0, step=0.0001, max = 1, description='HLIX_MZ_LEN_F:')
    hlix_mz_len_f_widget.style.description_width = 'initial'
    
    hlix_over_tol_widget = widgets.BoundedFloatText(value=1e-9, min=0, step=0.00000000001, max = 1, description='HLIX_OVER_TOL:')
    hlix_over_tol_widget.style.description_width = 'initial'
    
    hlix_min_nmmer_widget = widgets.BoundedIntText(value=15, min = 0, max = 1000000000, description='HLIX_MIN_NMMER:')
    hlix_min_nmmer_widget.style.description_width = 'initial'
    
    a_bprop_widget = widgets.BoundedFloatText(value=0.5, max=1, min=0, description='A_BPROP:')
    a_bprop_widget.style.description_width = 'initial'
    
    a_max_p_branch_widget = widgets.BoundedFloatText(value=5, max=100000000, min=0, description='A_MAX_P_BRANCH:')
    a_max_p_branch_widget.style.description_width = 'initial'
    
    mt_rad_widget = widgets.BoundedFloatText(value=100.5, min = 0,  max = 100000000,description='MT_RAD:')
    mt_rad_widget.style.description_width = 'initial'
    mt_rad_widget.disabled = True
    
    mt_nunits_widget = widgets.BoundedIntText(value=13, max= 100000, min = 0,description='MT_NUNITS:')
    mt_nunits_widget.style.description_width = 'initial'
    mt_nunits_widget.disabled = True

    file_destination_widget = FileChooser(os.getcwd(), title = "Select the path and the name for the new file")
    
    save_button = widgets.Button(description="Save")
    
    # Etiquetas para los widgets
    labels = {
        'hlix_type_widget': 'Filament type: actin or microtubule',
        'hlix_mmer_rad_widget': 'monomer radius in Angstroms',
        'hlix_pmer_l_widget': 'Distance in Angstroms between two structural units in a polymer / monomer radius',
        'hlix_pmer_occ_widget': '% occupancy',
        'hlix_min_p_len_widget': 'Minimum polymer persistence length',
        'hlix_hp_len_widget': 'Length in Angstroms to complete an inner turn of structural units on a center-line curve',
        'hlix_mz_len_widget': 'monomer length in Angstroms in z-axis',
        'hlix_mz_len_f_widget': 'Minimum helical curve elevation in Z-axis as a fraction of circular component',
        'hlix_over_tol_widget': '% overlapping tolerance',
        'hlix_min_nmmer_widget': 'Minimum number of structural units, positive integer',
        'a_bprop_widget': 'Branching probability for filament networks',
        'a_max_p_branch_widget': 'Maximum branching probability for filament networks',
        'mt_rad_widget': 'Microtubule radius Z-ring in Angstroms',
        'mt_nunits_widget': 'Microtubule number of units por ring, positive integer'
    }
    
    widgets_and_labels = [
        widgets.HBox([hlix_type_widget, widgets.Label(value=labels['hlix_type_widget'])]),
        widgets.HBox([hlix_mmer_rad_widget, widgets.Label(value=labels['hlix_mmer_rad_widget'])]),
        widgets.HBox([hlix_pmer_l_widget, widgets.Label(value=labels['hlix_pmer_l_widget'])]),
        widgets.HBox([hlix_pmer_occ_widget, widgets.Label(value=labels['hlix_pmer_occ_widget'])]),
        widgets.HBox([hlix_min_p_len_widget, widgets.Label(value=labels['hlix_min_p_len_widget'])]),
        widgets.HBox([hlix_hp_len_widget, widgets.Label(value=labels['hlix_hp_len_widget'])]),
        widgets.HBox([hlix_mz_len_widget, widgets.Label(value=labels['hlix_mz_len_widget'])]),
        widgets.HBox([hlix_mz_len_f_widget, widgets.Label(value=labels['hlix_mz_len_f_widget'])]),
        widgets.HBox([hlix_over_tol_widget, widgets.Label(value=labels['hlix_over_tol_widget'])]),
        widgets.HBox([hlix_min_nmmer_widget, widgets.Label(value=labels['hlix_min_nmmer_widget'])]),
        widgets.HBox([a_bprop_widget, widgets.Label(value=labels['a_bprop_widget'])]),
        widgets.HBox([a_max_p_branch_widget, widgets.Label(value=labels['a_max_p_branch_widget'])]),
        widgets.HBox([mt_rad_widget, widgets.Label(value=labels['mt_rad_widget'])]),
        widgets.HBox([mt_nunits_widget, widgets.Label(value=labels['mt_nunits_widget'])])
    ]
    
    # Mostrar los contenedores HBox en un VBox
    display(widgets.VBox(widgets_and_labels), file_destination_widget, save_button)
    return hlix_type_widget, hlix_mmer_rad_widget, hlix_pmer_l_widget, hlix_pmer_occ_widget, hlix_min_p_len_widget, hlix_hp_len_widget, hlix_mz_len_widget, hlix_mz_len_f_widget, hlix_over_tol_widget, hlix_min_nmmer_widget, a_bprop_widget, a_max_p_branch_widget, mt_rad_widget, mt_nunits_widget, file_destination_widget, save_button

    
def update_widgets_based_on_hlix(hlix_type, a_bprop_widget, a_max_p_branch_widget, mt_rad_widget, mt_nunits_widget):
    """
    Update helix widgets
    :param hlix_type: filament type
    :param a_bprop_widget: widget to update
    :param a_max_p_branch_widget: widget to update
    :param mt_rad_widget: widget to update
    :param mt_nunits_widget:  widget to update
    """
    if hlix_type == 'actin':
        a_bprop_widget.disabled = False
        a_max_p_branch_widget.disabled = False
        mt_rad_widget.disabled = True
        mt_nunits_widget.disabled = True
    elif hlix_type == 'mt':
        a_bprop_widget.disabled = True
        a_max_p_branch_widget.disabled = True
        mt_rad_widget.disabled = False
        mt_nunits_widget.disabled = False



def update_values_based_on_hlix(hlix_type, hlix_mmer_rad, hlix_pmer_l, hlix_pmer_occ, hlix_min_p_len, hlix_hp_len, hlix_mz_len,
                                hlix_mz_len_f, hlix_over_tol, hlix_min_nmmer, a_bprop, a_max_p_branch, mt_rad, mt_nunits):
    """
    Update helix values based on filament type.
    :param hlix_type: Filament type ('actin' or 'mt').
    :param hlix_mmer_rad: Monomer radius.
    :param hlix_pmer_l: Distance between two structural units in a polymer / monomer radius.
    :param hlix_pmer_occ: Percentage occupancy as a range [min, max].
    :param hlix_min_p_len: Minimum polymer persistence length.
    :param hlix_hp_len: Length to complete an inner turn of structural units on a center-line curve.
    :param hlix_mz_len: Monomer length in z-axis.
    :param hlix_mz_len_f: Minimum helical curve elevation in Z-axis as a fraction of circular component.
    :param hlix_over_tol: Overlapping tolerance percentage.
    :param hlix_min_nmmer: Minimum structural units (positive integer).
    :param a_bprop: Branching probability for filament networks.
    :param a_max_p_branch: Maximum branching probability for filament networks.
    :param mt_rad: Microtubule radius Z-ring.
    :param mt_nunits: Microtubule ring number of monomers (positive integer).
    """
    if hlix_type == 'actin':
        hlix_mmer_rad.value = 25
        hlix_pmer_l.value = 1.2
        hlix_pmer_occ.value = [0.001, 0.003]
        hlix_min_p_len.value = 177000
        hlix_hp_len.value = 720
        hlix_mz_len.value = 50
        hlix_mz_len_f.value = 0.2
        hlix_over_tol.value = 1e-9
        hlix_min_nmmer.value = 15
        a_bprop.value = 0.5
        a_max_p_branch.value = 5
        mt_rad.value = 100.5
        mt_nunits.value = 13
    else:
        hlix_mmer_rad.value = 40
        hlix_pmer_l.value = 1.2
        hlix_pmer_occ.value = [0.3, 1]
        hlix_min_p_len.value = 5200e4 
        hlix_hp_len.value = 216000 
        hlix_mz_len.value = 50
        hlix_mz_len_f.value = 0.4
        hlix_over_tol.value = 1e-9
        hlix_min_nmmer.value = 15
        a_bprop.value = 0.5
        a_max_p_branch.value = 5
        mt_rad.value = 100.5
        mt_nunits.value = 13


def widgets_helix():
    """
    Widgets to visualize filament
    :return: tuple with widget created
    """
    custom_width = '550px'
    title = widgets.Label(value='Select file you want to visualize')
    select_file_button = widgets.Button(description="Select file")
    dropdown_1 = widgets.Dropdown(value=None, layout={'width': custom_width})
    
    exec_button = widgets.Button(description="Exec")
    display(widgets.VBox([title, widgets.HBox([select_file_button, dropdown_1]),  exec_button]))

    return select_file_button, exec_button, dropdown_1 


def widgets_add_app_files():
    """widgets to add membrane files
    :return: tuple with all widgets created
    """
    # Widgets
    custom_width = '550px'
    
    label = widgets.Label(value='Types of structures:')
    
    select_file_button_membrane = widgets.Button(description="Select membrane files")
    select_file_button_membrane.layout.width = '200px'
    dropdown_membrane = widgets.Dropdown(options=['Select paths...'], layout={'width': custom_width})
    hbox_membrane = widgets.HBox([select_file_button_membrane, dropdown_membrane])

    select_file_button_helix = widgets.Button(description="Select filament files")
    select_file_button_helix.layout.width = '200px'
    dropdown_helix = widgets.Dropdown(options=['Select paths...'], layout={'width': custom_width})
    hbox_helix = widgets.HBox([select_file_button_helix, dropdown_helix])

    select_file_button_proteins = widgets.Button(description="Select proteins files")
    select_file_button_proteins.layout.width = '200px'
    dropdown_proteins = widgets.Dropdown(options=['Select paths...'], layout={'width': custom_width})
    hbox_proteins = widgets.HBox([select_file_button_proteins, dropdown_proteins])

    select_file_button_mproteins = widgets.Button(description="Select membrane protein files")
    select_file_button_mproteins.layout.width = '200px'
    dropdown_mproteins = widgets.Dropdown(options=['Select paths...'], layout={'width': custom_width})
    #dropdown_mproteins = widgets.Dropdown( value=None, layout={'width': custom_width})
    hbox_mproteins = widgets.HBox([select_file_button_mproteins, dropdown_mproteins])


    # Display
    display(widgets.VBox([label, hbox_membrane, hbox_helix, hbox_proteins, hbox_mproteins]))

    return select_file_button_membrane, dropdown_membrane, select_file_button_helix, dropdown_helix, select_file_button_proteins, dropdown_proteins, select_file_button_mproteins, dropdown_mproteins


def widgets_show_files(list1, list2, list3, list4):
    """
    Widgets to show files with Dropdowns.
    :param list1: List of options for Dropdown 1.
    :param list1: List of options for Dropdown 2.
    :param list1: List of options for Dropdown 3.
    :param list1: List of options for Dropdown 4.
    
    """
    description_width = 'initial'  # Puedes ajustar este valor según tus necesidades
    custom_width = '550px'
    
    dropdown_1 = widgets.Dropdown(options=list1, description='Membrane files:', style={'description_width': description_width}, layout={'width': custom_width})
    dropdown_2 = widgets.Dropdown(options=list2, description='Helix files:', style={'description_width': description_width}, layout={'width': custom_width})
    dropdown_3 = widgets.Dropdown(options=list3, description='Protein files:', style={'description_width': description_width}, layout={'width': custom_width})
    dropdown_4 = widgets.Dropdown(options=list4, description='Membrane protein files:', style={'description_width': description_width}, layout={'width': custom_width})

    
    display(dropdown_1, dropdown_2, dropdown_3, dropdown_4)
    return dropdown_1, dropdown_2, dropdown_3, dropdown_4


def widgets_exec_app():
    """
    Widgets to exe the app
    :return: tuple with widget created
    """
    widget_out_dir = FileChooser(os.getcwd(), show_only_dirs=True, title="Select where you want to save the output files:")
    ntomos_widget = widgets.IntText(value=1,description='N_TOMOS (number of tomograms in the dataset):')
    ntomos_widget.style.description_width = 'initial'
    ntomos_widget.layout.width = '380px'

    voi_shape1 = widgets.BoundedIntText(value=400, min = 0, max= 100000, description='VOI_SHAPE (Tomogram shape, voxels):')
    voi_shape1.style.description_width = 'initial'
    voi_shape1.layout.width = '310px'
    voi_shape2 = widgets.BoundedIntText(value=400, min = 0, max= 100000 )
    voi_shape2.layout.width = '80px'
    voi_shape3 = widgets.BoundedIntText(value=100, min = 0, max= 100000 )
    voi_shape3.layout.width = '80px'

    hbox_voi_shape = widgets.HBox([voi_shape1, voi_shape2, voi_shape3])
    
    voi_size_widget = widgets.BoundedFloatText(value=10, min=0, description='VOI_VOXEL_SIZE (Voxel size, voxels/A):')
    voi_size_widget.style.description_width = 'initial'
    voi_size_widget.layout.width = '350px'
    
    mmer_tries_widget = widgets.BoundedFloatText(value=20, min=0, description='MMER_TRIES (Maximun number of tries for monomers):')
    mmer_tries_widget.style.description_width = 'initial'
    mmer_tries_widget.layout.width = '400px'
    
    pmer_tries_widget = widgets.BoundedFloatText(value=100, min=0, max= 1000, description='PMER_TRIES (Maximun number of tries for polymers):')
    pmer_tries_widget.style.description_width = 'initial'
    pmer_tries_widget.layout.width = '390px'

    surf_dec_widget = widgets.BoundedFloatText(value=0.9, min=0, description='SURF_DEC (Decimation for surface respresentation, [0, 1]):')
    surf_dec_widget.style.description_width = 'initial'
    surf_dec_widget.layout.width = '410px'

    malign_mn_widget = widgets.BoundedFloatText(value=1, min=0, description='MALIGN_MN (Micrograph miss-alginment mean, pixels):')
    malign_mn_widget.style.description_width = 'initial'
    malign_mn_widget.layout.width = '400px'

    malign_mx_widget = widgets.BoundedFloatText(value=1.5, min=0, description='MALIGN_MX (Micrograph miss-alginment max, pixels):')
    malign_mx_widget.style.description_width = 'initial'
    malign_mx_widget.layout.width = '400px'

    malign_sg_widget = widgets.BoundedFloatText(value=0.2, min=0, description='MALIGN_SG (Micrograph miss-alginment sigma, pixels):')
    malign_sg_widget.style.description_width = 'initial'
    malign_sg_widget.layout.width = '400px'

    detector_snr_widget_low = widgets.BoundedFloatText(value=1, min=0, step=0.0001, description='DETECTOR_SNR (Micrographs SNR range):')
    detector_snr_widget_low.style.description_width = 'initial'
    detector_snr_widget_low.layout.width = '330px'
    detector_snr_widget_high = widgets.BoundedFloatText(value=2, min=0, step=0.0001)
    detector_snr_widget_high.style.description_width = 'initial'
    detector_snr_widget_high.layout.width = '80px'

    hbox_snr = widgets.HBox([detector_snr_widget_low, detector_snr_widget_high])

    widget_min = widgets.BoundedIntText(value=-60, min = -1000, max= 100000, description='TILT_ANGS (Degrees; start, end, step):')
    widget_min.style.description_width = 'initial'
    widget_min.layout.width = '290px'
    widget_max = widgets.BoundedIntText(value=61, min = -1000, max= 100000 )
    widget_max.layout.width = '60px'
    widget_paso = widgets.BoundedIntText(value=3, min = -1000,max= 100000 )
    widget_paso.layout.width = '60px'

    hbox_tilt_angs= widgets.HBox([widget_min, widget_max, widget_paso])
  

    voi_off_widget_1 = widgets.BoundedIntText(value=4, description='VOI_OFF (Empty halo, voxels):')
    voi_off_label = widgets.Label('this is the Start and End in voxels for the effective volume')
    voi_off_widget_1.style.description_width = 'initial'
    voi_off_widget_1.layout.width = '250px'
    voi_off_widget_2 = widgets.BoundedIntText(value=396, max = 100000, min = 1)
    voi_off_widget_2.layout.width = '70px'
    voi_off_widget_3 = widgets.BoundedIntText(value=4, max = 100000, min = 1 )
    voi_off_widget_3.layout.width = '70px'
    voi_off_widget_4 = widgets.BoundedIntText(value=396,max = 100000, min = 1)
    voi_off_widget_4.layout.width = '70px'
    voi_off_widget_5 = widgets.BoundedIntText(value=4, max = 100000, min = 1 )
    voi_off_widget_5.layout.width = '70px'
    voi_off_widget_6 = widgets.BoundedIntText(value=96, max = 100000, min = 1)
    voi_off_widget_6.layout.width = '70px'

    exec_button = widgets.Button(description="Exec")

    hbox_voi_off= widgets.HBox([voi_off_widget_1,voi_off_widget_2,voi_off_widget_3,voi_off_widget_4,voi_off_widget_5,voi_off_widget_6, voi_off_label])
    
    display(widget_out_dir, ntomos_widget, hbox_voi_shape, hbox_voi_off, voi_size_widget,  mmer_tries_widget, pmer_tries_widget,
             surf_dec_widget, malign_mn_widget, malign_mx_widget, malign_sg_widget,
            hbox_snr, hbox_tilt_angs, exec_button)

    return widget_out_dir, ntomos_widget, voi_shape1, voi_shape2, voi_shape3, voi_off_widget_1, voi_off_widget_2, voi_off_widget_3, voi_off_widget_4, \
           voi_off_widget_5, voi_off_widget_6,voi_size_widget,mmer_tries_widget, \
           pmer_tries_widget, surf_dec_widget, malign_mn_widget, malign_mx_widget, malign_sg_widget,  \
           detector_snr_widget_low, detector_snr_widget_high, widget_min, widget_max, widget_paso, exec_button


def widgets_change_order(lists):
    """
    Change files order
    :return: tuple with all widgets created
    """
    widgets_list = []
    for list in lists:
        if len(list) > 1:
            dropdown_widget = widgets.Dropdown(options=list, value=list[0], layout={'width': '550px'})
        else:
             dropdown_widget = widgets.Dropdown(options=list, layout={'width': '550px'})
        up_button = widgets.Button(description="↑ Up selected file")
        down_button = widgets.Button(description="↓ Down selected file")
        widgets_list.append((up_button, down_button,dropdown_widget))


    if widgets_list:
        # Crear una cuadrícula vertical 4x1
        display(widgets.VBox([widgets.HBox(tupla) for tupla in widgets_list]))
        return tuple(widgets_list)
    else:
        print("No hay listas válidas para mostrar widgets.")
        return ()

    
def update_dropdown(*args):
    """
    Update the options and selected values for multiple Dropdown widgets.
    :param args: Pairs of lists and Dropdown widgets.
    """
    for i in range(0, len(args), 2):
        current_list = args[i]
        current_dropdown = args[i + 1]

        if len(current_list) > 0:
            current_dropdown.options = current_list
            current_dropdown.value = current_list[0]


def replace_options(new_options, dropdown):
    """
    Replace the options of a Dropdown widget.
    :param dropdown: Dropdown widget to be updated.
    :param new_options: List of new options.
    """
    dropdown.options = []

    if len(new_options) > 0:
        dropdown.options = new_options
        dropdown.value = new_options[0]

