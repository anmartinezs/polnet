from vtkmodules.util.numpy_support import vtk_to_numpy
import numpy as np
import vtk
from .utilities import insert_maxis
from polnet import lio, polymer, utils
import os
from scipy.ndimage import rotate, affine_transform
from .tk_utilities import *


keys = {
    'j': "Interactor style: Joystick Actor (j)",
    't': "Interactor style: Joystick Camera (t)",
    'a': "Interactor style: Trackball Actor (a)",
    'c': "Interactor style: Trackball Camera (c)",
    's': "Position saved (s)"
}


def create_poly_mrc(path):
    """
    Create a poly mrc mapper
    :param path: path to mrc file
    :return tuple with vtk objects created
    """
    # tomo_vti = vtk.vtkMRCReader()
    # tomo_vti.SetFileName(path)
    tomo_np = lio.load_mrc(path, mmap=False,no_saxes=False)
    # Normalization
    # print('MRC loaded')
    tomo_np = utils.lin_map(tomo_np, lb=0, ub=1)
    tomo_vti = lio.numpy_to_vti(tomo_np)
    del tomo_np
    # print('VTKimage generated')
    
    iso = vtk.vtkContourFilter()
    iso.SetInputData(tomo_vti)
    iso.SetValue(0, 0.1)  # Ajusta el valor inicial del isovalor (puedes cambiarlo)
    
    iso_mapper = vtk.vtkPolyDataMapper()
    iso_mapper.SetInputConnection(iso.GetOutputPort())
    iso_mapper.ScalarVisibilityOff()
    
    iso_actor = vtk.vtkActor()
    iso_actor.SetMapper(iso_mapper)

    # print('Iso-map constructed')

    return tomo_vti, iso, iso_mapper, iso_actor


def create_window(x,y):
    """ 
    Create a vtk window
    :param x: the width of the window
    :param y: the heigth of the window
    return: tuple with the render, render window and interactor
    """
    ren = vtk.vtkRenderer()
    ren.SetBackground(0,0,0)

    ren_win = vtk.vtkRenderWindow()
    ren_win.AddRenderer(ren)
    ren_win.SetSize(x, y)
    screen_width = ren_win.GetScreenSize()[0]
    screen_height = ren_win.GetScreenSize()[1]
    window_width = ren_win.GetSize()[0]
    window_height = ren_win.GetSize()[1]
    center_x = (screen_width - window_width) // 2
    center_y = (screen_height - window_height) // 2
    ren_win.SetPosition(center_x, center_y)

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(ren_win)
    return ren, ren_win, iren


def visualize_pdb(path, x, y):
    """
    Load an input PDB file and show it
    :param path: the input PDB file
    :param x: the width of the window
    :param y: the heigth of the window
    """
    colors = vtk.vtkNamedColors()
    ren, ren_win, iren = create_window(x,y)

    pdb = vtk.vtkPDBReader()
    pdb.SetFileName(path)
    pdb.SetHBScale(1.0)
    pdb.SetBScale(1.0)
    pdb.Update()

    resolution = pow(300000.0 / pdb.GetNumberOfAtoms(), 0.5)
    resolution = min(max(resolution, 4), 20)

    sphere = vtk.vtkSphereSource()
    sphere.SetCenter(0, 0, 0)
    sphere.SetRadius(1)
    sphere.SetThetaResolution(int(resolution))
    sphere.SetPhiResolution(int(resolution))

    glyph = vtk.vtkGlyph3D()
    glyph.SetInputConnection(pdb.GetOutputPort())
    glyph.SetOrient(1)
    glyph.SetColorMode(1)
    glyph.SetScaleMode(2)
    glyph.SetScaleFactor(0.25)
    glyph.SetSourceConnection(sphere.GetOutputPort())

    atomMapper = vtk.vtkPolyDataMapper()
    atomMapper.SetInputConnection(glyph.GetOutputPort())
    atomMapper.UseLookupTableScalarRangeOff()
    atomMapper.ScalarVisibilityOn()
    atomMapper.SetScalarModeToDefault()

    atom = vtk.vtkLODActor()
    atom.SetMapper(atomMapper)
    atom.GetProperty().SetRepresentationToSurface()
    atom.GetProperty().SetInterpolationToGouraud()
    atom.GetProperty().SetAmbient(0.1)
    atom.GetProperty().SetDiffuse(0.7)
    atom.GetProperty().SetSpecular(0.5)
    atom.GetProperty().SetSpecularPower(80)
    atom.GetProperty().SetSpecularColor(colors.GetColor3d("White"))

    ren.AddActor(atom)

    tube = vtk.vtkTubeFilter()
    tube.SetInputConnection(pdb.GetOutputPort())
    tube.SetNumberOfSides(int(resolution))
    tube.CappingOff()
    tube.SetRadius(0.2)
    tube.SetVaryRadius(0)
    tube.SetRadiusFactor(10)

    bondMapper = vtk.vtkPolyDataMapper()
    bondMapper.SetInputConnection(tube.GetOutputPort())
    bondMapper.UseLookupTableScalarRangeOff()
    bondMapper.ScalarVisibilityOff()
    bondMapper.SetScalarModeToDefault()

    bond = vtk.vtkLODActor()
    bond.SetMapper(bondMapper)
    bond.GetProperty().SetRepresentationToSurface()
    bond.GetProperty().SetInterpolationToGouraud()
    bond.GetProperty().SetAmbient(0.1)
    bond.GetProperty().SetDiffuse(0.7)
    bond.GetProperty().SetSpecular(0.5)
    bond.GetProperty().SetSpecularPower(80)
    bond.GetProperty().SetSpecularColor(colors.GetColor3d("White"))

    ren.AddActor(bond)
    base_name = os.path.splitext(os.path.basename(path))[0]
    ren_win.SetWindowName(base_name.upper())
    iren.Initialize()
    ren_win.Render()
    iren.Start()
    

def create_slider(iren):
    """
    Load an input PDB file and show it
    :param iren: window interactor
    :return: tuple with the sliders
    """
    # Create iso surface slider
    slider_rep = vtk.vtkSliderRepresentation2D()
    slider_rep.SetMinimumValue(0.0)  # Maximun value
    slider_rep.SetMaximumValue(1.0)  # Minimun value
    slider_rep.SetValue(0.1)  # Default value
    slider_rep.GetPoint1Coordinate().SetCoordinateSystemToNormalizedDisplay()
    slider_rep.GetPoint1Coordinate().SetValue(0.1, 0.1)
    slider_rep.GetPoint2Coordinate().SetCoordinateSystemToNormalizedDisplay()
    slider_rep.GetPoint2Coordinate().SetValue(0.9, 0.1)
    slider_rep.SetTitleText("Iso Value")

    slider_widget = vtk.vtkSliderWidget()
    slider_widget.SetInteractor(iren)
    slider_widget.SetRepresentation(slider_rep)
    slider_widget.SetAnimationModeToAnimate()
    return slider_rep, slider_widget


def select_isosurface(path, x, y):
    """
    Visualize protein and select isosurface
    :param path: the input PDB file
    :param x: the width of the window
    :param y: the heigth of the window
    :return: isosurface selected value
    """
    ren, ren_win, iren = create_window(x,y)

    tomo_vti, iso, iso_mapper, iso_actor = create_poly_mrc(path)

    iso_actor.GetProperty().SetColor(1.0, 1.0, 0.7)  # Color amarillo
    outline = vtk.vtkOutlineFilter()
    outline.SetInputData(tomo_vti)

    outline_mapper = vtk.vtkPolyDataMapper()
    outline_mapper.SetInputConnection(outline.GetOutputPort())

    outline_actor = vtk.vtkActor()
    outline_actor.SetMapper(outline_mapper)

    ren.AddActor(outline_actor)
    ren.AddActor(iso_actor)

    iren.Initialize()

    slider_rep, slider_widget = create_slider(iren)
    
    def update_iso_surface(obj, event):
        value = slider_widget.GetRepresentation().GetValue()
        iso.SetValue(0, value)
        ren_win.Render()

    slider_widget.AddObserver("InteractionEvent", update_iso_surface)

    slider_widget.On()

    base_name = os.path.splitext(os.path.basename(path))[0]
    ren_win.SetWindowName(base_name.upper())
    iren.Start()

    return  slider_widget.GetRepresentation().GetValue()


def protein_to_axis(membrane_path, axis_path, x, y, v_size, outpath):
    """
    Align membrane protein in axis
    :param membrane_path: path to membrane protein
    :param axis_path: path to ref axis
    :param x: the width of the window
    :param y: the height of the window
    :param v_size: voxel size
    :param outpath: path to save 
    """
    
    ren, ren_win, iren = create_window(x, y)

    tomo_vti, iso, iso_mapper, iso_actor = create_poly_mrc(membrane_path)
    iso_actor.GetProperty().SetColor(1.0, 1.0, 0.7)  # Color amarillo
   
    tomo_vti2, iso2, iso_mapper2, iso_actor2 = create_poly_mrc(axis_path)
    iso_actor2.GetProperty().SetColor(0.6, 0.6, 0.6)  # Color gris claro
    iso_actor2.PickableOff()

    
    actor_matrix_initial = iso_actor.GetMatrix()
    transform_initial = vtk.vtkTransform()
    transform_initial.SetMatrix(actor_matrix_initial)
    
    transform_filter_initial = vtk.vtkTransformFilter()
    #transform_filter_initial.SetInputConnection(tomo_vti.GetOutputPort())
    transform_filter_initial.SetInputData(tomo_vti)
    transform_filter_initial.SetTransform(transform_initial)
    transform_filter_initial.Update()
    print("Initial Rotation: ", transform_initial.GetOrientation())
    
    # Calculo centro de masas
    center_mass_o = vtk.vtkCenterOfMass()
    center_mass_o.SetInputConnection(iso.GetOutputPort())
    center_mass_o.Update()
    c_o = center_mass_o.GetCenter()
    print("Center of mass vtkCenterOfMass:", c_o)
    
    ren.AddActor(iso_actor)
    ren.AddActor(iso_actor2)

    # Create axis bounding box
    outline_filter_axis = vtk.vtkOutlineFilter()
    outline_filter_axis.SetInputConnection(iso_actor2.GetMapper().GetInputConnection(0, 0))  # Use the same mapper as membrane
    
    outline_mapper_axis = vtk.vtkPolyDataMapper()
    outline_mapper_axis.SetInputConnection(outline_filter_axis.GetOutputPort())
    
    # Crea un actor para el bounding box del eje
    outline_actor_axis = vtk.vtkActor()
    outline_actor_axis.SetMapper(outline_mapper_axis)
    outline_actor_axis.GetProperty().SetColor(0.5, 0.5, 0.5)  
    outline_actor_axis.PickableOff()  # disable mouse selected
    ren.AddActor(outline_actor_axis)

    
    def update(obj, event):
        iren.Render()
        key = obj.GetKeySym().lower()
        text = keys.get(key, "No function for this key")
        text_actor.SetInput(text)
        if text.lower() == "position saved (s)":
            
            actor_matrix = iso_actor.GetMatrix()
            transform = vtk.vtkTransform()
            transform.SetMatrix(actor_matrix)

            # Calculate center of mass
            transform_filter = vtk.vtkTransformFilter()
            transform_filter.SetInputData(tomo_vti)
            #transform_filter.SetInputConnection(tomo_vti.GetOutputPort())
            transform_filter.SetTransform(transform)
            transform_filter.Update()
            
            center_mass_o.SetInputConnection(transform_filter.GetOutputPort())
            center_mass_o.Update()
            
            c_m = center_mass_o.GetCenter()
            angles = transform.GetOrientation()

            print("Center of mass o ", c_o)
            print("Center of mass m " , c_m)
            print("Rotation ", angles)
        
            rotation_transform = vtk.vtkTransform()
            rotation_transform.RotateZ(angles[2])
            rotation_transform.RotateX(angles[0])
            rotation_transform.RotateY(angles[1])
            rotation_transform.Inverse()

            resliced = vtk.vtkImageReslice()
            resliced.SetInputData(tomo_vti)
            resliced.SetAutoCropOutput(True)
            resliced.SetResliceTransform(rotation_transform)
            resliced.SetInterpolationModeToLinear()
            resliced.Update()

            image_data = resliced.GetOutput()
            image_array = vtk_to_numpy(image_data.GetPointData().GetScalars()).reshape(image_data.GetDimensions(), order='F')

            #eje_array = np.zeros(shape=tomo_vti2.GetOutput().GetDimensions(), dtype=np.float32)
            eje_array = np.zeros(shape=tomo_vti2.GetDimensions(), dtype=np.float32)
            p = insert_maxis(image_array, eje_array, np.asarray(c_m), outpath, os.path.splitext(os.path.split(membrane_path)[1])[0])
            window_align_to_axis(p)
 
    
    iren.AddObserver("KeyPressEvent", update)
    
    text_actor = vtk.vtkTextActor()
    text_actor.SetInput('Interactor style: joystick actor')
    text_actor.SetPosition(0.0, 0.1)
    
    text_representation = vtk.vtkTextRepresentation()
    text_representation.GetPositionCoordinate().SetValue(0.01, 0.92)
    text_representation.GetPosition2Coordinate().SetValue(0.35, 0.08)
    
    text_widget = vtk.vtkTextWidget()
    text_widget.SetRepresentation(text_representation)
    text_widget.SetInteractor(iren)
    text_widget.SetTextActor(text_actor)
    text_widget.SelectableOff()
    
    text_widget.On()
    
    ren_win.SetWindowName("Align membrane protein")
    ren_win.Render()
    iren.Start()


def visualize_helix(data, x, y, path):
    """
    Create a vtk window to visualize helix
    :param data: VTK PolyData object
    :param x: the width of the window
    :param y: the height of the window
    :param path: path to take the name for the window 
    """
    ren, ren_win, iren = create_window(x, y)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(data)
    
    scalar_range = data.GetScalarRange()
    mapper.SetScalarRange(scalar_range)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    ren.AddActor(actor)

    outline_filter = vtk.vtkOutlineFilter()
    outline_filter.SetInputData(data)

    outline_mapper = vtk.vtkPolyDataMapper()
    outline_mapper.SetInputConnection(outline_filter.GetOutputPort())

    outline_actor = vtk.vtkActor()
    outline_actor.SetMapper(outline_mapper)
    outline_actor.GetProperty().SetColor(0.5, 0.5, 0.5)  #   # Color del bounding box 1.0, 1.0, 0.7

    ren.AddActor(outline_actor)

    base_name = os.path.splitext(os.path.basename(path))[0]
    ren_win.SetWindowName(base_name.upper())
    iren.Initialize()
    ren_win.Render()
    iren.Start()




