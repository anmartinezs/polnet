# PolNet


[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://github.com/anmartinezs/polnet/blob/main/LICENSE)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/anmartinezs)


Python package for generating synthetic datasets of the cellular context for Cryo-Electron Tomography.

## Installation

### Requirements

- **IMOD** must be installed on the system since PolNet calls to some of its standalone commands: https://bio3d.colorado.edu/imod/doc/guide.html
- Miniconda or Anaconda with Python 3.
- Git.
- IMOD can be used for MRC files visualization. Paraview can be used for VTK (.vtp) files visualization. Pandas is recommended for managing the CSV files.

### Installation procedure
Here is how to get it installed:

1. Download PolNet source code:
    ```console
    git clone https://github.com/anmartinezs/polnet.git
    cd polnet
    ```
2. Create a conda virtual environment
    ```console
    conda create --name polnet pip
    conda activate polnet
    ```
   
3. Install PolNet package with its requirements:
    ```console
    pip install -e .
    ```
**For developers** who do not want to install PolNet in the virtual environment as a package, you can only install 
the requirements by:

    pip install -r requirements.txt

You can check all requirements in the **requirements.txt** file (JAX is optional).

The installation has been tested in Ubuntu 22.04 and Windows 10 and 11.

## Usage

First open **Jupyter** by running:

```console
jupyter notebook
```
If you are not familiar with jupyter notebooks, first get started with https://docs.jupyter.org/en/latest/running.html 

To generate a synthetic dataset run on **Jupyter** next notebook: **gui/gen_dataset.ipynb**

To create you own structural models next Jupyter notebooks are available:
    
1. Membranes:  **gui/create_membrane_models.ipynb**
2. Filaments:  **gui/create_filament_models.ipynb**
3. Macromolecules: 
   - Atomic model (PDB) to electron density map (MRC): **gui/atomic_to_density.ipynb**
   - Only for membrane bound macromolecules: **gui/align_membrane_proteins.ipynb**
   - Models: **gui/create_macromolecule_models.ipynb**

Exemplary videos for using the GUI are avilable at [Zenodo](https://zenodo.org/records/10735525).

**Important note**: all Jupyter notebooks are thoroughly self-documented in order to guide the user in the process. In addition, they contain graphic objects and default setting to facilitate the process.

## For developers

### Package description
* **polnet**: python package with the Python implemented functionality for generating the synthetic data.
* **gui**: set of Jupyter notebooks with Graphic User Interface (GUI).
  * **core**: functionality required by the notebooks.
* **scripts**: python scripts for generating different types of synthetic datasets. Folders:
  + **data_gen**: scripts for data generation.
    * **deprecated**: contains 
    some scripts for evaluations carried out during the software development, they are not prepared for external users
    because some hardcoded paths need to be modified.
      * **templates**: scripts for building the structural units for macromolecules (requires the installation **EMAN2**). Their usage is strongly deprecated, now GUI notebooks include all functionality.
  + **csv**: scripts for postprocessing the CSV generated files.
  + **data_prep**: script to convert the generated dataset in [nn-UNet](https://github.com/MIC-DKFZ/nnUNet) format.
* **tests**: unit tests for functionalities in **polnet**. The script **tests/test_transformations.py** requires to generate at 
least 1 output tomo with the script **scripts/all_features.py** and modified the hardcoded input paths, that is because
the size of the input data avoid to upload them to the repository.
* **data**: contains input data, mainly macromolecules densities and configuration input files, that con be used to simulate tomograms. These are the default input, an user can add/remove/modify these input data using the notebooks in **GUI**.
  * **in_10A**: input models for macromolecules at 10A voxel size.
  * **in_helix**: input models for helical structures.
  * **in_mbsx**: input models for membrane structures.
  * **tempaltes**: atomic models and density maps used by macromolecular models.
* **docs**:
  * API documentation. 
  * A PDF with the suplementary material for [1] with the next tables:
    + Glossary of acronyms by order of appearance in the main text.
    + Glossary mathematical symbols defined in the main text organized by scope
    + Table Variables used by the input files to model the generators.
    + Table with the structures used to simulate the cellular context. 

### Code documentation

The API documentation for polnet Python package is available in [docs/apidoc/index.html](http://htmlpreview.github.io/?https://raw.githubusercontent.com/anmartinezs/polnet/main/docs/apidoc/index.html) 


## Main publication (Citation)

[1] Martinez-Sanchez A.*, and Lamm L., Jasnin M. and Phelippeau H. (2024) "Simulating the cellular context in synthetic datasets for cryo-electron tomography" *IEEE Transactions on Medical Imaging* [10.1109/TMI.2024.3398401](https://doi.org/10.1109/TMI.2024.3398401)



