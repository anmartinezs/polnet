# PolNet
Python package for generating synthetic datasets of the cellular context for Cryo-Electron Tomography.

## Installation

Here is how to get it installed:

    git clone anmartinezs/polnet
    cd polnet
    conda create --name polnet pip
    conda activate polnet
    pip install -e .

* You can check all requirements in the requirements.txt file (JAX is optional).
* Install **IMOD**, PolNet calls to some of its standalone commands: https://bio3d.colorado.edu/imod/doc/guide.html

### Installation with EMAN2

If you want to generate you own macromolecular templates you need to **install EMAN2 first** because its script **e2pdb2mrc.py** is needed.
In such case the installation steps are:

1. Install [EMAN2](https://blake.bcm.edu/emanwiki/EMAN2) (we strongly recommend v2.91 since it was tested).
2. Activate EMAN2 conda environment:
    ```
    conda activate eman-install-path/or/env-name
    ```
3. Install **PolNet** using EMAN's conda environment:
    ```
    git clone anmartinezs/polnet
    cd polnet
    pip install -r requirements.txt
    ```

### Installation notes

* **mrcfile**: python package manual installation with conda may need adding the channel conda-forge, see instructions in: https://pypi.org/project/mrcfile/
* **jax**: python package installation with GPU support is required for faster data generation: https://github.com/google/jax#installation
* **pynrrd**: this python package is only needed to run the script for preparing nn-UNet training datasets

## Content

* **polnet**: python package with the Python implemented functionality for generating the synthetic data.
* **scripts**: python scripts for generating different types of synthetic datasets. Folders:
  + **data_gen**: scripts for data generation.
    * **deprecated**: contains 
    some scripts for evaluations carried out during the software development, they are not prepared for external users
    because some hardcoded paths need to be modified.
  + **templates**: scripts for building the structural units for macromolecules (requires installation with **EMAN2**)
  + **csv**: scripts for postprocessing the CSV generated files.
  + **data_prep**: script to convert the generated dataset in [nn-UNet](https://github.com/MIC-DKFZ/nnUNet) format.
* **tests**: unit tests for functionalities in **polnet**. The script **tests/test_transformations.py** requires to generate at 
least 1 output tomo with the script **scripts/all_features.py** and modified the hardcoded input paths, that is because
the size of the input data avoid to upload them to the repository.
* **docs**: contains a PDF with the suplementary material for [1] with the next tables:
  + Glossary of acronyms by order of appearance in the main text.
  + Glossary mathematical symbols defined in the main text organized by scope
  + Table Variables used by the input files to model the generators.
  + Table with the structures used to simulate the cellular context. 

## Generating a generic synthetic dataset

If you want to generate a generic synthetic dataset for cryo-ET without any specific requirement, you must run the script
**scripts/data_gen/all_features.py**:

    cd scripts
    python data_gen/all_features.py

It is set up to generate a dataset with all different types of structures including membranes, microtubules, actin networks, 18 different cytosolic and 8 membrane-bound macromolecules.
Input configuration files for geometry and transformation generators, as well as the macromolecular models are available
under **data** folder.

A PDF with a text description of the default settings used in the scripts **scripts/data_gen/all_features.py** is in **docs/default_settings.pdf**. The mathematical 
symbols and acronyms correspond with the definitions in the **main publication** document.

## Main publication (Citation)

[1] Martinez-Sanchez A.*, Jasnin M., Phelippeau H. and Lamm L. "Simulating the cellular context in synthetic datasets for cryo-electron tomography" *bioRxiv*



