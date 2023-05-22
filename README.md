# PolNet
Python package for generating synthetic datasets of the cellular context for Cryo-Electron Tomography.

## Installation

* Get the python code package from github: https://github.com/anmartinezs/polnet.git
* Install all requirements (using pip or conda): see requirements.txt file (JAX is optional).
* Install **IMOD**, PolNet may call to some of its standalone commands: https://bio3d.colorado.edu/imod/doc/guide.html
* **EMAN2** is required to generate macromolecular templates.

### Installation notes

* **mrcfile**: python package installation with conda may need adding the channel conda-forge, see instructions in: https://pypi.org/project/mrcfile/
* **jax**: python package installation with GPU support is required for faster data generation: https://github.com/google/jax#installation

## Content

* **polnet**: python package with the functionality for generating the synthetic data.
* **scripts**: python scripts for generating different types of synthetic datasets. IMPORTANT NOTE: for the moment some paths are hardcoded, so you must edit them to successfully run them on your computer.
* **tests**: unit tests for functionalities in **polnet**. IMPORTANT NOTE: for the moment some paths are hardcoded, so you must edit them to successfully run them on your computer.

## Generating a generic synthetic dataset

If you want to generate a generic synthetic dataset for cryo-ET without any specific requirement, you must run the script
**scripts/data_gen/all_features.py**, it is set up to generate a dataset with all different types of structures including
membranes, microtubules, actin networks, 18 different cytosolic and 8 membrane-bound macromolecules.

Input configuration files for geometry and transformation generators, as well as the macromolecular models are available
under **data** folder.

A PDF with a text description of the setting used in the scripts **scripts/data_gen/all_features.py** is in **docs/default_settings.pdf**. The mathematical 
symbols and acronyms correspond with the definitions in the **main publication** document.

## Main publication (Citation)

[1] Martinez-Sanchez A., Jasnin M. and Phelipeau H. "Simulating the cellular context in synthetic datasets for cryo-electron tomography" *In preparation*



