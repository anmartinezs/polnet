# PolNet
Python package for generating synthetic datasets of the cellular context for Cryo-Electron Tomography.

## Installation

* Get the python code package from github: https://github.com/anmartinezs/polnet.git
* Install all requirements (using pip or conda): see requirements.txt file (JAX is optional).
* Install **IMOD**, PolNet may call to some of its standalone commands: https://bio3d.colorado.edu/imod/doc/guide.html
* **EMAN2** is required to generate macromolecular templates as the script **e2pdb2mrc.py** is used.

### Installation notes

* **mrcfile**: python package installation with conda may need adding the channel conda-forge, see instructions in: https://pypi.org/project/mrcfile/
* **jax**: python package installation with GPU support is required for faster data generation: https://github.com/google/jax#installation

To install it:

    git clone anmartinezs/polnet
    cd polnet
    conda create --name polnet pip
    pip install -e .

## Content

* **polnet**: python package with the functionality for generating the synthetic data.
* **scripts**: python scripts for generating different types of synthetic datasets. Subfolder **deprecated** contains 
some scripts for evaluations carried out during the software development, they are not prepared for external users
because some hardcoded paths need to be modified.
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

A PDF with a text description of the setting used in the scripts **scripts/data_gen/all_features.py** is in **docs/default_settings.pdf**. The mathematical 
symbols and acronyms correspond with the definitions in the **main publication** document.

## Main publication (Citation)

[1] Martinez-Sanchez A.*, Jasnin M., Phelippeau H. and Lamm L. "Simulating the cellular context in synthetic datasets for cryo-electron tomography" *bioRxiv*



