# PolNet
Python package for generating synthetic datasets of the cellular context for Cryo-Electron Tomography.

## Installation

* Get the python code package from github: https://github.com/anmartinezs/polnet.git
* Install all requirements (using pip or conda): see requirements.txt file
* Install **IMOD**, PolNet may call to some of its standalone commands: https://bio3d.colorado.edu/imod/doc/guide.html

### Installation notes

* **mrcfile**: python package installation with conda may need adding the channel conda-forge, see instructions in: https://pypi.org/project/mrcfile/
* **jax**: python package installation with GPU support is required for faster data generation: https://github.com/google/jax#installation

