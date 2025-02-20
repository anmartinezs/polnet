Scripts exploiting the functionality of polnet module.

* **csv**: manipulating output CSV files with the particles (or motifs) lists
* **templates**: scripts for preprocessing input templates for macromolecules.
* **data_gen**: scripts for synthetic data generation
  * The script **all_features.py** can generate tomograms with all structural features available. It is set up to generate
    a dataset with all different types of structures including membranes, microtubules, actin networks, 18 different
    cytosolic and 8 membrane-bound macromolecules. Input configuration files for geometry and transformation generators,
    as well as the macromolecular models are available under **data** folder.
* **deprecated**: contains
some scripts for evaluations carried out during the software development, they are not prepared for external users
because some hardcoded paths need to be modified.

# Run polnet with docker

## Docker build
bash create_docker.sh

## Launch data generation
bash run_docker.sh --config path/to/config.yaml --out_dir /path/to/out_dir

The data generation is divided into **two steps**. The first one creates the digital sample and the second one simulates an acquisition. You can simulate different acquisitions for the same digital sample. All parameters are defined in a yaml file. 

NB :
* examples of config.yaml files are here: ./scripts/config_acquisition.yaml and ./scripts/config_sample.yaml 
* a dir will be created in out_dir (corresponding to the name written in config.yaml) and all data will
be saved there. If you simulate different acquisitions for the same digital sample please change the name of your acquisition, otherwise it will overwrite the previous one.
