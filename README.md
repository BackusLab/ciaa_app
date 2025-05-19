# CIAA: Integrated Proteomics and Structural Modeling for Understanding Cysteine Reactivity with Iodoacetamide Alkyne

Cysteine residues play key roles in protein structure and function and can serve as targets for chemical probes and even drugs. Chemoproteomic studies have revealed that heightened cysteine reactivity towards electrophilic probes, such as iodoacetamide alkyne (IAA), is indicative of likely residue functionality. However, while the cysteine coverage of chemoproteomic studies has increased substantially, these methods still only provide a partial assessment of proteome-wide cysteine reactivity, with cysteines from low abundance proteins and tough-to-detect peptides still largely refractory to chemoproteomic analysis. Here we integrate cysteine chemoproteomic reactivity datasets with structure-guided computational analysis to delineate key structural features of proteins that favor elevated cysteine reactivity towards IAA.

## About
[https://chemrxiv.org/engage/chemrxiv/article-details/678722edfa469535b90fe695](https://chemrxiv.org/engage/api-gateway/chemrxiv/assets/orp/resource/item/67872b91fa469535b910ab1e/original/supporting-information.pdf)

## File Descriptions

### Input Files
* `data`: directory that contains all the csv files
* `isotop_pdb.csv`: input raw dataset
* `isotop_training_set.csv`: input training dataset
* `isotop_test_set.csv`: input test dataset
* `isotop_alphafold_pdb_set.csv`: input alphafold validation dataset
* `isotop_alphafold_denovo_set.csv`: input alphafold validation dataset

### Analysis Files
* `isotop_residue_function.csv`: input for isotop_residue_function_barplot.ipynb
* `list_found_peptides.csv`: input for isotop_pdb_filtering.ipynb
* `final_selection.csv`: input for isotop_pdb_filtering.ipynb
* `isotop_training_nonredundant_complete_final_identifiers.csv`: input for isotop_pdb_filtering.ipynb
* `isotop_pdb_meric_state.csv`: input for isotop_pdb_filtering.ipynb
* `isotop_descriptors.csv`: input dataset for isotop_descriptor_pearson_correlations.ipynb

### Scripts
* `make_dataset.py`: Prepare PDB structure list using input files
* `download_pdbs.py`: Download PDB files from input pdb_files.txt
* `reduce_pdbs.sh`: Clean PDBs of heteroatoms, ligands, and add hydrogens
* `get_neighbor_graphs.py`: Find neighbors of cysteines in PDB structure list
* `get_simple_descriptors.py`: Calculate raw structural descriptors of cysteines
* `get_simple_descriptors.py`: One hot encode raw structural descriptors of cysteines

### Jupyter Notebooks
* `isotop_residue_function_barplot.ipynb`: Analysis for classifying cysteines based on function
* `isotop_pdb_filtering.ipynb`: Analysis for filtering PDB structures
* `isotop_descriptor_pearson_correlations.ipynb`: Analysis for computing pearson correlation coefficients 
* `ciaa_hyperreactive_cysteines_model.ipynb`: Analysis for developing the CIAA model based on input files

### Output Files
* `ciaa_results`: directory that contains all the csv files and png images from the model 
* `ciaa_rf_classifier_model.pkl`: output from ciaa_hyperreactive_cysteines_model.ipynb

## Requirements
* Python 3
* Numpy
* Scipy
* Pandas
* Matplotlib
* Scikit-learn (ML models)
* MDAnalysis (used for calculating HB)
* Prody (for building biological units; http://www.bahargroup.org/prody/manual/getprody.html)
* Modeller (for building incomplete sidechains)
* Reduce (for adding hydrogens; http://kinemage.biochem.duke.edu/software/reduce/)
* DSSP (for SS & SASA of cysteines; https://swift.cmbi.umcn.nl/gv/dssp/index.html)
* propka

## Installation

```bash
conda create -n cysteine_reactivity python=3
conda activate cysteine_reactivity
conda install -c conda-forge -c salilab numpy scipy pandas matplotlib scikit-learn matplotlib \
  pymol-open-source mdanalysis modeller propka notebook
```

## Usage

### Create Descriptors of Cysteine Microenvironments from PDB Structures

### Step 0: Create List of PDBs from IsoTOP 2025 (Optional)
```
python3 scripts/make_dataset.py
```

### Step 1: Create List of PDBs
Create input file with a list of PDB accessions called pdb_files.txt and store in the data folder (example: data/pdb_files.txt)

### Step 2: Download PDBs
```
python3 scripts/download_pdbs.py -i 'pdb_files.txt'
```

### Step 3: Clean PDBs of Heteroatoms and Ligands
```
bash ../scripts/reduce_pdbs.sh
```

### Step 4: Get Neighbors of Cysteines
```
python3 scripts/get_neighbor_graphs.py
```

### Step 5: Create Raw Descriptors
Update /Users/{user_name}/anaconda3/bin/mkdssp to access local installation of dssp

```
python3 scripts/get_simple_descriptors.py
```

### Step 6: One Hot Encode Raw Descriptors
```
python3 scripts/get_simple_descriptors.py
```

### Step 7: Map PDB Residues to [SIFTS](https://github.com/lmboat/protein_structure_annotations) (optional)

## Citation
Boatner, L., Eberhardt, J., Shikwana, F., Holcomb, M., Lee, P., Houk, K., Forli, S. & Backus, K. (2025). CIAA: Integrated Proteomics and Structural Modeling for Understanding Cysteine Reactivity with Iodoacetamide Alkyne.
