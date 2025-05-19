# CIAA: Integrated Proteomics and Structural Modeling for Understanding Cysteine Reactivity with Iodoacetamide Alkyne

Cysteine residues play key roles in protein structure and function and can serve as targets for chemical probes and even drugs. Chemoproteomic studies have revealed that heightened cysteine reactivity towards electrophilic probes, such as iodoacetamide alkyne (IAA), is indicative of likely residue functionality. However, while the cysteine coverage of chemoproteomic studies has increased substantially, these methods still only provide a partial assessment of proteome-wide cysteine reactivity, with cysteines from low abundance proteins and tough-to-detect peptides still largely refractory to chemoproteomic analysis. Here we integrate cysteine chemoproteomic reactivity datasets with structure-guided computational analysis to delineate key structural features of proteins that favor elevated cysteine reactivity towards IAA.

## About
[https://chemrxiv.org/engage/chemrxiv/article-details/678722edfa469535b90fe695](https://chemrxiv.org/engage/api-gateway/chemrxiv/assets/orp/resource/item/67872b91fa469535b910ab1e/original/supporting-information.pdf)

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
* dssp (for SS & SASA of cysteines; https://swift.cmbi.umcn.nl/gv/dssp/index.html)
* propka

## Installation

```bash
conda create -n cysteine_reactivity python=3
conda activate cysteine_reactivity
conda install -c conda-forge -c salilab numpy scipy pandas matplotlib scikit-learn matplotlib \
  pymol-open-source mdanalysis modeller propka notebook
```
