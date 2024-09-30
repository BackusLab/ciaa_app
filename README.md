# ciaa_app


An insilico structure-based approach to identify IAA reactive cysteines.

## Requirements
* Python 3
* Numpy
* Scipy
* Pandas
* Matplotlib
* Scikit-learn (ML models)
* MDAnalysis (used for calculating HB)
* Prody (for building biological units)
* Modeller (for building incomplete sidechains)
* dssp (for SS & SASA of cysteines; https://swift.cmbi.umcn.nl/gv/dssp/index.html)
* stride (for calculating secondary structure; http://webclu.bio.wzw.tum.de/stride/)
* propka

## Installation

```bash
conda create -n cysteine_reactivity python=3
conda activate cysteine_reactivity
conda install -c conda-forge -c salilab numpy scipy pandas matplotlib scikit-learn matplotlib \
  pymol-open-source mdanalysis prody modeller propka notebook
```
