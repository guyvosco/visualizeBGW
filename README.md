# visualizeBGW

*A Python toolkit and GUI for exploring BerkeleyGW outputs*

## Overview

**visualizeBGW** provides an extensible Python visualization library for [BerkeleyGW](https://berkeleygw.org/), a widely used many-body perturbation theory package for computing quasiparticle energies and excitonic properties.
This library streamlines the workflow of researchers who want to quickly explore, visualize, and understand their BerkeleyGW calculations without writing ad-hoc scripts.

---

## Motivation

BerkeleyGW offers robust tools for large-scale GW and BSE calculations, but exploring its outputs typically requires custom scripts and manual data extraction. BerkeleyGW generates complex datasets stored across multiple `HDF5` and `.dat` files, making direct inspection time-consuming and error-prone.
visualizeBGW provides a unified interface and visualization toolkit that accelerates data analysis, improves reproducibility, and simplifies scientific interpretation.
With this library, users will be able to interactively explore:

* Wavefunction character and orbital projections
* Band structures from mean-field and quasiparticle corrections
* Components of the dielectric matrix
* Optical absorption spectra (RPA and BSE)
* Exciton contributions

---

## Project Structure

```
visualizeBGW/
          │
          ├─ src/visualizeBGW/             # Core visualization library
          │                ├─ io/             # Input/output parsers
          │                ├─ analysis/       # Data manipulation utilities
          │                ├─ plotting/       # Matplotlib-based plotting functions
          │                └─ gui/            # GUI implementation
          │                    └─ pages/         # Pages classes
          ├─ examples/                     # QE and BerkeleyGW calculations inputs
          │        ├─ Si/                     # Bulk Silicon
          │        └─ MoS2/                   # Transition metal dichalcogenide monolayer
          ├─ tests/
          ├─ README.md
          └─ requirements.txt
```

### Features:

#### Data Parsing
* Parse and process BerkeleyGW outputs
  - chi_converge.dat
  - eqp.dat
  - eigenvalues.dat
* Read and analyse BerkelyGW HDF5 formats
  - WFN.h5
  - epsmat.h5
  - eigenvectors.h5

#### Export files
* POSCAR - stractual geometry file format
* XSF - 3D electron densities data

#### Core Calculations
* Electron densities in real-space
* High-symmetry paths
* Broaden absorption spectrum
* Excitonic component weights

#### Visualization Tools
* Crystal structure visualization
* Mean-field wavefunction projection in real-space
* DFT and GW-corrected quasiparticle band structure
* Convergence of $\chi$ with respect to empty orbitals
* Convergence of the Coulomb-hole sum value vs. the number of bands included in the sum
* *head* of the dielectric function $\epsilon(\omega)$
* RPA and BSE absorption spectra
* Exciton band contributions and $k$-resolved weights

#### GUI
* Tabs for each visualization module
* File browser for loading BerkeleyGW output files
* Interactive parameter selection
* Export structure visualization data

---

## Dependencies

* **Python 3**
* **h5py** for HDF5 reading
* **NumPy / SciPy** for numerical processing
* **ase** for atomic properties data
* **matplotlib** for plotting
* **PyVista** for 3D visualiztion
* **PySide6** for the GUI
* **pytest** for testing
