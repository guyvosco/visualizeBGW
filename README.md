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

* Band structures from mean-field and quasiparticle corrections
* Wavefunction character and orbital projections
* Components of the dielectric matrix
* Optical absorption spectra (RPA and BSE)
* Exciton contributions in real and reciprocal space

---

## Project Structure (Draft)

```
visualizeBGW/
│
├─ src/                 # Core visualization library
│   ├─ io/                 # Input/output parsers
│   ├─ analysis/           # Data manipulation utilities
│   ├─ plotting/           # Matplotlib-based plotting functions
│   └─ gui/                # GUI implementation
├─ examples/
├─ tests/
├─ README.md
└─ requirements.txt
```

### Features (Planned)

#### Data Parsing
* Parse and process `.dat` outputs
  - [ ] chi_converge.dat
  - [ ] eqp.dat
  - [ ] eigenvalues.dat
* Read BerkeleyGW `HDF5` outputs
  - [ ] WFN.h5
  - [ ] epsmat.h5
  - [ ] eigenvectors.h5

#### Core Calculations
  - [ ] High-symmetry paths
  - [ ] Band-structure interpolation
  - [ ] Excitonic component weights

#### Visualization Tools
  - [ ] DFT and GW-corrected quasiparticle band structure
  - [ ] Convergence of $\chi$ with respect to empty orbitals
  - [ ] *head* of the dielectric function $\epsilon(\omega)$
  - [ ] RPA and BSE absorption spectra
  - [ ] Exciton band contributions and $k$-resolved weights

#### GUI (Tkinter / PyQt / Streamlit – TBD)
* File browser for loading BerkeleyGW output folders
* Tabs for each visualization module
* Interactive parameter selection (energy windows, $k$-paths, etc.)
* Export plots to PNG/PDF

### Stretch Goals (Optional)
- [ ] 3D crystal structure visualization
- [ ] Real-space wavefunction projections
- [ ] Real-space exciton visualization

---

## Dependencies

* **Python 3**
* **h5py** for HDF5 reading
* **NumPy / SciPy** for numerical processing
* **matplotlib** for plotting
* **PyQt6 / Tkinter / Streamlit** for the GUI (final choice TBD)
* **pytest** for testing
