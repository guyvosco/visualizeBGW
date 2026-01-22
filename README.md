# visualizeBGW

*A Python toolkit and GUI for exploring BerkeleyGW outputs*

## Overview

**visualizeBGW** provides an extensible Python visualization library for [BerkeleyGW](https://berkeleygw.org/), a widely used many-body perturbation theory package for computing quasiparticle energies and excitonic properties.
This library streamlines the workflow of researchers who want to quickly explore, visualize, and understand their BerkeleyGW calculations without writing ad-hoc scripts.

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
* Read and analyse BerkeleyGW HDF5 datasets
  - WFN.h5
  - epsmat.h5
  - eigenvectors.h5

#### Export files
* POSCAR - structural geometry file format
* XSF - 3D electron densities data

#### Core Calculations
* Electron densities in real-space
* High-symmetry paths
* broaden absorption spectra
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

## Dependencies and Installation

* **Python 3.9+**
* **h5py** for HDF5 reading
* **NumPy / SciPy** for numerical processing
* **ase** for atomic properties data
* **matplotlib** for plotting
* **PyVista** for 3D visualization
* **PySide6** for the GUI
* **pytest** for testing

### Installation

#### Option 1: Conda (recommended)

For a fully working setup, use the provided Conda environment:

```
conda env create -f environment.yml
conda activate visualizebgw
```

Option 2: pip installation

You can also install the package directly from the project directory:
> If you plan to use only the Python library (without the GUI), install:
> ```
> pip install .
> ```
> For the full GUI and 3D visualization features, install with the optional extras:
> ```
> pip install visualizeBGW[gui]
> ```

## Usage

1. Using the Python API

    The library can be used directly in scripts or notebooks. All functionality is available through the top-level modules:
    * visualizeBGW.io – Readers for BerkeleyGW text and HDF5 files
    * visualizeBGW.analysis – Numerical processing (exciton components, absorption spectra, wavefunctions projection, etc.)
    * visualizeBGW.plotting – Matplotlib and PyVista visualizations.

    See [Doc](./documentation.md) for full available methods description and the attached [examples](./examples/).

2. Launching the GUI:

    After installation, run:
    ```
    visualizeBGW-gui
    ```

    This opens the multi-page interface where you can:
    * DFT and structure
      - View and export crystal structures
      - Plot and export wavefunction projections
    * GW
      - Explore band structures
      - Analyze $\chi$ and COHSEX convergence
      - Plot dielectric functions
    * BSE
      - Plot absorption spectra (with/without electron–hole interactions)
      - Inspect exciton eigenvectors and k-space components

    Each page includes file selectors, parameter inputs, and plotting buttons.
    3D views open in separate PyVista windows, ensuring the Qt event loop remains stable.

---

## Notes & Limitations

* visualizeBGW is tested on macOS and Linux with Python 3.9–3.13. Windows may work but PyVista/VTK support varies.
* HDF5 parsing assumes BerkeleyGW’s modern file schema. Older versions may require adjustments.
* Very large wavefunction grids can be memory-intensive, PyVista visualization is optimized but not memory-free.
* Multiprocessing-based 3D visualization may behave differently on Windows vs. Linux/macOS.

