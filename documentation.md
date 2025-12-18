# documentation.md

API Reference for visualizeBGW

---

## visualizeBGW.io

### io.berkeleygw_readers

#### read_h5_group(group)
```
Input:
    group (h5py.Group) — HDF5 group object

Output:
    Nested dictionary of NumPy arrays

Description:
    Recursively extracts a full HDF5 group into a Python dictionary.
    ! Note: Use only for small to medium-sized datasets! Useful for reading metadata.
```

#### read_structure(path)
```
Input:
    path (str) — path to WFN.h5

Output:
    structure (dict) — alat, lattice vectors, atomic species, fractional positions
    mf_header (dict) — parsed metadata from the HDF5 file

Description:
    Reads structural and metadata information from BerkeleyGW’s WFN.h5.
```

#### read_chi_converge(path)
```
Input:
    path (str) — path to chi_converge.dat

Output:
    dict of convergence data per q-point

Description:
    Parses BerkeleyGW chi_converge.dat and returns convergence information.
```

#### read_ch_converge(path)
```
Input:
    path (str) — path to ch_converge.dat

Output:
    dict with Coulomb-hole term vs number of bands

Description:
    Reads the COHSEX Coulomb-hole convergence output file.
```

#### read_eqp(path)
```
Input:
    path (str) — path to eqp.dat

Output:
    kpts (np.ndarray) — array of k-points
    emf (np.ndarray) — array of mean-field energies
    eqp (np.ndarray) — array of quasiparticle energies

Description:
    Reads BerkeleyGW eqp.dat file and extracts k-points, mean-field energies, and quasiparticle energies.
```

#### read_eigenvalues_noeh(path)
```
Input:
    path (str) — path to eigenvalues_noeh.dat

Output:
    eigs (np.ndarray) — array of eigenvalues without electron-hole interaction
    dipole (np.ndarray) — array of complex dipole matrix elements

Description:
    Reads BerkeleyGW eigenvalues_noeh.dat file and extracts eigenvalues and dipole matrix elements.
```

#### read_eigenvalues(path)
```
Input:
    path (str) — path to eigenvalues.dat

Output:
    eigs (np.ndarray) — array of eigenvalues with electron-hole interaction
    dipole (np.ndarray) — array of complex dipole matrix elements

Description:
    Reads BerkeleyGW eigenvalues.dat file and extracts eigenvalues and dipole matrix elements.
```

### io.berkeleygw_exports

#### export_structure(path, structure)
```
Input:
    path (str) — directory
    structure (dict) — lattice & atoms

Output:
    Writes structure.vasp in POSCAR format

Description:
    Exports structural geometry to a VASP-compatible file.
```

#### export_xsf(xsf_dir, rho, mf_header)
```
Input:
    xsf_dir (str) — directory
    rho (ndarray) — scalar field
    mf_header (dict) — structural metadata

Output:
    Writes projwfn.xsf file

Description:
    Saves the 3D real-space wavefunction density into XSF format for XCrySDen.
```

--- 

## visualizeBGW.analysis

### analysis.berkeleygw_processing

#### p_in_kpts(kpts, p, tol)
```
Input:
    kpts (np.ndarray): of shape (N, 3), list of k-points
    p (list or np.ndarray): of shape (3,), point to find in kpts
    tol (float): tolerance for matching points - default 1e-6

Output:
    idx (int): index of point p in kpts
```

#### p2p(kpts, p1, p2)
```
Input:
    kpts (np.ndarray): of shape (N, 3), list of k-points
    p1 (list or np.ndarray): of shape (3,), starting point of the line segment
    p2 (list or np.ndarray): of shape (3,), ending point of the line segment
    tol (float): tolerance for distance to the line segment - default 1e-6

Output:
    p2p (np.ndarray): of shape (M,), indices of k-points that lie on the line segment from p1 to p2, ordered along the segment
```

#### get_kpts_in_plane(kpts, quantity, dims):
```
Input:
    kpts (np.ndarray): of shape (N, 3), list of k-points
    quantity (np.ndarray): of shape (N,), quantity defined at each k-point
    dims (list or np.ndarray): of shape (2,), indices of the two dimensions defining the plane (e.g., [0, 1] for kx-ky plane)

Output:
    kpts_in_plane (np.ndarray): of shape (M, 2), unique k-points projected onto the specified 2D plane
    quantity_in_plane (np.ndarray): of shape (M,), summed quantity at each unique k-point in the plane

Description:
    Projects k-points onto the specified 2D plane and sums the quantity for k-points that map to the same point in that plane.
```

#### calc_rho(wfn_file, bnd, kpt)
```
Input:
    wfn_file (str): path to the wavefunction HDF5 file
    bnd (int): band index
    kpt (list or np.ndarray): of shape (3,), k-point at which to calculate the charge density
    
Output:
    rho (np.ndarray): charge density grid
    mf_header (dict): metadata from the wavefunction file

Description:
    Computes the mean-field wavefunction density on a 3D grid in real space.
```

#### get_grid_and_levels(rho, structure, n_levels, p_low, p_high)
```
Input:
    rho (np.ndarray): 3D charge density grid
    structure (dict): structure information containing lattice vectors
    n_levels (int): number of isosurface levels to generate - default 1
    p_low (float): lower percentile for isosurface levels - default 0.98
    p_high (float): upper percentile for isosurface levels - default 0.995

Output:
    grid (pyvista.StructuredGrid): structured grid with charge density data
    levels (np.ndarray): array of isosurface levels

Description:
    Creates a structured grid from the charge density and computes isosurface levels based on specified percentiles.
```

#### get_eps_head(epsmat_file):
```
Input:
    epsmat_file (str): path to the epsilon matrix HDF5 file

Output:
    epshead (list): head of the dielectric function at each q-point
    qlen (list): length of each q-point in the bdot metric

Description:
    Extracts the head of the dielectric function and the lengths of q-points from the epsilon
```

#### get_absorption_spectra(eigs, dipole, broadening, function, de):
```
Input:
    eigs (np.ndarray): of shape (N,), excitation energies
    dipole (np.ndarray): of shape (N,), dipole matrix elements
    broadening (float): broadening parameter
    function (str): type of broadening function ('lorentzian', 'voigt', or 'gaussian')
    de (float): energy step for the spectra - default 0.001 eV

Output:
    energy (np.ndarray): energy grid
    spectra (np.ndarray): normalized absorption spectra

Description:
    Computes the absorption spectra using the specified broadening function.
```

#### get_eigenvectors_components(eigenvectors_file, index):
```
Input:
    eigenvectors_file (str): path to the eigenvectors HDF5 file
    index (int): index of the exciton state to analyze

Output:
    Ak (np.ndarray): of shape (Nk,), k-point contribution to the exciton state
    Ac (np.ndarray): of shape (Nc,), conduction band contribution to the exciton state
    Av (np.ndarray): of shape (Nv,), valence band contribution to the exciton state
    exciton_header (dict): metadata from the exciton file
    mf_header (dict): metadata from the mean-field file

Description:
    Decomposes the exciton eigenvector into contributions from k-points, conduction bands, and valence bands.
```

---

## visualizeBGW.plotting

### plotting.berkeleygw_plots

#### plot_structure(pl, structure, times_X, times_Y, times_Z, scaling, resolution):
```
Inputs:
    pl (pyvista.Plotter): PyVista plotter object to add the structure to.
    structure (dict): Dictionary containing structure information
    times_X (tuple): Number of times to repeat the structure in the negative and positive X directions - default is (0, 0).
    times_Y (tuple): Number of times to repeat the structure in the negative and positive Y directions - default is (0, 0).
    times_Z (tuple): Number of times to repeat the structure in the negative and positive Z directions - default is (0, 0).
    scaling (float): Scaling factor for atomic radii - default is 0.15.
    resolution (int): Resolution of the spheres representing atoms - default is 16.

Outputs:
    Renders the atomic structure within the given PyVista plotter.

Description:
    Displays the atomic structure using spheres for atoms and a transparent box for the unit cell.
```

#### plot_proj_wfn(pl, structure, rho, times_X, times_Y, times_Z, scaling, resolution)
```
Inputs:
    pl (pyvista.Plotter): PyVista plotter object to add the structure to.
    structure (dict): Dictionary containing structure information
    rho (ndarray): 3D array representing the projected wavefunction density.
    times_X (tuple): Number of times to repeat the structure in the negative and positive X directions - default is (0, 0).
    times_Y (tuple): Number of times to repeat the structure in the negative and positive Y directions - default is (0, 0).
    times_Z (tuple): Number of times to repeat the structure in the negative and positive Z directions - default is (0, 0).
    scaling (float): Scaling factor for atomic radii - default is 0.15.
    resolution (int): Resolution of the spheres representing atoms - default is 16.

Outputs:
    Renders the atomic structure and isosurfaces of the projected wavefunction density within the given PyVista plotter.

Description:
    Displays the atomic structure using spheres for atoms and isosurfaces representing the projected wavefunction density.
```

#### plot_chi_converge(ax, chi_converge_data, G)
```
Inputs:
    ax (matplotlib.axes.Axes): Matplotlib Axes object to plot on.
    chi_converge_data (dict): Dictionary containing convergence data for chi.
    G (str): '0' for chi(0,0) or 'Gmax' for chi(Gmax,Gmax).

Outputs:
    Plots the convergence of chi with respect to the number of conduction bands.

Description:
    Generates a plot showing how chi converges as the number of conduction bands increases, including an extrapolated value.
```

#### plot_ch_convergence(ax, ch_converge_data, quantity)
```
Inputs:
    ax (matplotlib.axes.Axes): Matplotlib Axes object to plot on.
    ch_converge_data (dict): Dictionary containing convergence data for Coulomb-Hole term.
    quantity (str): 'vbm', 'cbm', or 'diff' to specify which quantity to plot.

Outputs:
    Plots the convergence of the specified Coulomb-Hole quantity with respect to the number of bands.

Description:
    Generates a plot showing how the Coulomb-Hole term converges as the number of bands increases, including an extrapolated value and error bars.
```

#### plot_bandstructure(ax, kpts, emf, eqp, hsp, nv, use_emf, use_qp)
```
Inputs:
    ax (matplotlib.axes.Axes): Matplotlib Axes object to plot on.
    kpts (ndarray): Array of k-points.
    emf (ndarray): Mean-field energies.
    eqp (ndarray): Quasiparticle energies.
    hsp (list): list of High-symmetry points labels and their coordinates.
    nv (int): Number of valence bands.
    use_emf (bool): Whether to plot mean-field energies - default is True.
    use_qp (bool): Whether to plot quasiparticle energies - default is True.

Outputs:
    Plots the band structure along the specified k-point path.

Description:
    Generates a band structure plot showing mean-field and/or quasiparticle energies along a defined k-point path, with high-symmetry point labels.
```

#### plot_eps_head(ax, qlen, epshead, plot_inverse)
```
Inputs:
    ax (matplotlib.axes.Axes): Matplotlib Axes object to plot on.
    qlen (ndarray): Array of |q| values.
    epshead (ndarray): Array of epsilon head values.
    plot_inverse (bool): Whether to plot the inverse of epsilon - default is True.

Outputs:
    Plots the static dielectric function or its inverse as a function of |q|.

Description:
    Generates a plot showing either the static dielectric function eps(q) or its inverse eps^{-1}(q) against the magnitude of the wavevector |q|.
```

#### plot_absorption_spectra(ax, energy, spectra, xmin, xmax, noeh_energy, noeh_spectra)
```
Inputs:
    ax (matplotlib.axes.Axes): Matplotlib Axes object to plot on.
    energy (ndarray): Array of energy values.
    spectra (ndarray): Array of absorption spectra values.
    xmin (float): Minimum x-axis limit.
    xmax (float): Maximum x-axis limit.
    noeh_energy (ndarray): Array of energy values without electron-hole interactions - default is None.
    noeh_spectra (ndarray): Array of absorption spectra values without electron-hole interactions - default is None.

Outputs:
    Plots the absorption spectra with and without electron-hole interactions.

Description:
    Generates a plot showing the absorption spectra as a function of energy, optionally including a comparison without electron-hole interactions.
```

#### plot_eigenvector_components(ax, index, Ak, Ac, Av, exciton_header, mf_header, eigs, dipole, scaling)
```
Inputs:
    ax (matplotlib.axes.Axes): Matplotlib Axes object to plot on.
    index (int): Exciton index.
    Ak (ndarray): Array of exciton eigenvector amplitudes.
    Ac (ndarray): Array of conduction band contributions.
    Av (ndarray): Array of valence band contributions.
    exciton_header (dict): Dictionary containing exciton calculation parameters.
    mf_header (dict): Dictionary containing mean-field calculation parameters.
    eigs (float): Exciton eigenvalue.
    dipole (ndarray): Exciton dipole moment.
    scaling (float): Scaling factor for plotting - default is 300.0.

Outputs:
    Plots the exciton eigenvector components in k-space and band contributions.

Description:
    Generates a plot showing the exciton eigenvector components in k-space along with the contributions from valence and conduction bands.
```

#### plot_eigenvector_k_3D(pl, Ak, exciton_header, mf_header, scaling, resolution):
```
Inputs:
    pl (pyvista.Plotter): PyVista plotter object to add the 3D plot to.
    Ak (ndarray): Array of exciton eigenvector amplitudes.
    exciton_header (dict): Dictionary containing exciton calculation parameters.
    mf_header (dict): Dictionary containing mean-field calculation parameters.
    scaling (float): Scaling factor for the size of the k-point spheres - default is 0.05.
    resolution (int): Resolution of the spheres representing k-points - default is 8.

Outputs:
    Renders a 3D plot of the exciton eigenvector components in k-space within the given PyVista plotter.

Description:
    Displays the exciton eigenvector components in k-space using spheres for k-points, colored by electron and hole contributions, along with the first Brillouin zone.
```

---

## visualizeBGW.gui

### MainWindow

Entry point for the GUI.
Manages page navigation using a QStackedWidget.

### pages/*

Each page corresponds to one visualization module:
- CrystalStructurePage
- WavefunctionProjectionPage
- ChiConvergencePage
- CHConvergencePage
- BandstructurePage
- DielectricFunctionPage
- AbsorptionSpectraPage
- ExcitonComponentsPage

Each page provides:
- File selectors
- Validation
- Plotting actions
- Calls to analysis and plotting modules
- Optional PyVista window spawning
- 