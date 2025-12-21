import numpy as np
import pyvista as pv
from scipy.spatial import Voronoi, ConvexHull
from ase.data import covalent_radii
from ase.data.colors import cpk_colors
from ..analysis.berkeleygw_processing import (
    p_in_kpts,
    p2p,
    get_kpts_in_plane,
    get_grid_and_levels,
)


def plot_structure(
    pl,
    structure,
    times_x=(0, 0),
    times_y=(0, 0),
    times_z=(0, 0),
    scaling=0.15,
    resolution=16,
):
    """
    Inputs:
        pl (pyvista.Plotter): PyVista plotter object to add the structure to.
        structure (dict): Dictionary containing structure information
        times_x (tuple): Number of times to repeat the structure in the negative and positive X directions - default is (0, 0).
        times_y (tuple): Number of times to repeat the structure in the negative and positive Y directions - default is (0, 0).
        times_z (tuple): Number of times to repeat the structure in the negative and positive Z directions - default is (0, 0).
        scaling (float): Scaling factor for atomic radii - default is 0.15.
        resolution (int): Resolution of the spheres representing atoms - default is 16.

    Outputs:
        Renders the atomic structure within the given PyVista plotter.

    Description:
        Displays the atomic structure using spheres for atoms and a transparent box for the unit cell.
    """
    nat = structure["atom_positions"].shape[0]

    a1, a2, a3 = structure["lattice_vectors"]
    cell_points = np.array(
        [[0.0, 0.0, 0.0], a1, a2, a3, a1 + a2, a1 + a3, a2 + a3, a1 + a2 + a3]
    )
    faces = np.hstack(
        [
            [4, 0, 1, 4, 2],
            [4, 3, 5, 7, 6],
            [4, 0, 1, 5, 3],
            [4, 0, 2, 6, 3],
            [4, 1, 4, 7, 5],
            [4, 2, 4, 7, 6],
        ]
    )
    cell_mesh = pv.PolyData(cell_points, faces)
    pl.add_mesh(
        cell_mesh, color="white", opacity=0.15, show_edges=True, edge_color="black"
    )

    try:
        per_atom_radii = scaling * covalent_radii[structure["atom_types"]]
    except (IndexError, TypeError, ValueError):
        per_atom_radii = np.full(nat, scaling, dtype=float)

    for Z, pos, r in zip(
        structure["atom_types"], structure["atom_positions"], per_atom_radii
    ):
        sphere = pv.Sphere(
            radius=float(r),
            center=pos,
            theta_resolution=resolution,
            phi_resolution=resolution,
        )
        for x in range(-times_x[0], times_x[1] + 1):
            for y in range(-times_y[0], times_y[1] + 1):
                for z in range(-times_z[0], times_z[1] + 1):
                    shift = x * a1 + y * a2 + z * a3
                    shifted_sphere = sphere.translate(shift, inplace=False)
                    pl.add_mesh(
                        shifted_sphere, color=cpk_colors[Z], smooth_shading=True
                    )


def plot_proj_wfn(
    pl,
    structure,
    rho,
    times_x=(0, 0),
    times_y=(0, 0),
    times_z=(0, 0),
    scaling=0.15,
    resolution=16,
):
    """
    Inputs:
        pl (pyvista.Plotter): PyVista plotter object to add the structure to.
        structure (dict): Dictionary containing structure information
        rho (ndarray): 3D array representing the projected wavefunction density.
        times_x (tuple): Number of times to repeat the structure in the negative and positive X directions - default is (0, 0).
        times_y (tuple): Number of times to repeat the structure in the negative and positive Y directions - default is (0, 0).
        times_z (tuple): Number of times to repeat the structure in the negative and positive Z directions - default is (0, 0).
        scaling (float): Scaling factor for atomic radii - default is 0.15.
        resolution (int): Resolution of the spheres representing atoms - default is 16.

    Outputs:
        Renders the atomic structure and isosurfaces of the projected wavefunction density within the given PyVista plotter.

    Description:
        Displays the atomic structure using spheres for atoms and isosurfaces representing the projected wavefunction density.
    """
    plot_structure(
        pl, structure, times_x, times_y, times_z, scaling=scaling, resolution=resolution
    )

    grid, levels = get_grid_and_levels(
        np.sqrt(rho), structure, n_levels=resolution, p_low=0.98, p_high=0.995
    )
    iso = grid.contour(isosurfaces=levels, scalars="rho")

    for x in range(-times_x[0], times_x[1] + 1):
        for y in range(-times_y[0], times_y[1] + 1):
            for z in range(-times_z[0], times_z[1] + 1):
                shift = (
                    x * structure["lattice_vectors"][0]
                    + y * structure["lattice_vectors"][1]
                    + z * structure["lattice_vectors"][2]
                )
                iso_shifted = iso.copy()
                iso_shifted.translate(shift, inplace=True)
                pl.add_mesh(
                    iso_shifted, opacity=0.1, cmap="Wistia", show_scalar_bar=False
                )


def plot_chi_converge(ax, chi_converge_data, G):
    """
    Inputs:
        ax (matplotlib.axes.Axes): Matplotlib Axes object to plot on.
        chi_converge_data (dict): Dictionary containing convergence data for chi.
        G (str): '0' for chi(0,0) or 'Gmax' for chi(Gmax,Gmax).

    Outputs:
        Plots the convergence of chi with respect to the number of conduction bands.

    Description:
        Generates a plot showing how chi converges as the number of conduction bands increases, including an extrapolated value.
    """
    ncbands = chi_converge_data["ncbands"]
    chi = chi_extrap = ylabel = None
    if G == "0":
        chi = chi_converge_data["chi(0,0)"]
        chi_extrap = chi_converge_data["chi(0,0)-extrap"]
        ylabel = r"$\chi(0,0)$"
    else:
        chi = chi_converge_data["chi(G_max,G_max)"]
        chi_extrap = chi_converge_data["chi(G_max,G_max)-extrap"]
        ylabel = r"$\chi(G_{\max},G_{\max})$"
    if chi is None or chi_extrap is None:
        raise ValueError("Invalid G specified for chi convergence plot.")

    ax.plot(ncbands, chi, "o-")
    ax.axhline(y=chi_extrap, color="r", linestyle="--")
    ax.set_xlabel("Number of Conduction Bands")
    ax.set_xlim([0, max(ncbands)])
    ax.set_ylabel(ylabel)
    ax.set_title(r"$\chi$ Convergence")
    ax.grid(True)
    ax.figure.tight_layout()


def plot_ch_convergence(ax, ch_converge_data, quantity):
    """
    Inputs:
        ax (matplotlib.axes.Axes): Matplotlib Axes object to plot on.
        ch_converge_data (dict): Dictionary containing convergence data for Coulomb-Hole term.
        quantity (str): 'vbm', 'cbm', or 'diff' to specify which quantity to plot.

    Outputs:
        Plots the convergence of the specified Coulomb-Hole quantity with respect to the number of bands.

    Description:
        Generates a plot showing how the Coulomb-Hole term converges as the number of bands increases, including an extrapolated value and error bars.
    """
    nbands = y = y_extrap = error = ylabel = None
    if quantity == "vbm":
        nbands = ch_converge_data["nbands"]
        y = ch_converge_data["CH(vbm)"]
        y_extrap = ch_converge_data["CH(vbm)-extrap"]
        error = ch_converge_data["CH(vbm)-error"]
        ylabel = r"$\sum_{CH}(vbm)$"
    elif quantity == "cbm":
        nbands = ch_converge_data["nbands"]
        y = ch_converge_data["CH(cbm)"]
        y_extrap = ch_converge_data["CH(cbm)-extrap"]
        error = ch_converge_data["CH(cbm)-error"]
        ylabel = r"$\sum_{CH}(cbm)$"
    elif quantity == "diff":
        nbands = ch_converge_data["nbands"]
        y = ch_converge_data["diff"]
        y_extrap = ch_converge_data["diff-extrap"]
        error = ch_converge_data["diff-error"]
        ylabel = r"$\sum_{CH}$"

    if nbands is None or y is None or y_extrap is None or error is None:
        raise ValueError(
            "Invalid quantity specified for Coulomb-Hole convergence plot."
        )

    ax.plot(nbands, y, "o-")
    ax.axhline(y=y_extrap, color="r", linestyle="--")
    ax.fill_between(
        [0, max(nbands)], y_extrap - error, y_extrap + error, color="r", alpha=0.2
    )
    ax.set_xlabel("Number of Bands")
    ax.set_xlim([0, max(nbands)])
    ax.set_ylabel(ylabel)
    ax.figure.tight_layout()


def plot_bandstructure(ax, kpts, emf, eqp, hsp, nv, use_emf=True, use_qp=True):
    """
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
    """
    hsp_labels = [hsp[0][0]]
    hsp_points = np.array([0], dtype=int)
    path = np.array([], dtype=int)
    for i in range(1, len(hsp)):
        hsp_labels.append(hsp[i][0])
        point = hsp[i][1]
        point0 = hsp[i - 1][1]
        _ = p_in_kpts(kpts, point0)
        segment = p2p(kpts, point0, point)[:-1]
        path = np.concatenate((path, segment))
        hsp_points = np.concatenate((hsp_points, [len(path)]))
    path = np.concatenate((path, [p_in_kpts(kpts, hsp[-1][1])]))
    if use_emf:
        ax.plot(
            np.arange(len(path)),
            emf[path, :] - np.max(emf[:, nv - 1]),
            color="gray",
            linestyle="--",
            label="MF",
        )
    if use_qp:
        ax.plot(
            np.arange(len(path)),
            eqp[path, :] - np.max(eqp[:, nv - 1]),
            color="b",
            label="GW QP",
        )
    if use_emf and use_qp:
        handles, labels = ax.get_legend_handles_labels()
        unique = {}
        for h, l in zip(handles, labels):
            if l not in unique:
                unique[l] = h
        ax.legend(unique.values(), unique.keys())
    ax.set_xticks(hsp_points)
    ax.set_xticklabels(hsp_labels)
    ax.set_xlabel("$k$-Path")
    ax.set_xlim([0, len(path) - 1])
    ax.set_ylabel("$E - E_{vbm}$ (eV)")
    ax.grid(True)
    ax.figure.tight_layout()


def plot_eps_head(ax, qlen, epshead, plot_inverse=True):
    """
    Inputs:
        ax (matplotlib.axes.Axes): Matplotlib Axes object to plot on.
        qlen (ndarray): Array of |q| values.
        epshead (ndarray): Array of epsilon head values.
        plot_inverse (bool): Whether to plot the inverse of epsilon - default is True.

    Outputs:
        Plots the static dielectric function or its inverse as a function of |q|.

    Description:
        Generates a plot showing either the static dielectric function eps(q) or its inverse eps^{-1}(q) against the magnitude of the wavevector |q|.
    """
    if plot_inverse:
        ax.plot(qlen, epshead, marker="o", linestyle="-")
        ax.set_ylabel(r"$\epsilon_{(0,0)}^{-1}(q)$")
    else:
        ax.plot(qlen, 1 / epshead, marker="o", linestyle="-")
        ax.set_ylabel(r"$\epsilon_{(0,0)}(q)$")
    ax.set_xlabel("|q|")
    ax.set_title("Static Dielectric Function vs |q|")
    ax.grid(True)
    ax.figure.tight_layout()


def plot_absorption_spectra(
    ax, energy, spectra, xmin, xmax, noeh_energy=None, noeh_spectra=None
):
    """
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
    """
    ax.plot(energy, spectra)
    if noeh_spectra is not None and noeh_energy is not None:
        ax.plot(noeh_energy, noeh_spectra, label="w/o e-h interactions", linestyle="--")
        ax.legend()
    ax.set_xlabel("Energy (eV)")
    if xmin is None or xmax is None:
        xmin0, xmax0 = ax.get_xlim()
        xmin = xmin if xmin is not None else xmin0
        xmax = xmax if xmax is not None else xmax0
    ax.set_xlim([xmin, xmax])
    ax.set_xlim([xmin, xmax])
    ax.set_ylabel("Absorption Spectra (arb. units)")
    ax.set_title("BSE Absorption Spectra")
    ax.grid(True)
    ax.figure.tight_layout()


def plot_eigenvector_components(
    ax, index, Ak, Ac, Av, exciton_header, mf_header, eigs, dipole, scaling=300.0
):
    """
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
    """
    nc = exciton_header["params"]["nc"]
    nv = exciton_header["params"]["nv"]
    kpts = exciton_header["kpoints"]["kpts"]
    Qshift = exciton_header["kpoints"]["exciton_Q_shifts"]
    bvec = mf_header["crystal"]["bvec"]
    kgrid = mf_header["kpoints"]["kgrid"]

    dims = [dim for dim in [0, 1, 2] if dim != 2 - np.argmin(np.flip(kgrid))]
    kpts_e = np.dot(kpts, bvec)
    kpts_e_2d, Ak_e_2d = get_kpts_in_plane(kpts_e, Ak[0, :], dims)

    kpts_h = np.dot(kpts - Qshift[0, :], bvec)
    kpts_h_2d, Ak_h_2d = get_kpts_in_plane(kpts_h, Ak[0, :], dims)

    fig = ax.figure
    ax.set_visible(False)
    gs = ax.get_subplotspec().subgridspec(1, 2, width_ratios=[1.0, 0.25], wspace=0.35)

    axk = fig.add_subplot(gs[0])
    axk.scatter(kpts_e_2d[:, 0], kpts_e_2d[:, 1], s=1.0, color="gray", alpha=0.25)
    axk.scatter(
        kpts_e_2d[:, 0], kpts_e_2d[:, 1], s=Ak_e_2d * scaling, color="red", alpha=0.5
    )
    axk.scatter(
        kpts_h_2d[:, 0], kpts_h_2d[:, 1], s=Ak_h_2d * scaling, color="blue", alpha=0.5
    )
    if 0 in dims:
        axk.set_xlabel("$k_x$")
        if 1 in dims:
            axk.set_ylabel("$k_y$")
        else:
            axk.set_ylabel("$k_z$")
    else:
        axk.set_xlabel("$k_y$")
        axk.set_ylabel("$k_z$")

    axb = fig.add_subplot(gs[1])
    y_valence = np.arange(nv)
    y_conduction = nv + 1.5 + np.arange(nc)
    valence_labels = []
    for j in range(nv):
        offset = nv - 1 - j
        valence_labels.append("VB" if offset == 0 else f"VB-{offset}")
    conduction_labels = []
    for i in range(nc):
        conduction_labels.append("CB" if i == 0 else f"CB+{i}")

    for y in y_valence:
        axb.hlines(y, 0.0, 1.0, color="black", linewidth=1)
    for y in y_conduction:
        axb.hlines(y, 0.0, 1.0, color="black", linewidth=1)
    bands_scale = 300.0
    axb.scatter(
        np.full(nv, 0.5),
        np.flip(y_valence),
        s=Av[0, :] * bands_scale,
        color="blue",
        alpha=0.75,
    )
    axb.scatter(
        np.full(nc, 0.5), y_conduction, s=Ac[0, :] * bands_scale, color="red", alpha=0.75
    )

    yticks = np.concatenate([y_valence, y_conduction])
    axb.set_yticks(yticks)
    axb.set_yticklabels(valence_labels + conduction_labels)
    axb.set_xticks([])
    axb.set_xlim(0.0, 1.0)
    axb.set_ylim(-0.5, y_conduction[-1] + 0.5)
    axb.spines["top"].set_visible(False)
    axb.spines["right"].set_visible(False)
    axb.spines["bottom"].set_visible(False)
    axb.spines["left"].set_visible(False)
    axb.set_ylabel("Bands")

    axk.set_title(f"Exciton Index {index}, {eigs:.3f} eV - O.S. {np.abs(dipole):.3f}")
    fig.tight_layout()


def plot_eigenvector_k_3d(
    pl, Ak, exciton_header, mf_header, scaling=0.05, resolution=8
):
    """
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
    """
    kpts = exciton_header["kpoints"]["kpts"]
    Qshift = exciton_header["kpoints"]["exciton_Q_shifts"]
    bvec = mf_header["crystal"]["bvec"]

    coords = np.array(
        [
            i * bvec[0] + j * bvec[1] + k * bvec[2]
            for i in range(-1, 2)
            for j in range(-1, 2)
            for k in range(-1, 2)
        ]
    )
    origin_idx = np.where(np.all(np.isclose(coords, 0.0), axis=1))[0][0]
    vor = Voronoi(coords)
    region_index = vor.point_region[origin_idx]
    region = [v for v in vor.regions[region_index] if v != -1]
    bz_vertices = vor.vertices[region]
    hull = ConvexHull(bz_vertices)
    faces = np.hstack([[3] + list(tri) for tri in hull.simplices])
    bz_mesh = pv.PolyData(bz_vertices, faces)
    pl.add_mesh(bz_mesh, color="white", opacity=0.2, show_edges=False)

    radii = Ak[0, :] / np.max(Ak[0, :]) * scaling + (scaling * 0.1)
    base_sphere = pv.Sphere(
        radius=1.0, theta_resolution=resolution, phi_resolution=resolution
    )

    kpts_e = kpts @ bvec
    e_points = kpts_e
    e_cloud = pv.PolyData(e_points)
    e_cloud["radius"] = radii
    e_glyphs = e_cloud.glyph(geom=base_sphere, scale="radius", factor=1.0, orient=False)
    pl.add_mesh(
        e_glyphs, scalars="radius", opacity=0.5, cmap="Reds", show_scalar_bar=False
    )

    kpts_h = (kpts - Qshift[0, :]) @ bvec
    h_points = kpts_h
    h_cloud = pv.PolyData(h_points)
    h_cloud["radius"] = radii
    h_glyphs = h_cloud.glyph(geom=base_sphere, scale="radius", factor=1.0, orient=False)
    pl.add_mesh(
        h_glyphs, scalars="radius", opacity=0.5, cmap="Blues", show_scalar_bar=False
    )

    pl.add_axes()
