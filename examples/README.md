# Quantum ESPRESSO and BerkeleyGW Workflow

This folder contains complete example workflows demonstrating how to prepare, run, and analyze calculations using **Quantum ESPRESSO (QE)** and **BerkeleyGW (BGW)**.

Each example directory includes:

- QE input files (`scf.in`, `wfn.in`, `pw2bgw.in`, etc.)
- BGW input files (`epsilon.inp`, `sigma.inp`, `kernel.inp`, `absorption.inp`)
- Shell scripts for running the calculations (optional — edit to match your machine)
- `notebook.ipynb` showing how to visualize the results using **visualizeBGW**

These folders provide a minimal example setup for running a QE–BGW workflow.
They serve as templates only and are not converged or production-level calculations.

For more details please refer BerkeleyGW [manual](http://manual.berkeleygw.org/4.0/overview-workflow/) or the [GitHub](https://github.com/BerkeleyGW/BerkeleyGW-examples) for more examples.
