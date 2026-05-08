# SeaMotions

SeaMotions is a Boundary Element Method (BEM) solver that models the interaction between ocean waves and floating bodies across a wide frequency range. The codebase serves both research and industrial design workflows by combining high-order hydrodynamics, multi-body formulations, and curated tooling to inspect loads, motions, and derived quantities.

## Capabilities

### Hydrodynamics Core
- BEM formulation for mono- and multi-body layouts, including Oscillating Water Columns (OWCs).
- Database generators for deep- and finite-depth conditions.
- Consistent hydrostatic restoring and initial-stability estimators that feed the dynamic solvers.

### First-Order Frequency-Domain Analysis
- Added-mass and radiation damping matrices for 6N rigid-body modes.
- Incident-wave diffraction, Froude-Krylov loads, and complex wave exciting forces.
- Field reconstruction for velocities, pressures, and free-surface elevations throughout the computational domain.

### Second-Order and Nonlinear Estimators
- Quadratic Transfer Functions (QTFs) via Pinkster and X. Chen (indirect) formulations.
- Mean-drift force prediction built from first-order solutions.

### Numerical Infrastructure and Tooling
- MPI-enabled execution paths for heavy frequency sweep workloads.
- Python-based utilities in `aux_tools/` to create, validate, and fit hydrodynamic databases.
- Ready-to-run examples under `examples/` and regression tests in `tests/` for continuous validation.

## Repository Layout
- `src/`: Core C++ solvers, interfaces, and numerical utilities.
- `examples/`: Representative hydrodynamic setups with ready-made meshes.
- `aux_tools/`: Python scripts to generate integrals, fit coefficient surfaces, and cross-check outputs.
- `tests/`: Catch2-based unit and integration tests guarding dispersion, hydrostatics, MPI, and IO behavior.
- `bin/`: Convenience launchers produced after building (`seamotions_freq`, `seamotions_stab`, etc.).

## Documentation Roadmap
- Build and install instructions: see [INSTALL.md](INSTALL.md).
- Parallel HDF5 compilation guide: see [PARALLEL_HDF5.md](PARALLEL_HDF5.md).
- Example setups and validation workflows: inspect the `examples/` and `tests/` folders.

## Command Line
Run the frequency solver by passing the case folder path. The RAM usage report is optional and writes CSV files into the case folder.

Example:
```bash
seamotions_freq --mem-report <case_folder>
```

Options:
- `-c, --case <path>`: case folder path (optional when passing as a positional argument)
- `--mem-report`: write RAM usage CSV reports
- `-h, --help`: show help

## Contact
For any software-related questions or bug reports, contact [Sergio Fernandez Ruano](https://ihcantabria.com/directorio-personal/sergio-fernandez-ruano/).