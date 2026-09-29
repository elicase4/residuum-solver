# Example catalog

Every example lives in its own directory as `config.yaml`; a scenario with both 2D and 3D
variants splits into `2d/`/`3d/` subdirectories rather than a filename prefix.

## `mesh` — mesh generation and import

| Config | Demonstrates |
|---|---|
| `examples/mesh/square/config.yaml` | `BlockMesh2D` generator, quad elements |
| `examples/mesh/cube/config.yaml` | `BlockMesh3D` generator, hex elements |
| `examples/mesh/ubend/config.yaml` | Gmsh `.msh` import → `.pmsh` conversion, unstructured geometry |
| `examples/mesh/ubend_fine/config.yaml` | Same geometry as `ubend`, finer resolution |

## `heateq` — steady

| Config | Demonstrates |
|---|---|
| `steady/constant_conductivity/2d/config.yaml` | Baseline linear steady solve, 2D, constant isotropic conductivity |
| `steady/constant_conductivity/3d/config.yaml` | Same, 3D hex mesh |
| `steady/anisotropic_conductivity/2d/config.yaml` | Anisotropic conductivity tensor — a genuinely nonsymmetric operator, solved with GMRES (CG assumes SPD) |
| `steady/anisotropic_conductivity/3d/config.yaml` | Same, 3D |
| `steady/file_mode/config.yaml` | File-driven (not expression-driven) IC/BC/source fields (`.pndf` format) |
| `steady/gmsh_constant_conductivity/config.yaml` | Unstructured Gmsh-imported geometry, monitors (boundary flux integration) |
| `steady/ubend_nonlinear/config.yaml` | Nonlinear (temperature-dependent) conductivity + Newton + GMRES, fine unstructured geometry |

## `heateq` — transient

| Config | Demonstrates |
|---|---|
| `transient/heating_block/config.yaml` | Backward Euler time integration, VTU time series output ([tutorial](heating-block-tutorial.md)) |
| `transient/ubend_nonlinear/config.yaml` | Transient + nonlinear (temperature-dependent conductivity) + Newton + unstructured geometry — the most complete demo end to end |

## Running an example

Configs reference other files (meshes, field data) by path relative to the *current directory*,
not the config file's own location — `cd` into the example's directory before invoking the binary:

```bash
cd examples/mesh/square && /path/to/build/bin/mesh config.yaml && cd -
cd examples/heateq/steady/constant_conductivity/2d && /path/to/build/bin/heateq config.yaml && cd -
```

Output (VTK/VTU, monitor CSVs, solver/driver logs) lands in an `output/` directory next to the
config, which must already exist (the app does not create it) — every example directory here
already has one; git-ignores its contents (`examples/**/output/`).
