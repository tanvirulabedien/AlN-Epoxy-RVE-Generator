# AlN/Epoxy RVE Geometry Generator for ANSYS FEM

A Python tool that generates **random close-packed (RCP) Representative Volume Element (RVE)** geometries of AlN particles in an Epoxy matrix, exported as named multi-body **STEP files** ready for direct import into **ANSYS DesignModeler**.

Developed for FEM-based thermal conductivity simulation of particulate polymer composites.

---

## Overview

Predicting the effective thermal conductivity of particle-filled polymer composites using FEM requires a realistic microstructure model. This tool automates the full geometry generation pipeline:

1. Computes RVE box dimensions based on target volume fraction and particle size
2. Places spheres randomly using **Random Sequential Addition (RSA)** or **Lubachevsky-Stillinger (LS) growth**
3. Performs Boolean subtraction to create the matrix body
4. Exports a **multi-body STEP file** with bodies pre-named for ANSYS

---

## Features

- Supports **three particle sizes**: 1 µm, 5 µm, and 80 µm diameter AlN particles
- Supports **four volume fractions**: 10%, 20%, 30%, and 40%
- RSA algorithm for VF ≤ 30% (fast and reliable)
- Lubachevsky-Stillinger growth algorithm for VF = 40% (overcomes RSA jamming limit)
- Bodies are **pre-labeled** in the STEP file:
  - `Epoxy_matrix` → the cube with spherical holes (matrix)
  - `AlN_particle_001`, `AlN_particle_002`, ... → individual particles
- Minimum 0.5 µm surface-to-surface gap between all spheres (prevents meshing failures)

---

## Generated Files

| File | VF | Spheres | RVE Size | Algorithm |
|---|---|---|---|---|
| `AlN_Epoxy_VF10_named.step` | 10% | 20 | 377.1 µm | RSA |
| `AlN_Epoxy_VF20_named.step` | 20% | 20 | 299.3 µm | RSA |
| `AlN_Epoxy_VF30_named.step` | 30% | 60 | 377.1 µm | RSA |
| `AlN_Epoxy_VF40_named.step` | 40% | 40 | 299.3 µm | LS growth |

---

## Requirements

- Python 3.8 or higher
- numpy
- gmsh

Install dependencies:

```bash
pip install numpy gmsh
```

---

## Usage

```bash
python generate_rcp_geometry.py
```

Output STEP files are written to the `results/` folder.

---

## How to Import into ANSYS DesignModeler

1. Open ANSYS Workbench → Geometry cell → DesignModeler
2. `File → Import External Geometry File` → select the `.step` file
3. Click **Generate**
4. In the model tree, bodies appear as `Epoxy_matrix` and `AlN_particle_001...NNN`
5. Shift-click to select all `AlN_particle_*` bodies → assign **AlN** material
6. Select `Epoxy_matrix` → assign **Epoxy** material
7. Right-click geometry → **Shared Topology → Share** (critical for conformal mesh at interfaces)

---

## Material Properties Used

| Material | Thermal Conductivity | Density | Specific Heat |
|---|---|---|---|
| Epoxy (matrix) | 0.20 W/(m·K) | 1200 kg/m³ | 1100 J/(kg·K) |
| AlN (particles) | 180 W/(m·K) | 3260 kg/m³ | 740 J/(kg·K) |

---

## Algorithm Details

### Random Sequential Addition (RSA)

Used for VF = 10%, 20%, and 30%.

```
1. Calculate box size L from target VF and number of spheres
2. Pick a random (x, y, z) point inside the box
3. Check: does this sphere overlap any placed sphere?
          surface-to-surface distance < 0.5 µm → reject
4. Check: is the sphere fully inside the box?
          center closer than R to any wall → reject
5. If both pass → place sphere
6. Repeat until target VF is reached
```

### Lubachevsky-Stillinger (LS) Growth

Used for VF = 40% because RSA jams at ~38% for monodisperse spheres.

```
1. Place all sphere CENTERS randomly (radius = 1 µm, no overlaps)
2. Grow all spheres simultaneously toward target radius (40 µm)
3. At each growth step, push apart any overlapping spheres
4. Stop when all spheres reach full size
```

### Boolean Subtraction (gmsh + OpenCASCADE)

```
RVE Box  −  All Spheres  =  Matrix body (with holes)
                           + Sphere bodies (kept separate)
```

The OpenCASCADE kernel (same as used in Salome and FreeCAD) ensures watertight conformal geometry.

---

## Validation

FEM results from this geometry can be compared against established Effective Medium Approximation (EMA) models:

- **Maxwell model** — accurate at low VF (< 20%)
- **Hashin-Shtrikman bounds** — rigorous upper and lower bounds
- **Lewis-Nielsen model** — accounts for random close packing limit (φ_max = 0.64)
- **Bruggeman model** — self-consistent, suitable for medium-to-high VF

### Recommended References

- Ngo, I-L. et al. (2017). *A modified Hashin-Shtrikman model for predicting thermal conductivity of polymer composites with randomly distributed hybrid fillers.* Composites Part B.
- Kumari, P. et al. (2010). *A computational and experimental investigation on thermal conductivity of particle reinforced epoxy composites.* Computational Materials Science. (**uses ANSYS, same setup**)
- Progelhof, R.C. et al. (1976). *Methods for predicting the thermal conductivity of composite systems.* Polymer Engineering & Science.

---

## Project Structure

```
AlN-Epoxy-RVE-Generator/
│
├── generate_rcp_geometry.py     # Main geometry generation script
├── README.md                    # This file
├── requirements.txt             # Python dependencies
│
├── results/                     # Generated STEP files
│   ├── AlN_Epoxy_VF10_named.step
│   ├── AlN_Epoxy_VF20_named.step
│   ├── AlN_Epoxy_VF30_named.step
│   └── AlN_Epoxy_VF40_named.step
│
└── figures/                     # Geometry visualizations (optional)
```

---

## How to Cite

If you use this code in your research, please cite:

```
[Your Name] ([Year]). AlN/Epoxy RVE Geometry Generator for ANSYS FEM.
GitHub. https://github.com/YOURUSERNAME/AlN-Epoxy-RVE-Generator
```

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## Contact

For questions or collaboration, please open an **Issue** on this repository.
