"""
================================================================================
AlN/Epoxy RVE Geometry Generator for ANSYS FEM
================================================================================
Generates random close-packed (RCP) Representative Volume Element (RVE)
geometries of AlN particles in an Epoxy matrix.

Exports named multi-body STEP files directly importable into ANSYS DesignModeler.
Bodies are pre-labeled:
    - "Epoxy_matrix"        → cube with spherical holes (matrix)
    - "AlN_particle_001"    → first AlN sphere
    - "AlN_particle_002"    → second AlN sphere
    - ...

Algorithms:
    - Random Sequential Addition (RSA)      → used for VF = 10%, 20%, 30%
    - Lubachevsky-Stillinger (LS) growth    → used for VF = 40%

Usage:
    python generate_rcp_geometry.py

Requirements:
    pip install numpy gmsh

Output:
    results/AlN_Epoxy_VF10_named.step
    results/AlN_Epoxy_VF20_named.step
    results/AlN_Epoxy_VF30_named.step
    results/AlN_Epoxy_VF40_named.step
    results/sphere_coordinates_all_VF.csv

Author : [Your Name]
Date   : 2024
License: MIT
================================================================================
"""

import numpy as np
import gmsh
import os
import re
import csv

# ── Configuration ─────────────────────────────────────────────────────────────

# Particle sizes (radius in µm) — change this to switch between 1 µm, 5 µm, 80 µm
PARTICLE_DIAMETER_UM = 80.0          # µm  (change to 1.0 or 5.0 for other sizes)
R_NOM = PARTICLE_DIAMETER_UM / 2     # radius in µm

# Volume fractions to generate
VOLUME_FRACTIONS = [0.10, 0.20, 0.30, 0.40]

# Target number of spheres per RVE (20 is statistically sufficient)
# For VF = 30% we use 60 spheres (larger RVE) to help RSA reach 30%
TARGET_N = {
    0.10: 20,
    0.20: 20,
    0.30: 60,
    0.40: 40,
}

# Minimum surface-to-surface gap between spheres (µm)
# Prevents mesh failures at touching faces
MIN_GAP = 0.5

# Wall clearance: sphere centers kept at least this far from box walls (µm)
WALL_CLEARANCE = R_NOM + 0.3

# Output directory
OUTPUT_DIR = "results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Random seeds (one per VF for reproducibility)
SEEDS = {0.10: 101, 0.20: 202, 0.30: 303, 0.40: 404}


# ── Helper: sphere volume ─────────────────────────────────────────────────────

def sphere_volume(r):
    """Volume of a sphere with radius r."""
    return (4/3) * np.pi * r**3


# ── Packing: RSA ──────────────────────────────────────────────────────────────

def rsa_pack(vf, n_target, r, seed):
    """
    Random Sequential Addition (RSA) packing.

    Places spheres one at a time at random positions, rejecting any that
    overlap an existing sphere or breach the wall clearance.

    Parameters
    ----------
    vf       : float  — target volume fraction (e.g. 0.10)
    n_target : int    — target number of spheres
    r        : float  — sphere radius (µm)
    seed     : int    — random seed for reproducibility

    Returns
    -------
    centers  : np.ndarray, shape (N, 3) — sphere center coordinates (µm)
    L        : float                    — RVE box side length (µm)
    achieved : float                    — actual achieved volume fraction
    """
    rng     = np.random.default_rng(seed)
    min_d   = 2 * r + MIN_GAP          # minimum center-to-center distance
    wall    = r + 0.3                  # wall clearance for centers
    L       = (n_target * sphere_volume(r) / vf) ** (1/3)  # box size
    lo, hi  = wall, L - wall

    centers = []

    for _ in range(5_000_000):         # max attempts
        if len(centers) >= n_target:
            break
        p = rng.uniform(lo, hi, 3)    # random candidate position

        # Check overlap with all placed spheres
        if centers:
            arr  = np.array(centers)
            dist2 = np.sum((arr - p) ** 2, axis=1)
            if np.any(dist2 < min_d ** 2):
                continue               # overlap → reject

        centers.append(p.tolist())    # no overlap → accept

    centers  = np.array(centers) if centers else np.empty((0, 3))
    achieved = len(centers) * sphere_volume(r) / L ** 3
    return centers, L, achieved


# ── Packing: Lubachevsky-Stillinger (LS) growth ───────────────────────────────

def ls_pack(vf, n_target, seed):
    """
    Lubachevsky-Stillinger (LS) growth algorithm.

    All sphere centers are placed randomly at a tiny initial radius, then
    grown simultaneously toward the target radius. At each growth step,
    overlapping pairs are pushed apart. This overcomes the RSA jamming
    limit (~38%) and can reach 40%+ VF.

    Parameters
    ----------
    vf       : float — target volume fraction
    n_target : int   — number of spheres
    seed     : int   — random seed

    Returns
    -------
    centers  : np.ndarray, shape (N, 3)
    L        : float
    achieved : float
    """
    rng    = np.random.default_rng(seed)
    r_fin  = R_NOM
    L      = (n_target * sphere_volume(r_fin) / vf) ** (1/3)
    wall   = r_fin + 0.3

    # Place all centers randomly (radius=0 → no overlap check needed)
    centers = rng.uniform(wall, L - wall, (n_target, 3))

    # Grow radii from tiny → target in n_steps steps
    r_start = 1.0       # µm
    n_steps = 400
    radii   = np.linspace(r_start, r_fin, n_steps)

    for r in radii:
        min_d = 2 * r + MIN_GAP

        # Resolve overlaps iteratively at this radius
        for _ in range(300):
            moved = False
            for i in range(n_target):
                for j in range(i + 1, n_target):
                    diff = centers[i] - centers[j]
                    dist = np.linalg.norm(diff)
                    if dist < min_d:
                        # Push the two spheres apart along their connecting axis
                        if dist < 1e-9:
                            diff = rng.normal(0, 1, 3)
                            dist = np.linalg.norm(diff)
                        push    = (min_d - dist) / 2.0 + 0.01
                        delta   = diff / dist * push
                        centers[i] += delta
                        centers[j] -= delta
                        # Keep spheres inside the box
                        centers[i] = np.clip(centers[i], wall, L - wall)
                        centers[j] = np.clip(centers[j], wall, L - wall)
                        moved = True
            if not moved:
                break   # converged — no more overlaps at this radius

    # Final validation: remove any spheres still overlapping after growth
    min_d  = 2 * r_fin + MIN_GAP
    final  = [centers[0].tolist()]
    for i in range(1, n_target):
        arr = np.array(final)
        if not np.any(np.sum((arr - centers[i]) ** 2, axis=1) < min_d ** 2):
            final.append(centers[i].tolist())

    final    = np.array(final)
    achieved = len(final) * sphere_volume(r_fin) / L ** 3
    return final, L, achieved


# ── Geometry: build STEP via gmsh ─────────────────────────────────────────────

def build_step_raw(centers, L, vf_pct, r=None):
    """
    Use gmsh (OpenCASCADE kernel) to:
      1. Create the RVE box
      2. Create all spheres
      3. Boolean-subtract spheres from box (matrix body)
      4. Keep spheres as separate bodies
      5. Export raw STEP file

    Parameters
    ----------
    centers : np.ndarray — sphere center coordinates (µm)
    L       : float      — box side length (µm)
    vf_pct  : int        — volume fraction percentage (for filename)
    r       : float      — sphere radius (µm), defaults to R_NOM

    Returns
    -------
    tmp_path : str — path to the raw (unnamed) STEP file
    """
    if r is None:
        r = R_NOM

    tmp_path = os.path.join(OUTPUT_DIR, f"_tmp_VF{vf_pct:02d}.step")

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)   # suppress console output
    gmsh.model.add(f"AlN_Epoxy_VF{vf_pct}")

    # RVE box (will become the Epoxy matrix after subtraction)
    box = gmsh.model.occ.addBox(0, 0, 0, L, L, L)

    # AlN spheres
    sphere_tags = [
        gmsh.model.occ.addSphere(cx, cy, cz, r)
        for cx, cy, cz in centers
    ]

    # Boolean cut: box − spheres → matrix body; keep spheres (removeTool=False)
    gmsh.model.occ.cut(
        [(3, box)],
        [(3, t) for t in sphere_tags],
        removeTool=False
    )

    gmsh.model.occ.synchronize()
    gmsh.write(tmp_path)
    gmsh.finalize()

    return tmp_path


# ── Post-processing: rename PRODUCT entities in STEP ─────────────────────────

def rename_step_bodies(in_path, out_path):
    """
    Post-process the STEP file to rename PRODUCT entities.

    ANSYS DesignModeler reads the first argument of each PRODUCT() entity
    as the body name shown in the model tree. This function renames them to:
        Body 0 → "Epoxy_matrix"
        Body 1 → "AlN_particle_001"
        Body 2 → "AlN_particle_002"
        ...

    Parameters
    ----------
    in_path  : str — path to raw gmsh STEP file
    out_path : str — path to write renamed STEP file

    Returns
    -------
    n_bodies : int — total number of bodies renamed
    """
    with open(in_path) as f:
        content = f.read()

    # PRODUCT entries look like:
    # PRODUCT('Open CASCADE STEP translator 7.8 1',
    #   'Open CASCADE STEP translator 7.8 1', ...
    pattern = r"(PRODUCT\()'([^']+)'(,\s*\n\s*)'([^']+)'"
    matches  = list(re.finditer(pattern, content))

    new_content = content
    offset = 0
    for i, m in enumerate(matches):
        new_name = "Epoxy_matrix" if i == 0 else f"AlN_particle_{i:03d}"
        old = m.group(0)
        new = m.group(1) + f"'{new_name}'" + m.group(3) + f"'{new_name}'"
        pos = m.start() + offset
        new_content = new_content[:pos] + new + new_content[pos + len(old):]
        offset += len(new) - len(old)

    with open(out_path, "w") as f:
        f.write(new_content)

    return len(matches)


# ── Export: sphere coordinates CSV ───────────────────────────────────────────

def write_csv(all_results):
    """
    Write a CSV file with sphere center coordinates for all volume fractions.
    Useful for documentation, post-processing, or re-creating geometry.
    """
    csv_path = os.path.join(OUTPUT_DIR, "sphere_coordinates_all_VF.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "VF_percent", "sphere_index",
            "cx_um", "cy_um", "cz_um",
            "radius_um", "diameter_um", "RVE_side_um"
        ])
        for vf_pct, centers, L in all_results:
            for i, (cx, cy, cz) in enumerate(centers, start=1):
                writer.writerow([
                    vf_pct, i,
                    f"{cx:.4f}", f"{cy:.4f}", f"{cz:.4f}",
                    f"{R_NOM:.4f}", f"{R_NOM*2:.4f}", f"{L:.4f}"
                ])
    return csv_path


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 65)
    print("  AlN/Epoxy RVE Geometry Generator")
    print(f"  Particle diameter : {PARTICLE_DIAMETER_UM:.1f} µm")
    print(f"  Output directory  : {os.path.abspath(OUTPUT_DIR)}")
    print("=" * 65)

    all_results = []

    for vf in VOLUME_FRACTIONS:
        vf_pct   = int(round(vf * 100))
        n_target = TARGET_N[vf]
        seed     = SEEDS[vf]

        print(f"\n── VF = {vf_pct}% ──────────────────────────────────────────")

        # ── Step 1: Place spheres ──────────────────────────────────────────────
        if vf <= 0.30:
            print(f"  Algorithm : RSA  (n_target={n_target})")
            centers, L, achieved = rsa_pack(vf, n_target, R_NOM, seed)
        else:
            print(f"  Algorithm : Lubachevsky-Stillinger growth  (n_target={n_target})")
            centers, L, achieved = ls_pack(vf, n_target, seed)

        n_placed = len(centers)
        print(f"  Placed    : {n_placed} spheres")
        print(f"  Achieved VF: {achieved*100:.2f}%  (target: {vf_pct}%)")
        print(f"  RVE size  : {L:.2f} µm × {L:.2f} µm × {L:.2f} µm")

        # ── Step 2: Build geometry and export raw STEP ─────────────────────────
        print(f"  Building geometry (gmsh OpenCASCADE)...")
        tmp_path = build_step_raw(centers, L, vf_pct)

        # ── Step 3: Rename bodies in STEP for ANSYS ────────────────────────────
        out_path = os.path.join(OUTPUT_DIR, f"AlN_Epoxy_VF{vf_pct:02d}_named.step")
        n_bodies = rename_step_bodies(tmp_path, out_path)
        os.remove(tmp_path)   # clean up temp file

        sz = os.path.getsize(out_path)
        print(f"  Bodies    : {n_bodies} (1 Epoxy_matrix + {n_bodies-1} AlN_particle_*)")
        print(f"  Output    : {out_path}  ({sz // 1024} kB)")

        all_results.append((vf_pct, centers, L))

    # ── Write coordinates CSV ──────────────────────────────────────────────────
    csv_path = write_csv(all_results)
    print(f"\n  Coordinates CSV : {csv_path}")

    print("\n" + "=" * 65)
    print("  DONE")
    print("  Import any .step file into ANSYS DesignModeler:")
    print("  File → Import External Geometry File → Generate")
    print("  Then assign materials by body name in the model tree.")
    print("=" * 65)


if __name__ == "__main__":
    main()
