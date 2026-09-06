# /// script
# requires-python = ">=3.10, <3.13"
# dependencies = [
#     "pymol-open-source-whl",
# ]
# ///
"""
PyMOL rendering script for ProteinMPNN + ESMFold validation pipeline.
Renders: original structure, predicted structure, and superposition overlay.

Usage:
    uv run render_structures.py <original.pdb> <predicted.pdb> <output_dir>
"""

import os
import sys

# Set environment variable for headless rendering (MUST be before pymol import)
os.environ["PYOPENGL_PLATFORM"] = "osmesa"

import pymol  # pytype: disable=import-error
pymol.pymol_argv = ["pymol", "-cq"]
pymol.finish_launching()

from pymol import cmd  # pytype: disable=import-error

def main():
    if len(sys.argv) != 4:
        print("Usage: uv run render_structures.py <original.pdb> <predicted.pdb> <output_dir>")
        sys.exit(1)

    original_pdb = sys.argv[1]
    predicted_pdb = sys.argv[2]
    output_dir = sys.argv[3]

    os.makedirs(output_dir, exist_ok=True)

    # ── Verify files exist ────────────────────────────────────────────────
    for f in [original_pdb, predicted_pdb]:
        if not os.path.exists(f):
            print(f"ERROR: File not found: {f}")
            cmd.quit()
            sys.exit(1)

    # ══════════════════════════════════════════════════════════════════════
    # 1. Render Original Structure
    # ══════════════════════════════════════════════════════════════════════
    print("=== Rendering Original Structure ===")
    cmd.load(original_pdb, "original")
    n_atoms = cmd.count_atoms("original")
    if n_atoms == 0:
        print("ERROR: No atoms loaded from original PDB")
        cmd.quit()
        sys.exit(1)
    print(f"  Loaded {n_atoms} atoms")

    cmd.show("cartoon", "original")
    cmd.color("green", "original and ss h")   # helices
    cmd.color("yellow", "original and ss s")  # sheets
    cmd.color("gray", "original and (ss l+'')")  # loops
    cmd.orient("original")
    cmd.set("ray_opaque_background", 1)
    cmd.png(os.path.join(output_dir, "original_structure.png"),
            width=1200, height=900, dpi=150)
    print(f"  Saved: {output_dir}/original_structure.png")

    cmd.delete("all")

    # ══════════════════════════════════════════════════════════════════════
    # 2. Render Predicted (ESMFold) Structure
    # ══════════════════════════════════════════════════════════════════════
    print("\n=== Rendering Predicted Structure ===")
    cmd.load(predicted_pdb, "predicted")
    n_atoms = cmd.count_atoms("predicted")
    if n_atoms == 0:
        print("ERROR: No atoms loaded from predicted PDB")
        cmd.quit()
        sys.exit(1)
    print(f"  Loaded {n_atoms} atoms")

    cmd.show("cartoon", "predicted")
    # Color by pLDDT (B-factor) — ESMFold stores pLDDT in B-factor column
    cmd.spectrum("b", "red_white_blue", "predicted")
    cmd.orient("predicted")
    cmd.set("ray_opaque_background", 1)
    cmd.png(os.path.join(output_dir, "predicted_structure_plddt.png"),
            width=1200, height=900, dpi=150)
    print(f"  Saved: {output_dir}/predicted_structure_plddt.png")

    cmd.delete("all")

    # ══════════════════════════════════════════════════════════════════════
    # 3. Superposition — Original vs Predicted
    # ══════════════════════════════════════════════════════════════════════
    print("\n=== Rendering Superposition ===")
    cmd.load(original_pdb, "original")
    cmd.load(predicted_pdb, "predicted")

    # Try align first (sequence-based), fallback to cealign
    try:
        result = cmd.align("predicted", "original")
        if result[1] < 20:
            raise ValueError("Poor sequence alignment")
        method = "align"
        rmsd = result[0]
        n_aligned = result[1]
        print(f"  Method: align (sequence-based)")
        print(f"  RMSD: {rmsd:.3f} Å over {n_aligned} atoms")
    except Exception:
        result = cmd.cealign("original", "predicted")
        method = "cealign"
        rmsd = result["RMSD"]
        n_aligned = result["alignment_length"]
        print(f"  Method: cealign (structure-based fallback)")
        print(f"  RMSD: {rmsd:.3f} Å over {n_aligned} residues")

    cmd.show("cartoon", "all")
    cmd.color("cyan", "original")
    cmd.color("salmon", "predicted")
    cmd.orient()
    cmd.set("ray_opaque_background", 1)
    cmd.png(os.path.join(output_dir, "superposition.png"),
            width=1200, height=900, dpi=150)
    print(f"  Saved: {output_dir}/superposition.png")

    # Also save session for interactive exploration
    cmd.save(os.path.join(output_dir, "session.pse"))
    print(f"  Saved: {output_dir}/session.pse")

    print(f"\n=== Done! All renders saved to {output_dir}/ ===")
    cmd.quit()


if __name__ == "__main__":
    main()
