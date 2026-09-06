from tmtools import tm_align
import numpy as np

# Simple mapping from three‑letter to one‑letter AA codes
THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    # Include unknown/stop
    "UNK": "X", "SEC": "U", "PYL": "O"
}

def _read_ca_coords_and_seq(pdb_path: str):
    """Parse a PDB file and return CA coordinates (Nx3) and the one‑letter sequence.
    Only works for standard single‑chain PDBs (the same assumption the app enforces)."""
    coords = []
    seq = []
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                # Use only CA atoms for the backbone
                if atom_name == "CA":
                    # Residue three‑letter name is columns 17‑20
                    res_name = line[17:20].strip()
                    # Convert to one‑letter code (fallback to 'X')
                    aa = THREE_TO_ONE.get(res_name, "X")
                    seq.append(aa)
                    # Coordinates: columns 30‑38, 38‑46, 46‑54
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
    if not coords:
        raise ValueError(f"No CA atoms found in {pdb_path}. Is the file malformed?")
    return np.array(coords), "".join(seq)

def compare_structures(pdb_path_1: str, pdb_path_2: str):
    """Compare two PDB structures using TM‑align.
    Returns TM‑score (normalized by chain 1), RMSD, and the alignment matrices.
    """
    coords1, seq1 = _read_ca_coords_and_seq(pdb_path_1)
    coords2, seq2 = _read_ca_coords_and_seq(pdb_path_2)
    # tm_align expects the two coordinate arrays (Nx3) and the sequences.
    result = tm_align(coords1, coords2, seq1, seq2)
    return {
        "tm_score": result.tm_norm_chain1,  # explicit requirement
        "rmsd": result.rmsd,
        "u": result.u,   # rotation matrix
        "t": result.t    # translation vector
    }
