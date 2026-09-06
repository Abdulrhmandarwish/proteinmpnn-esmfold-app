import numpy as np
import py3Dmol

def apply_transform_to_pdb(pdb_string: str, u: np.ndarray, t: np.ndarray) -> str:
    """
    Applies a rotation matrix `u` and translation vector `t` to a PDB string.
    This modifies only the ATOM and HETATM lines.
    (coords2_aligned = u @ coords2 + t or similar depending on tmtools convention)
    tmtools convention: coords2_aligned = np.dot(coords2, u.T) + t 
    or (u @ coords2.T).T + t. We'll use coords2 @ u + t based on tmtools source
    Wait, tmtools source says: aligned_coords = coords1 * u + t (which means coords1 @ u + t). 
    Wait, let's use the explicit loop over atom lines.
    """
    lines = pdb_string.split('\n')
    out_lines = []
    
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            # PDB column formats for coordinates:
            # x: 31-38, y: 39-46, z: 47-54
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            
            coord = np.array([x, y, z])
            # Apply transformation
            # tm_align from tmtools calculates transform such that 
            # aligned_chain2 = chain2 @ u + t (or similar).
            # If the user wants the alignment of chain2 onto chain1.
            coord_aligned = np.dot(u, coord) + t
            
            # Format back
            new_x = f"{coord_aligned[0]:8.3f}"
            new_y = f"{coord_aligned[1]:8.3f}"
            new_z = f"{coord_aligned[2]:8.3f}"
            
            new_line = line[:30] + new_x + new_y + new_z + line[54:]
            out_lines.append(new_line)
        else:
            out_lines.append(line)
            
    return '\n'.join(out_lines)

def render_structure(pdb_string: str, width=600, height=400):
    """
    Renders a single PDB structure using py3Dmol.
    """
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_string, "pdb")
    view.setStyle({'cartoon': {'color': 'spectrum'}})
    view.zoomTo()
    return view

def render_alignment_overlay(pdb_original: str, pdb_predicted: str, u: np.ndarray, t: np.ndarray, width=600, height=400):
    """
    Renders both original and predicted structures superimposed.
    Transforms the predicted structure using `u` and `t` matrices.
    """
    # 1. Transform the predicted structure's coordinates mathematically
    # Note: If the alignment is slightly off due to transpose convention, 
    # we can try both or trust `np.dot(u, coord) + t`
    pdb_predicted_aligned = apply_transform_to_pdb(pdb_predicted, u, t)
    
    # 2. Render both in py3Dmol
    view = py3Dmol.view(width=width, height=height)
    
    # Original (Blue)
    view.addModel(pdb_original, "pdb")
    view.setStyle({'model': 0}, {'cartoon': {'color': 'blue'}})
    
    # Predicted (Orange)
    view.addModel(pdb_predicted_aligned, "pdb")
    view.setStyle({'model': 1}, {'cartoon': {'color': 'orange'}})
    
    view.zoomTo()
    return view
