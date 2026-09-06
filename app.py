import os
import subprocess
import tempfile
import torch
import numpy as np
import streamlit as st
from stmol import showmol

from utils.mpnn_utils import design_sequence
from utils.esmfold_utils import predict_structure
from utils.comparison_utils import compare_structures
from utils.viz_utils import render_structure, render_alignment_overlay

# ──────────────────────────────────────────────────────────────────────────────
st.set_page_config(page_title="Protein Design Pipeline", layout="wide")
st.title("🧬 ProteinMPNN + ESMFold Self-Validating Design")

st.sidebar.header("⚙️ Settings")
model_choice = st.sidebar.selectbox(
    "ProteinMPNN Model Variant",
    ["vanilla", "soluble", "ca"],
    format_func=lambda x: {
        "vanilla": "Vanilla (General Purpose)",
        "soluble": "Soluble (Favors Solubility)",
        "ca": "CA-only (Uses less structural info)"
    }[x]
)
device = "cuda" if torch.cuda.is_available() else "cpu"
st.sidebar.info(f"Device: **{device.upper()}**")

# ──────────────────────────────────────────────────────────────────────────────
# PyMOL render helper
# ──────────────────────────────────────────────────────────────────────────────
def run_pymol_render(orig_pdb: str, pred_pdb: str, out_dir: str) -> bool:
    script = os.path.join(os.path.dirname(__file__), "render_structures.py")
    uv = os.path.expanduser("~/.local/bin/uv")
    if not os.path.exists(uv):
        uv = "uv"
    try:
        r = subprocess.run(
            [uv, "run", script, orig_pdb, pred_pdb, out_dir],
            capture_output=True, text=True, timeout=300,
            env={**os.environ, "PYOPENGL_PLATFORM": "osmesa"}
        )
        return r.returncode == 0
    except Exception:
        return False

# ──────────────────────────────────────────────────────────────────────────────
# Session state initialization
# ──────────────────────────────────────────────────────────────────────────────
if "pdb_string" not in st.session_state:
    st.session_state.pdb_string = None
    st.session_state.input_source = None
    st.session_state.results = None  # will hold all pipeline results

# ──────────────────────────────────────────────────────────────────────────────
# Input
# ──────────────────────────────────────────────────────────────────────────────
st.markdown("---")
col_upload, col_sample = st.columns([3, 1])
with col_upload:
    uploaded_file = st.file_uploader("Upload a PDB file", type=["pdb"])
with col_sample:
    st.markdown("<br>", unsafe_allow_html=True)
    if st.button("🧪 Use Sample (6MRR)"):
        sample_path = os.path.join(os.path.dirname(__file__), "sample_6MRR.pdb")
        if os.path.exists(sample_path):
            with open(sample_path, "r") as f:
                st.session_state.pdb_string = f.read()
            st.session_state.input_source = "6MRR"
            st.session_state.results = None  # clear old results
        else:
            st.error("sample_6MRR.pdb not found.")

if uploaded_file is not None:
    new_pdb = uploaded_file.getvalue().decode("utf-8")
    if new_pdb != st.session_state.pdb_string:
        st.session_state.pdb_string = new_pdb
        st.session_state.input_source = uploaded_file.name
        st.session_state.results = None

pdb_string_input = st.session_state.pdb_string
input_source = st.session_state.input_source

if pdb_string_input is not None:
    # ── Chain Selection & Filtering ──────────────────────────────────────
    chains = sorted(list(set(
        line[21] for line in pdb_string_input.split('\n')
        if (line.startswith("ATOM") or line.startswith("HETATM")) and len(line) > 21
    )))

    if not chains:
        st.error("No valid chains found in the PDB file.")
        st.stop()

    selected_chain = chains[0]
    if len(chains) > 1:
        st.warning(f"Multiple chains detected: {chains}")
        selected_chain = st.selectbox("Select which chain to design:", chains)

    # Filter the PDB string to ONLY include the selected chain
    filtered_lines = []
    for line in pdb_string_input.split('\n'):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            if len(line) > 21 and line[21] == selected_chain:
                filtered_lines.append(line)
        else:
            # Keep headers, TER, END lines
            filtered_lines.append(line)
    
    active_pdb_string = '\n'.join(filtered_lines)

    st.success(f"✅ Loaded **{input_source}** — Designing Chain **{selected_chain}**")

    # ── Run Pipeline button ──────────────────────────────────────────────
    if st.button("🚀 Run Pipeline", type="primary"):
        tmp_orig = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
        tmp_orig.write(active_pdb_string.encode("utf-8"))
        tmp_orig.close()

        try:
            # 1. ProteinMPNN
            with st.spinner("🔧 ProteinMPNN designing sequence..."):
                designed_seq, _recovery, _orig_seq = design_sequence(
                    tmp_orig.name, model_name=model_choice, device=device
                )

            # 2. ESMFold API
            with st.spinner("🧬 ESMFold refolding via API..."):
                predicted_pdb_string = predict_structure(designed_seq, device=device)

            tmp_pred = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
            tmp_pred.write(predicted_pdb_string.encode("utf-8"))
            tmp_pred.close()

            # 3. TM-align
            with st.spinner("📐 Comparing structures..."):
                comp = compare_structures(tmp_orig.name, tmp_pred.name)

            # 4. PyMOL render
            pymol_dir = tempfile.mkdtemp(prefix="pymol_")
            pymol_images = {}
            with st.spinner("🎨 PyMOL rendering..."):
                if run_pymol_render(tmp_orig.name, tmp_pred.name, pymol_dir):
                    for name in ["superposition.png", "original_structure.png",
                                 "predicted_structure_plddt.png"]:
                        path = os.path.join(pymol_dir, name)
                        if os.path.exists(path):
                            with open(path, "rb") as f:
                                pymol_images[name] = f.read()
                    pse_path = os.path.join(pymol_dir, "session.pse")
                    if os.path.exists(pse_path):
                        with open(pse_path, "rb") as f:
                            pymol_images["session.pse"] = f.read()

            # Store everything in session_state so it survives reruns
            st.session_state.results = {
                "designed_seq": designed_seq,
                "tm_score": comp["tm_score"],
                "rmsd": comp["rmsd"],
                "u": comp["u"],
                "t": comp["t"],
                "predicted_pdb": predicted_pdb_string,
                "model": model_choice,
                "pymol_images": pymol_images,
            }

        finally:
            os.remove(tmp_orig.name)
            if 'tmp_pred' in locals():
                os.remove(tmp_pred.name)

    # ══════════════════════════════════════════════════════════════════════
    # SHOW RESULTS (from session_state — survives download clicks)
    # ══════════════════════════════════════════════════════════════════════
    if st.session_state.results is not None:
        r = st.session_state.results
        tm_score = r["tm_score"]
        rmsd = r["rmsd"]
        designed_seq = r["designed_seq"]
        predicted_pdb = r["predicted_pdb"]

        st.markdown("---")
        st.header("📊 Results")

        # ── Metrics ──────────────────────────────────────────────────────
        c1, c2, c3 = st.columns(3)
        c1.metric("🎯 TM-Score", f"{tm_score:.4f}")
        c2.metric("📏 RMSD", f"{rmsd:.2f} Å")
        c3.metric("⚙️ Model", r["model"].capitalize())

        if tm_score >= 0.5:
            st.success(f"TM-score ≥ 0.5 → **Same fold!** The designed sequence refolds correctly.")
        else:
            st.warning(f"TM-score < 0.5 → The fold may differ from the original.")

        if rmsd < 2.0:
            st.info(f"RMSD = {rmsd:.2f} Å (< 2 Å) → Structures are nearly identical.")

        # ── Designed Sequence ────────────────────────────────────────────
        st.subheader("🧬 Designed Sequence (ProteinMPNN)")
        st.code(designed_seq, language=None)
        st.caption(f"Length: {len(designed_seq)} residues")

        # ── Structural Alignment (py3Dmol) ───────────────────────────────
        st.subheader("🏗️ Structural Alignment")
        st.caption("Blue = Original backbone  •  Orange = Predicted from designed sequence")
        view_align = render_alignment_overlay(
            active_pdb_string, predicted_pdb,
            r["u"], r["t"],
            width=800, height=500
        )
        showmol(view_align, height=500, width=800)

        # ── PyMOL Renders ────────────────────────────────────────────────
        pymol_images = r.get("pymol_images", {})
        if pymol_images:
            st.subheader("🖼️ PyMOL Renders")

            if "superposition.png" in pymol_images:
                st.markdown("**Superposition** — Cyan = Original, Salmon = Predicted")
                st.image(pymol_images["superposition.png"], use_container_width=True)

            pc1, pc2 = st.columns(2)
            with pc1:
                if "original_structure.png" in pymol_images:
                    st.markdown("**Original backbone**")
                    st.image(pymol_images["original_structure.png"], use_container_width=True)
            with pc2:
                if "predicted_structure_plddt.png" in pymol_images:
                    st.markdown("**Predicted (pLDDT coloring)**")
                    st.image(pymol_images["predicted_structure_plddt.png"], use_container_width=True)

            if "session.pse" in pymol_images:
                st.download_button("📦 Download PyMOL Session (.pse)",
                                   pymol_images["session.pse"],
                                   file_name="session.pse")

        # ── Downloads ────────────────────────────────────────────────────
        st.markdown("---")
        st.subheader("💾 Download")

        fasta = f">designed_{r['model']}_{input_source}\n{designed_seq}\n"
        results_txt = (
            f"TM-score: {tm_score:.4f}\n"
            f"RMSD: {rmsd:.4f} Å\n"
            f"Model: {r['model']}\n\n"
            f"Designed sequence:\n{designed_seq}\n"
        )

        d1, d2, d3 = st.columns(3)
        with d1:
            st.download_button("📄 Designed Sequence (FASTA)", fasta,
                               file_name="designed_sequence.fasta", mime="text/plain")
        with d2:
            st.download_button("📄 Results Summary", results_txt,
                               file_name="results.txt", mime="text/plain")
        with d3:
            st.download_button("📄 Predicted PDB", predicted_pdb,
                               file_name="predicted.pdb", mime="chemical/x-pdb")
