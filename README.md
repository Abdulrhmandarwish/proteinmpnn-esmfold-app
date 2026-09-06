# 🧬 ProteinMPNN + ESMFold Self-Validating Design

A Streamlit web app that demonstrates a **self-validating protein design pipeline**:

1. **Design** — Upload a backbone PDB, and [ProteinMPNN](https://github.com/dauparas/ProteinMPNN) designs a new amino acid sequence for it
2. **Refold** — [ESMFold](https://esmatlas.com/) predicts the 3D structure from the designed sequence (via API)
3. **Compare** — TM-align compares the original vs. refolded structure (TM-score, RMSD)
4. **Visualize** — Interactive 3D structural alignment + PyMOL renders

## 🚀 Live Demo

[![Open in Streamlit](https://proteinmpnn-esmfold-app-lwyo4735xvzefccepyphqf.streamlit.app/)

## 📦 Quick Start (Local)

```bash
git clone https://github.com/YOUR_USERNAME/proteinmpnn-esmfold-app.git
cd proteinmpnn-esmfold-app

# Create virtual environment
python3 -m venv venv
source venv/bin/activate

# Install dependencies (CPU-only PyTorch)
pip install -r requirements.txt

# Run the app
streamlit run app.py
```

## 🔬 How It Works

| Step | Tool | What it does |
|------|------|-------------|
| 1 | **ProteinMPNN** | Takes a protein backbone (3D coordinates) and designs a new amino acid sequence that should fold into that shape |
| 2 | **ESMFold API** | Predicts the 3D structure of the designed sequence (no local model download needed) |
| 3 | **TM-align** | Compares the original and predicted structures, returning TM-score and RMSD |

### Key Metrics

- **TM-score** ≥ 0.5 → same fold topology (success!)
- **RMSD** < 2 Å → structures are nearly identical
- If RMSD is very small, the original and predicted structures overlap almost perfectly in the alignment view

## 📁 Project Structure

```
proteinmpnn-esmfold-app/
├── app.py                          # Streamlit UI
├── requirements.txt                # Python dependencies
├── sample_6MRR.pdb                 # Sample PDB for testing
├── render_structures.py            # PyMOL rendering script
├── model_weights/
│   ├── vanilla_model_weights/      # ProteinMPNN vanilla weights
│   ├── soluble_model_weights/      # ProteinMPNN soluble weights
│   └── ca_model_weights/           # ProteinMPNN CA-only weights
└── utils/
    ├── mpnn_utils.py               # ProteinMPNN sequence design
    ├── esmfold_utils.py            # ESMFold API client
    ├── comparison_utils.py         # TM-align structural comparison
    ├── viz_utils.py                # py3Dmol 3D visualization
    └── protein_mpnn_utils.py       # Vendored ProteinMPNN core code
```

## 🛠️ Model Variants

| Variant | Description |
|---------|-------------|
| **Vanilla** | General-purpose sequence design |
| **Soluble** | Biased toward soluble protein sequences |
| **CA-only** | Uses only Cα atom positions (less structural info) |

## 📜 License

- ProteinMPNN: [MIT License](https://github.com/dauparas/ProteinMPNN/blob/main/LICENSE)
- ESMFold: Meta AI Research
- This app: MIT License
