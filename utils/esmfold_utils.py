import requests

ESMFOLD_API_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"

def predict_structure(sequence: str, device: str = "cpu") -> str:
    """
    Predicts the 3D structure of a protein sequence using the ESMFold API.
    No local model download needed — calls the free ESM Metagenomic Atlas API.
    Returns the predicted structure as a PDB-formatted string.
    """
    response = requests.post(
        ESMFOLD_API_URL,
        data=sequence,
        headers={"Content-Type": "text/plain"},
        timeout=300,
    )

    if response.status_code != 200:
        raise RuntimeError(
            f"ESMFold API returned status {response.status_code}: {response.text[:500]}"
        )

    pdb_string = response.text

    # Basic sanity check
    if "ATOM" not in pdb_string:
        raise RuntimeError("ESMFold API returned an invalid PDB (no ATOM records).")

    return pdb_string
