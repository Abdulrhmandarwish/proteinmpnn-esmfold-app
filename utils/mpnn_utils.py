import os
import torch
import numpy as np
import streamlit as st

# We expect protein_mpnn_utils.py to be in the same directory (vendored)
from utils.protein_mpnn_utils import ProteinMPNN, tied_featurize, parse_PDB

@st.cache_resource
def load_mpnn_model(model_name: str, device: str = "cpu"):
    """
    Loads one of the ProteinMPNN models (vanilla, soluble, ca).
    Cached via Streamlit to avoid reloading on every interaction.
    """
    hidden_dim = 128
    num_layers = 3 
    
    # Extract just the prefix for folder name
    model_folder_map = {
        "vanilla": "vanilla_model_weights",
        "soluble": "soluble_model_weights",
        "ca": "ca_model_weights"
    }
    folder = model_folder_map.get(model_name, "vanilla_model_weights")
    
    # We use v_48_020.pt which is standard
    checkpoint_path = os.path.join("model_weights", folder, "v_48_020.pt")
    
    if not os.path.exists(checkpoint_path):
        raise FileNotFoundError(f"Model weights not found at {checkpoint_path}. Ensure you have cloned the weights.")

    checkpoint = torch.load(checkpoint_path, map_location=device)
    
    # Specific adjustment: ca_only=True for the CA model variant
    ca_only = True if model_name == "ca" else False
    
    model = ProteinMPNN(
        ca_only=ca_only, 
        num_letters=21, 
        node_features=hidden_dim, 
        edge_features=hidden_dim, 
        hidden_dim=hidden_dim, 
        num_encoder_layers=num_layers, 
        num_decoder_layers=num_layers, 
        augment_eps=0.0, 
        k_neighbors=checkpoint['num_edges']
    )
    
    model.to(device)
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()
    
    return model

def design_sequence(pdb_path: str, model_name: str = "vanilla", device: str = "cpu"):
    """
    Designs a sequence for the given PDB file using the specified model.
    """
    model = load_mpnn_model(model_name, device)
    
    ca_only = True if model_name == "ca" else False
    
    # Parse the PDB
    pdb_dict_list = parse_PDB(pdb_path, ca_only=ca_only)
    
    # Featurize
    dataset_valid = [{"name": "input", "coords_chain_A": pdb_dict_list[0]['coords_chain_A'], "seq_chain_A": pdb_dict_list[0]['seq_chain_A']}]
    
    # Note: tied_featurize logic might need adjustment if multiple chains are handled, 
    # but we enforce single chain at the app level.
    # The current parse_PDB implementation in ProteinMPNN usually puts all chains in a dict.
    # To keep it generic, we featurize all chains returned by parse_PDB.
    
    # A standard tied_featurize call for ProteinMPNN:
    batch = tied_featurize(pdb_dict_list, device, None, None, None, None, ca_only=ca_only)
    
    X, S, mask, lengths, chain_M, chain_encoding_all, chain_list_list, visible_list_list, masked_list_list, masked_chain_length_list_list, chain_M_pos, omit_AA_mask, residue_idx, dihedral_mask, tied_pos_list_of_lists_list, pssm_coef, pssm_bias, pssm_log_odds_all, bias_by_res_all, tied_beta = batch
    
    # Generate sequence
    randn_1 = torch.randn(chain_M.shape, device=device)
    
    # Specific adjustment: omit_AAs_np and bias_AAs_np must be np.zeros(21)
    omit_AAs_np = np.zeros(21)
    bias_AAs_np = np.zeros(21)
    
    with torch.no_grad():
        # model.sample needs the specific numpy arrays AND bias_by_res from featurization
        sample_dict = model.sample(X, randn_1, S, chain_M, chain_encoding_all, residue_idx, mask=mask, temperature=0.1, omit_AAs_np=omit_AAs_np, bias_AAs_np=bias_AAs_np, chain_M_pos=chain_M_pos, omit_AA_mask=omit_AA_mask, bias_by_res=bias_by_res_all)
        
    sampled_seq_idx = sample_dict["S"]
    
    alphabet = 'ACDEFGHIKLMNPQRSTVWYX'
    
    # Decode the first sample (batch size = 1)
    # the output is indices.
    seq_out = sampled_seq_idx[0].cpu().numpy()
    
    # length might be padded. mask indicates valid residues
    valid_len = int(lengths[0].item())
    seq_out = seq_out[:valid_len]
    
    designed_sequence = "".join([alphabet[idx] for idx in seq_out])
    
    # Calculate recovery
    original_seq_idx = S[0].cpu().numpy()[:valid_len]
    recovery = np.mean(seq_out == original_seq_idx) * 100.0
    
    original_seq = "".join([alphabet[idx] for idx in original_seq_idx])
    
    return designed_sequence, recovery, original_seq
