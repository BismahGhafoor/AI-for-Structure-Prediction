#!/usr/bin/env python3

#Import Packages
import os
import sys
import subprocess
import pickle
import numpy as np
import pandas as pd
from Bio import PDB
import io
import time
import traceback
import logging
import csv

# Function to run a command and capture its output
def run_command(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    return stdout.decode().strip(), stderr.decode().strip(), process.returncode

# Function to activate conda environment
def activate_conda_env(env_name):
    activate_command = f"source $(conda info --base)/etc/profile.d/conda.sh && conda activate {env_name} && python -c 'import sys; print(sys.executable)'"
    python_path, error, return_code = run_command(activate_command)
    if return_code != 0:
        print(f"Error activating environment: {error}")
        sys.exit(1)
    return python_path

# Activate the AlphaFold environment and get the Python path
python_path = activate_conda_env("alphafold_v2.3.1")

# Run the rest of the script using the activated environment's Python
if sys.executable != python_path:
    os.execl(python_path, python_path, *sys.argv)

import af2plots
print("af2plots module imported successfully!")
print(sys.executable)

class DummyJaxArray:
    def __init__(self, *args, **kwargs):
        self.shape = args[0] if args else None
        self.dtype = kwargs.get('dtype', None)

    def __getattr__(self, name):
        return lambda *args, **kwargs: self

    def __array__(self):
        import numpy as np
        if self.shape:
            return np.zeros(self.shape)
        return np.array([])

class CustomUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module == "jax._src.device_array":
            return DummyJaxArray
        return super().find_class(module, name)

def custom_load_pickle(file):
    if isinstance(file, (str, bytes, os.PathLike)):
        with open(file, 'rb') as f:
            return CustomUnpickler(f).load()
    elif isinstance(file, io.IOBase):
        return CustomUnpickler(file).load()
    else:
        raise TypeError(f"Unsupported file type: {type(file)}")

# Monkey-patch the pickle.load function in af2plots
import af2plots.plotter
af2plots.plotter.pickle.load = custom_load_pickle

# Now proceed with the rest of the script
def reverse_and_scale_matrix(matrix: np.ndarray, pae_cutoff: float = 12.0) -> np.ndarray:
    scaled_matrix = (pae_cutoff - matrix) / pae_cutoff
    scaled_matrix = np.clip(scaled_matrix, 0, 1)
    return scaled_matrix

def process_alphafold_output(base_directory: str, model_numbers: int = 5, recycling_numbers: int = 5, Protein_1 = "A", Protein_2 = "B", pae_cutoff: float = 12.0) -> pd.DataFrame:
    series_list = []

    for model_num in range(1, model_numbers+1):
        for recycling_num in range(0, recycling_numbers):
            pdb_name = f"ranked_0.pdb"
            pdb_file_path = os.path.join(base_directory, pdb_name)

            parser = PDB.PDBParser(QUIET=True)
            structure = parser.get_structure("example", pdb_file_path)

            chain_lengths = {}

            for model in structure:
                for chain in model:
                    chain_id = chain.get_id()
                    chain_length = sum(1 for _ in chain.get_residues())
                    chain_lengths[chain_id] = chain_length
                    protein_a_len = chain_lengths.get('B', 0)

            pkl_name = f"result_model_1_multimer_v3_pred_0.pkl"
            pkl_file_path = os.path.join(base_directory, pkl_name)
            with open(pkl_file_path, 'rb') as f:
                d = custom_load_pickle(f)
            iptm = d.get('iptm')
            ptm = d.get('ptm')
            pae = d.get('predicted_aligned_error')
            plddt = np.mean(d.get('plddt'))
            confidence = d.get('ranking_confidence')

            thresholded_pae = np.where(pae < pae_cutoff, 1, 0)

            local_interaction_protein_a = np.count_nonzero(thresholded_pae[:protein_a_len, :protein_a_len])
            local_interaction_protein_b = np.count_nonzero(thresholded_pae[protein_a_len:, protein_a_len:])
            local_interaction_interface_1 = np.count_nonzero(thresholded_pae[:protein_a_len, protein_a_len:])
            local_interaction_interface_2 = np.count_nonzero(thresholded_pae[protein_a_len:, :protein_a_len])
            local_interaction_interface_avg = (local_interaction_interface_1 + local_interaction_interface_2)

            scaled_pae = reverse_and_scale_matrix(pae, pae_cutoff)
            selected_values_interaction1_score = scaled_pae[:protein_a_len, protein_a_len:][thresholded_pae[:protein_a_len, protein_a_len:] == 1]
            average_selected_interaction1_score = np.mean(selected_values_interaction1_score) if selected_values_interaction1_score.size > 0 else 0
            selected_values_interaction2_score = scaled_pae[protein_a_len:, :protein_a_len][thresholded_pae[protein_a_len:, :protein_a_len] == 1]
            average_selected_interaction2_score = np.mean(selected_values_interaction2_score) if selected_values_interaction2_score.size > 0 else 0
            average_selected_interaction_total_score = (average_selected_interaction1_score + average_selected_interaction2_score) / 2

            series_list.append(pd.Series({
                'Protein_1': Protein_1,
                'Protein_2': Protein_2,
                'LIS': round(average_selected_interaction_total_score, 3),
                'LIA': local_interaction_interface_avg,
                'ipTM': round(float(iptm), 3),
                'Confidence': round(float(iptm*0.8 + ptm*0.2), 3),
                'pTM': round(float(ptm), 3),
                'pLDDT': round(plddt, 2),
                'Model': model_num,
                'Recycle': recycling_num,
                'saved folder': os.path.dirname(pdb_file_path),
                'pdb': os.path.basename(pdb_file_path),
                'pkl': os.path.basename(pkl_file_path),
            }))

    result_df = pd.concat(series_list, axis=1).T

    # Calculate Best and Average LIS and LIA
    best_lis = result_df['LIS'].max()
    best_lia = result_df['LIA'].max()
    avg_lis = result_df['LIS'].mean()
    avg_lia = result_df['LIA'].mean()

    # Apply the new conditions for positive prediction
    is_positive = ((best_lis >= 0.100) and (best_lia >= 3432)) or ((avg_lis >= 0.060) and (avg_lia >= 1610))

    return result_df, is_positive

metrics_data = {
    'Metric': ['Average LIS', 'Best LIS', 'Average LIA', 'Best LIA', 'Average ipTM', 'Best ipTM', 'Average Confidence', 'Best Confidence', 'Average pDockQ', 'Best pDockQ', 'Average pDockQ2', 'Best pDockQ2'],
    'Optimal Threshold': [0.0734, 0.203, 1610.4, 3432, 0.322, 0.38, 0.3672, 0.432, 0.133427109, 0.148516258, 0.015093248, 0.02106924],
    'Specificity': [0.926011561, 0.919075145, 0.876300578, 0.855491329, 0.937572254, 0.823121387, 0.951445087, 0.85433526, 0.804624277, 0.865895954, 0.917919075, 0.895953757],
    'Sensitivity': [0.786516854, 0.730337079, 0.767790262, 0.775280899, 0.711610487, 0.734082397, 0.674157303, 0.685393258, 0.68164794, 0.666666667, 0.692883895, 0.670411985],
    'AUC': [0.910601632, 0.890710312, 0.888699097, 0.86613626, 0.891301336, 0.862501353, 0.859056959, 0.84053387, 0.79601221, 0.818267628, 0.841622827, 0.832053863],
    "Youden's Index": [0.712528415, 0.649412223, 0.64409084, 0.630772228, 0.649182741, 0.557203784, 0.62560239, 0.539728519, 0.486272218, 0.53256262, 0.61080297, 0.566365742]
}
metrics_data_df = pd.DataFrame(metrics_data)
metrics_data_df = metrics_data_df.round(3)

base_directory = '/scratch/alice/b/bg171/FinalProject/Predictions'  # Change this to your actual working directory
model_numbers = 5
recycling_numbers = 5
Protein_1 = "A"
Protein_2 = "B"
pae_cutoff = 30
lis_threshold = 0.103
lia_threshold = 3432.0

# List all directories starting with 'Q09472_1018-1840_and_'
subdirectories = [d for d in os.listdir(base_directory) if os.path.isdir(os.path.join(base_directory, d)) and d.startswith('Q09472_320-440_566-661_1018-1840_1660-1840_1900-2100_and_')]

positive_predictions_list = []
total_predictions_list = []

for subdirectory in subdirectories:
    subdirectory_path = os.path.join(base_directory, subdirectory)
    try:
        total_prediction, is_positive = process_alphafold_output(
            subdirectory_path, model_numbers, recycling_numbers,
            Protein_1, Protein_2, pae_cutoff
        )
        if not total_prediction.empty:
            total_predictions_list.append(total_prediction)
            if is_positive:
                positive_predictions_list.append(total_prediction)
    except Exception as e:
        print(f"Error processing directory {subdirectory}: {str(e)}")

# Combine all the results into single DataFrames
if total_predictions_list:
    total_predictions_df = pd.concat(total_predictions_list, ignore_index=True)
else:
    total_predictions_df = pd.DataFrame()

if positive_predictions_list:
    positive_predictions_df = pd.concat(positive_predictions_list, ignore_index=True)
else:
    positive_predictions_df = pd.DataFrame()

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def log_message(message):
    logging.info(message)

# Write DataFrames to an Excel file with three sheets
log_message("Writing results to csv file...")

try:
    if not positive_predictions_df.empty:
        positive_predictions_df.to_csv('Positive_PPI.csv', index=False)
        log_message(f"Wrote {len(positive_predictions_df)} rows to 'Positive_PPI.csv'")
    else:
        pd.DataFrame().to_csv('Positive_PPI.csv', index=False)
        log_message("Wrote empty DataFrame to 'Positive_PPI.csv'")

    if not total_predictions_df.empty:
        total_predictions_df.to_csv('Total_Prediction.csv', index=False)
        log_message(f"Wrote {len(total_predictions_df)} rows to 'Total_Prediction.csv'")
    else:
        pd.DataFrame().to_csv('Total_Prediction.csv', index=False)
        log_message("Wrote empty DataFrame to 'Total_Prediction.csv'")

    metrics_data_df.to_csv('Optimal_Thresholds.csv', index=False)
    log_message(f"Wrote {len(metrics_data_df)} rows to 'Optimal_Thresholds.csv'")

    log_message("CSV file generation complete.")
except Exception as e:
    log_message(f"Error writing CSV files: {str(e)}")
    log_message(traceback.format_exc())

log_message("Script finished.")
