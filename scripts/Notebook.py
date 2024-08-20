#!/usr/bin/env python3

import os
import sys
import subprocess
import json
import pickle
import numpy as np
import pandas as pd
from Bio import PDB
from absl import flags, app, logging

# Function to activate conda environment
def activate_conda_env(env_name):
    activate_command = f"source $(conda info --base)/etc/profile.d/conda.sh && conda activate {env_name} && python -c 'import sys; print(sys.executable)'"
    process = subprocess.Popen(activate_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, executable='/bin/bash')
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        print(f"Error activating environment: {stderr.decode()}")
        sys.exit(1)
    return stdout.decode().strip()

def count_subdirectories(file_list_path):
    with open(file_list_path, 'r') as file:
        lines = file.readlines()
    
    subdirectory_count = sum(1 for line in lines if line.strip().endswith(':'))
    return subdirectory_count

# Activate the AlphaFold environment and get the Python path
print("Activating the alphafold environment")
python_path = activate_conda_env("alphafold_v2.3.1")

# Re-run the script with the activated environment's Python if needed
if sys.executable != python_path:
    os.execl(python_path, python_path, *sys.argv)

# Now proceed with the rest of the script
print("Installing the packages")
print("Setting the flags")
flags.DEFINE_string("output_dir", '/home/bg171/Project/p300/Predictions', "directory where predicted models are stored")
flags.DEFINE_string("file_list_path", '/home/bg171/Project/p300/Predictions/file_list.txt', "Path to the file_list.txt")
flags.DEFINE_string("base_path", '/home/bg171/Project/p300/Predictions', "Base path where the files on CRYOEM are located")
flags.DEFINE_string("notebook_output_dir", '/home/bg171/Project/p300/AD_Correct_Predictions_Core2', "Directory where the notebook will be saved")
flags.DEFINE_float("cutoff", 10.0, "cutoff value of PAE for IDPs. i.e. only pae<cutoff is counted good")
flags.DEFINE_boolean("create_notebook", True, "Whether creating a notebook")
flags.DEFINE_integer("surface_thres", 50, "surface threshold for IDPs. must be integer")
flags.DEFINE_integer("pae_figsize", 50, "figsize of pae_plot, default is 50")
FLAGS = flags.FLAGS

def read_file_list(file_list_path):
    with open(file_list_path, 'r') as file:
        lines = file.readlines()
    
    file_groups = {}
    current_group = None

    for line in lines:
        line = line.strip()
        if line.endswith(':'):
            current_group = line[:-1].strip()
            file_groups[current_group] = []
        elif current_group:
            file_groups[current_group].append(line.strip())
    
    return file_groups

def check_files_exist(file_groups, base_path):
    valid_groups = {}
    skipped_groups = []
    for group, files in file_groups.items():
        print(f"Checking group: {group}")
        all_files_exist = True
        missing_files = []
        for file in files:
            file_path = os.path.join(base_path, group, file)
            if not os.path.exists(file_path):
                print(f"File not found: {file_path}")
                all_files_exist = False
                missing_files.append(file_path)
            else:
                print(f"File found: {file_path}")

        if all_files_exist:
            valid_groups[group] = files
        else:
            print(f"Skipping group {group} due to missing files: {missing_files}")
            skipped_groups.append(group)
    return valid_groups, skipped_groups
def check_files_exist(file_groups, base_path):
    valid_groups = {}
    skipped_groups = []
    for group, files in file_groups.items():
        print(f"Checking group: {group}")
        all_files_exist = True
        missing_files = []
        for file in files:
            file_path = os.path.join(base_path, group, file)
            if not os.path.exists(file_path):
                print(f"File not found: {file_path}")
                all_files_exist = False
                missing_files.append(file_path)
            else:
                print(f"File found: {file_path}")

        if all_files_exist:
            valid_groups[group] = files
        else:
            print(f"Skipping group {group} due to missing files: {missing_files}")
            skipped_groups.append(group)
    return valid_groups, skipped_groups

def collect_job_data(valid_groups, base_path):
    job_data = []
    missing_ranking_debug = 0
    for group in valid_groups:
        result_subdir = os.path.join(base_path, group.strip())
        ranking_debug_path = os.path.join(result_subdir, "ranking_debug.json")
        if os.path.isfile(ranking_debug_path):
            with open(ranking_debug_path, 'r') as f:
                ranking_data = json.load(f)
            best_model = ranking_data["order"][0]
            if "iptm+ptm" in ranking_data:
                iptm_ptm_score = ranking_data["iptm+ptm"][best_model]
                job_data.append((group, iptm_ptm_score))
        else:
            print(f"ranking_debug.json not found for {group}")
            missing_ranking_debug += 1
    print(f"Number of groups missing ranking_debug.json: {missing_ranking_debug}")
    return job_data

def obtain_seq_lengths(result_subdir):
    pdb_files = [f for f in os.listdir(result_subdir) if f.endswith('.pdb')]
    seq_lengths = []
    
    parser = PDB.PDBParser(QUIET=True)
    for pdb_file in pdb_files:
        structure = parser.get_structure('structure', os.path.join(result_subdir, pdb_file))
        for model in structure:
            for chain in model:
                seq_lengths.append(len(chain))
    
    return seq_lengths

def obtain_pae_and_iptm(result_subdir, best_model):
    model_number = int(best_model.split('_')[-1]) + 1
    file_path = os.path.join(result_subdir, f"result_model_{model_number}_multimer_v3_pred_0.pkl")
    if not os.path.exists(file_path):
        print(f"Error: {file_path} does not exist")
        return None, None
    
    with open(file_path, 'rb') as f:
        data = pickle.load(f)
    
    try:
        pae_mtx = data['predicted_aligned_error']
        iptm_score = data['iptm']
        return pae_mtx, iptm_score
    except KeyError:
        print(f"Error: {file_path} does not contain the required keys")
        return None, None

def examine_inter_pae(pae_mtx, seq_lengths, cutoff):
    if pae_mtx is None:
        return False

    old_length = 0
    for length in seq_lengths:
        new_length = old_length + length
        pae_mtx[old_length:new_length, old_length:new_length] = 50
        old_length = new_length
    check = np.where(pae_mtx < cutoff)[0].size != 0

    return check

def format_path_for_notebook(path):
    return './' + os.path.basename(path)
ef create_notebook(combo, output_dir, figsize):
    print("Creating notebook")
    import nbformat as nbf

    nb = nbf.v4.new_notebook()
    output_cells = []
    md_cell = nbf.v4.new_markdown_cell(
        "# A notebook to display all the predictions with good inter-pae scores for IDPs",
    )
    import_cell = nbf.v4.new_code_cell(
        "from analysis_pipeline.af2_3dmol import parse_results, parse_results_colour_chains"
    )
    disable_autosave_cell = nbf.v4.new_code_cell(f"%autosave 0")
    output_cells.append(md_cell)
    output_cells.append(disable_autosave_cell)
    output_cells.append(import_cell)

    import_cell = nbf.v4.new_code_cell(
        "from analysis_pipeline.utils import display_pae_plots"
    )
    output_cells.append(import_cell)

    # Add the new code cell
    new_code_cell = nbf.v4.new_code_cell('''
import af2plots
print("af2plots module imported successfully!")

from analysis_pipeline.af2_3dmol import parse_results, parse_results_colour_chains
from analysis_pipeline.utils import display_pae_plots
import sys
import pickle
import io
import os

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
''')
    output_cells.append(new_code_cell)

    base_dir = output_dir
    for i in range(combo.shape[0]):
        job = combo.iloc[i, 0]
        iptm_score = combo.iloc[i, -1]
        title_cell = nbf.v4.new_markdown_cell(f"## {job} with iptm: {iptm_score}")
        output_cells.append(title_cell)
        relative_path = f'./{job.strip()}'
        subtitle1 = nbf.v4.new_markdown_cell(f"### {job} PAE plots")
        output_cells.append(subtitle1)
        code_cell_1 = nbf.v4.new_code_cell(f"display_pae_plots('{relative_path}', figsize=({figsize}, {figsize}))")
        output_cells.append(code_cell_1)
        subtitle2 = nbf.v4.new_markdown_cell(f"### {job} coloured by plddt")
        output_cells.append(subtitle2)

        code_cell_2 = nbf.v4.new_code_cell(f"parse_results('{relative_path}')")
        output_cells.append(code_cell_2)
        subtitle3 = nbf.v4.new_markdown_cell(f"### {job} coloured by chains")
        output_cells.append(subtitle3)
        code_cell_3 = nbf.v4.new_code_cell(f"parse_results_colour_chains('{relative_path}')")
        output_cells.append(code_cell_3)
    nb["cells"].extend(output_cells)
    with open(os.path.join(output_dir, "output.ipynb"), "w") as f:
        nbf.write(nb, f)
    logging.info("A notebook has been successfully created.")
    return combo.shape[0]
def main(argv):
    # Count the number of subdirectories listed in file_list.txt
    subdirectory_count = count_subdirectories(FLAGS.file_list_path)
    print(f"Number of subdirectories listed in file_list.txt: {subdirectory_count}")

    file_groups = read_file_list(FLAGS.file_list_path)
    valid_groups, skipped_groups = check_files_exist(file_groups, FLAGS.base_path)
    
    # Log the number of valid and skipped groups
    print(f"Number of valid subdirectories: {len(valid_groups)}")
    print(f"Number of skipped subdirectories: {len(skipped_groups)}")

    job_data = collect_job_data(valid_groups, FLAGS.base_path)
    
    good_jobs = []
    iptm_ptm = []
    iptm = []
    count = 0
    total_jobs = len(job_data)
    for job, iptm_ptm_score in job_data:
        logging.info(f"Now processing {job}")
        print(f"Now processing {job}")
        count += 1
        result_subdir = os.path.join(FLAGS.base_path, job.strip())
        seq_lengths = obtain_seq_lengths(result_subdir)
        best_model = "result_model_1_multimer_v3_pred_0"  # Assuming best model is 1, modify as needed
        pae_mtx, iptm_score = obtain_pae_and_iptm(result_subdir, best_model)
        check = examine_inter_pae(pae_mtx, seq_lengths, cutoff=FLAGS.cutoff)
        if check:
            good_jobs.append(job)
            iptm_ptm.append(iptm_ptm_score)
            iptm.append(iptm_score)
            print(f"Job {job} added to the good jobs list and will be included in the notebook.")
        else:
            print(f"Job {job} does not meet the criteria and will not be included in the notebook.")
        logging.info(f"Done processing {job}. {count} out of {total_jobs} finished.")

    pi_score_df = pd.DataFrame()
    pi_score_df["jobs"] = good_jobs
    pi_score_df["iptm+ptm"] = iptm_ptm
    pi_score_df["iptm"] = iptm

    pi_score_df = pi_score_df.sort_values(by="iptm", ascending=False)
    if FLAGS.create_notebook:
        num_jobs_in_notebook = create_notebook(pi_score_df, FLAGS.notebook_output_dir, FLAGS.pae_figsize)
        print(f"Number of jobs in the notebook: {num_jobs_in_notebook}")

    print(f"Number of good jobs: {len(good_jobs)}")
    print(f"Skipped groups: {skipped_groups}")

if __name__ == "__main__":
    app.run(main)

