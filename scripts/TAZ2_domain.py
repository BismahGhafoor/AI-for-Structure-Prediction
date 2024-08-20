# import the packages
import os
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import rms
import matplotlib.pyplot as plt
import seaborn as sns

# Base directory containing all subdirectories
base_dir = "/scratch/alice/b/bg171/FinalProject/Predictions"
# Path to the reference structure
reference_pdb = "/scratch/alice/b/bg171/FinalProject/fold_p300_model_0.pdb"

# Load the reference structure
ref_universe = mda.Universe(reference_pdb)
# Select only C-alpha atoms in the TAZ2 domain (residues 1724 to 1840, segment ID 'A')
reference_atoms = ref_universe.select_atoms("segid A and resid 1724:1840 and name CA")

# Initialize lists to store results
rmsd_results = []
distance_results = []
contact_maps = {}

# Function to extract common residues between reference and target
def get_common_residues(ref_atoms, target_atoms):
    ref_resids = set(ref_atoms.residues.resids)
    target_resids = set(target_atoms.residues.resids)
    common_resids = ref_resids.intersection(target_resids)
    return common_resids

# Function to select common residues and handle duplicates
def select_common_residues(universe, segid, common_resids):
    # Select atoms based on common residue IDs
    selection_str = f"segid {segid} and resid {' '.join(map(str, common_resids))} and name CA"
    all_atoms = universe.select_atoms(selection_str)
    # Remove duplicates
    unique_atoms = mda.core.groups.AtomGroup([], universe)
    seen_residues = set()
    for atom in all_atoms:
        if atom.resid not in seen_residues:
            unique_atoms += atom
            seen_residues.add(atom.resid)
    return unique_atoms

# Function to compare structures using RMSD and additional metrics
def compare_structures(pdb_file, segid, resid_range):
    print(f"Processing {pdb_file}")
    u = mda.Universe(pdb_file)
    
    # Select only C-alpha atoms in the TAZ2 domain (residues 1724 to 1840, specific segment ID)
    taz2_atoms = u.select_atoms(f"segid {segid} and resid {resid_range} and name CA")
    print(f"Selected {len(taz2_atoms)} atoms from segid {segid}, residues {resid_range}")
    
    if len(taz2_atoms) == 0:
        print(f"No atoms selected for TAZ2 domain in {pdb_file}")
        return
    
    # Get common residues
    common_resids = get_common_residues(reference_atoms, taz2_atoms)
    if not common_resids:
        print(f"No common residues found between reference and {pdb_file}")
        return
    
    # Select only common residues and remove duplicates
    ref_common_atoms = select_common_residues(ref_universe, 'A', common_resids)
    target_common_atoms = select_common_residues(u, segid, common_resids)
    
    # Validate the number of atoms in the common selection
    if len(ref_common_atoms) != len(target_common_atoms):
        print(f"Atom mismatch after selecting common residues: reference ({len(ref_common_atoms)} atoms) vs. current ({len(target_common_atoms)} atoms)")
        # Print detailed atom information for debugging
        print("Reference atoms:")
        print(ref_common_atoms.residues.resids)
        print(ref_common_atoms.names)
        print("Target atoms:")
        print(target_common_atoms.residues.resids)
        print(target_common_atoms.names)
        return

    # Calculate RMSD to the reference structure based on common residues
    rmsd_value = rms.rmsd(target_common_atoms.positions, ref_common_atoms.positions)
    print(f"RMSD to reference: {rmsd_value:.2f} Å")
    
    # Store the RMSD result
    rmsd_results.append((os.path.basename(pdb_file), rmsd_value))

    # Calculate distance between key residues (example: 1750 and 1800)
    key_residue_pairs = [(1750, 1800), (1724, 1840)]
    for resid1, resid2 in key_residue_pairs:
        atom1 = u.select_atoms(f"segid {segid} and resid {resid1} and name CA")
        atom2 = u.select_atoms(f"segid {segid} and resid {resid2} and name CA")
        if len(atom1) > 0 and len(atom2) > 0:
            distance = np.linalg.norm(atom1.positions - atom2.positions)
            distance_results.append((os.path.basename(pdb_file), resid1, resid2, distance))
            print(f"Distance between residue {resid1} and {resid2}: {distance:.2f} Å")

    # Generate contact map
    all_residues = u.select_atoms(f"segid {segid} and resid {resid_range} and name CA")
    n_residues = len(all_residues)
    contact_map = np.zeros((n_residues, n_residues))

    for i, atom1 in enumerate(all_residues):
        for j, atom2 in enumerate(all_residues):
            if i < j:
                distance = np.linalg.norm(atom1.position - atom2.position)
                contact_map[i, j] = contact_map[j, i] = distance < 5.0  # Contact threshold of 5 Å

    contact_maps[os.path.basename(pdb_file)] = contact_map

# Calculate reference distances and contact map
ref_key_residue_pairs = [(1750, 1800), (1724, 1840)]
reference_distances = []
for resid1, resid2 in ref_key_residue_pairs:
    atom1 = reference_atoms.select_atoms(f"resid {resid1} and name CA")
    atom2 = reference_atoms.select_atoms(f"resid {resid2} and name CA")
    if len(atom1) > 0 and len(atom2) > 0:
        distance = np.linalg.norm(atom1.positions - atom2.positions)
        reference_distances.append((resid1, resid2, distance))
        print(f"Reference distance between residue {resid1} and {resid2}: {distance:.2f} Å")

# Generate contact map for the reference structure
n_residues_ref = len(reference_atoms)
contact_map_ref = np.zeros((n_residues_ref, n_residues_ref))

for i, atom1 in enumerate(reference_atoms):
    for j, atom2 in enumerate(reference_atoms):
        if i < j:
            distance = np.linalg.norm(atom1.position - atom2.position)
            contact_map_ref[i, j] = contact_map_ref[j, i] = distance < 5.0  # Contact threshold of 5 Å

contact_maps["reference"] = contact_map_ref

# Loop through each subdirectory and process the PDB files
for subdir, _, files in os.walk(base_dir):
    for file in files:
        if file.startswith("ranked_0"):
            pdb_file_path = os.path.join(subdir, file)
            try:
                # Specify the residue range and segment ID for the TAZ2 domain
                taz2_resid_range = "1724:1840"  # Updated range for TAZ2 domain
                segid = "B"  # Segment ID as per the screenshot

                # Compare structures using RMSD and additional metrics
                compare_structures(pdb_file_path, segid, taz2_resid_range)
            except Exception as e:
                print(f"Failed to process {pdb_file_path}: {e}")

# Create DataFrames from the results
rmsd_df = pd.DataFrame(rmsd_results, columns=["PDB File", "RMSD"])
distance_df = pd.DataFrame(distance_results, columns=["PDB File", "Residue 1", "Residue 2", "Distance"])
reference_distance_df = pd.DataFrame(reference_distances, columns=["Residue 1", "Residue 2", "Reference Distance"])

# Display the tables
print(rmsd_df)
print(distance_df)
print(reference_distance_df)

# Plot RMSD results
# Plot RMSD results using swarm plot
plt.figure(figsize=(4, 8))
sns.swarmplot(y=rmsd_df["RMSD"], color='skyblue')
plt.ylabel('RMSD (Å)')
plt.title('RMSD Distribution of TAZ2 Domain to Reference Structure')
plt.tight_layout()
plt.savefig("/scratch/alice/b/bg171/FinalProject/Predictions/rmsd_swarmplot.png")
plt.show()

# Plot Distance results
for (resid1, resid2), group in distance_df.groupby(['Residue 1', 'Residue 2']):
    plt.figure(figsize=(10, 6))
    plt.scatter(group.index, group["Distance"], color='skyblue')
    plt.axhline(y=reference_distance_df[(reference_distance_df['Residue 1'] == resid1) & (reference_distance_df['Residue 2'] == resid2)]["Reference Distance"].values[0], color='r', linestyle='--')
    plt.xlabel('Model Index')
    plt.ylabel('Distance (Å)')
    plt.title(f'Distance between residues {resid1} and {resid2} in TAZ2 Domain')
    plt.tight_layout()
    plt.savefig(f"/scratch/alice/b/bg171/FinalProject/Predictions/distance_{resid1}_{resid2}_scatter.png")
    plt.show()

# Plot contact maps
for pdb_file, contact_map in contact_maps.items():
    plt.figure(figsize=(10, 8))
    sns.heatmap(contact_map, cmap='Blues')
    plt.title(f'Contact Map of TAZ2 Domain - {pdb_file}')
    plt.xlabel('Residue Index')
    plt.ylabel('Residue Index')
    plt.tight_layout()
    plt.savefig(f"/scratch/alice/b/bg171/FinalProject/Predictions/contact_map_{pdb_file}.png")
    plt.show()

# Plot difference contact map between reference and ranked_0
contact_map_diff = contact_maps["ranked_0.pdb"] - contact_maps["reference"]
plt.figure(figsize=(10, 8))
sns.heatmap(contact_map_diff, cmap='coolwarm', center=0)
plt.title('Difference Contact Map (Bound - Reference)')
plt.xlabel('Residue Index')
plt.ylabel('Residue Index')
plt.tight_layout()
plt.savefig("/scratch/alice/b/bg171/FinalProject/Predictions/contact_map_diff_ranked_0.png")
plt.show()

# Save DataFrames to CSV files
rmsd_df.to_csv("/scratch/alice/b/bg171/FinalProject/Predictions/rmsd_results.csv", index=False)
distance_df.to_csv("/scratch/alice/b/bg171/FinalProject/Predictions/distance_results.csv", index=False)
reference_distance_df.to_csv("/scratch/alice/b/bg171/FinalProject/Predictions/reference_distance_results.csv", index=False)
