"""
Version 1.0,
written by Lu Peng, 2023.12.15
"""


import os
import csv
import glob
import shutil
import subprocess
import pandas as pd

def docking_score(i, N_for_docking, pdb_files, output_filename):
    """
    Process each peptide file and perform docking.

    Args:
    i (int): Identifier for the peptide.
    N_for_docking (int): Number of dockings to perform.
    pdb_files (list): List of pdb files for docking.
    """
    temp_folder = '_external_app'
    output_folder = f'Data_output/{output_filename}_appendix'
    input_folder = 'Data_input'

    # Ensure output_data folder exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Save current directory
    current_dir = os.getcwd()

    # Change to input_data directory for modpep
    os.chdir(temp_folder)

    for j in range(N_for_docking):
        fasta_file = os.path.join('../', temp_folder, f'peptide{i}_{j}.fasta')
        ss2_file = os.path.join('../', temp_folder, f'peptide{i}_{j}.ss2')
        models_file = os.path.join('../', output_folder, f'models{i}_{j}.pdb')

        # Generate pdb files using modpep software
        subprocess.run(['modpep', fasta_file, models_file, '-n', '1', '-L', './', '-h', ss2_file], check=True)

    # Change back to the original directory after modpep
    os.chdir(current_dir)

    # Perform molecular docking with hdock software
    for j in range(N_for_docking):
        for pdb_file in pdb_files:
            pdb_input = os.path.join(input_folder, pdb_file)
            models_file = os.path.join(output_folder, f'models{i}_{j}.pdb')
            hdock_output_file = f'Hdock{i}_{j}_{pdb_file[:-4]}.out'

            subprocess.run(['hdock', pdb_input, models_file, '-out', hdock_output_file], check=True)

            # Score the docking results using createpl software
            subprocess.run(['createpl', hdock_output_file, 'top1.pdb', '-nmax', '1', '-complex', '-models'], check=True)
            final_score_file = f'Score_{i}_{j}_{pdb_file[:-4]}.pdb'
            subprocess.run(['mv', 'model_1.pdb', final_score_file], check=True)

            # Move generated files to the output_data folder
            shutil.move(hdock_output_file, os.path.join(output_folder, hdock_output_file))
            shutil.move(final_score_file, os.path.join(output_folder, final_score_file))

def get_scores(i, output_filename):
    """
    Retrieve and calculate scores from pdb files and store them in a CSV file.

    Args:
    i (int): Identifier for the peptide.

    Returns:
    dict: A dictionary with ligand names and their scores.
    """
    output_folder = f'Data_output/{output_filename}_appendix'

    # Get pdb files from output_data folder
    pdb_files = [f for f in os.listdir(output_folder) if f.startswith(f"Score_{i}_") and f.endswith(".pdb")]
    pdb_files.sort()

    data = {}  # Dictionary to store data

    for pdb_file in pdb_files:
        with open(os.path.join(output_folder, pdb_file), "r") as f:
            lines = f.readlines()
            ligand_name, score = "", ""
            for line in lines:
                if line.startswith("REMARK Ligand:"):
                    ligand_name = line.split()[-1]
                elif line.startswith("REMARK Score:"):
                    score = line.split()[-1]
            if ligand_name and score:
                data.setdefault(ligand_name, {})[pdb_file.split("_")[-1][:-4]] = score

    # Calculate average scores
    for ligand_name, scores in data.items():
        avg_score = sum([float(score) for score in scores.values()]) / len(scores)
        scores["score"] = "{:.3f}".format(avg_score)

    # Write scores to a CSV file in output_data folder
    with open(os.path.join(output_folder, f"Get_Score{i}.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["REMARK Ligand"] + sorted(list(data.values())[0].keys()))
        for ligand_name, scores in data.items():
            writer.writerow([ligand_name] + [scores.get(pdb_file, "") for pdb_file in sorted(list(scores.keys()))])

    return data

def get_peptide_score(i, output_filename):
    """
    Get scores for peptides and sort them.

    Args:
    i (int): Identifier for the peptide.
    Returns:
    DataFrame: A DataFrame with peptides and their scores.
    """
    output_folder = f'Data_output/{output_filename}_appendix'

    # Get all CSV files in the output_data folder
    all_files = glob.glob(os.path.join(output_folder, "*.csv"))
    helix_file = os.path.join(output_folder, f"helixpool{i}.csv")
    score_file = os.path.join(output_folder, f"Get_Score{i}.csv")

    if helix_file in all_files and score_file in all_files:
        df_helix = pd.read_csv(helix_file)
        df_score = pd.read_csv(score_file)

        # Merge DataFrame objects
        df_merged = pd.concat([df_helix, df_score], axis=1)
        df_merged = df_merged.drop(columns=["REMARK Ligand"])
        df_merged.sort_values(by="score", inplace=True)

        # Save the merged DataFrame in the output_data folder
        merged_file_path = os.path.join(output_folder, f"peptide_score{i}.csv")
        df_merged.to_csv(merged_file_path, index=False)

    return df_merged

def create_dockingpool(i, N_putBack, output_filename):
    """
    Create a docking pool from sorted peptides.

    Args:
    i (int): Identifier for the peptide.
    N_putBack (int): Number of peptides to put back in the pool.

    Returns:
    list: A list of peptides in the docking pool.
    """
    output_folder = f'Data_output/{output_filename}_appendix'
    
    if N_putBack == 0:
        return []

    peptide_score_file = os.path.join(output_folder, f'peptide_score{i}.csv')
    docking_pool_file = os.path.join(output_folder, f'dockingpool{i}.csv')

    with open(peptide_score_file, 'r') as csvfile:
        dockingpool = []
        reader = csv.reader(csvfile)
        next(reader)  # Skip the header row
        for row in reader:
            dockingpool.append(row[0])
            if len(dockingpool) == N_putBack:
                break

    # Save the docking pool to a CSV file in the output_data folder
    with open(docking_pool_file, 'w') as f:
        f.write('dockingpool\n')
        for p in dockingpool:
            f.write(f"{p}\n")

    return dockingpool
