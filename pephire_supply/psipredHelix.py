"""
Version 1.0,
written by Lu Peng, 2023.12.15
"""


import os
import subprocess
import shutil
import pandas as pd
import numpy as np

def run_psipred(peptides, i, exeName='runpsipred_single'):
    """
    Runs PSIPRED software for helix prediction on a list of peptides.

    Args:
    peptides (list): List of peptide sequences.
    i (int): Identifier for the peptide batch.
    exeName (str): Name of the executable for PSIPRED.
    """
    # Specify the temp folder path
    temp_folder = '_external_app'
    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)

    # Create .fasta files for each peptide in the temp folder
    for j, peptide in enumerate(peptides):
        fasta_filename = os.path.join(temp_folder, f'peptide{i}_{j}.fasta')
        with open(fasta_filename, 'w') as f:
            f.write(f'>peptide{j}\n{peptide}')

    # Run PSIPRED on each .fasta file in the temp folder
    for filename in os.listdir(temp_folder):
        if filename.endswith('.fasta'):
            filepath = os.path.join(temp_folder, filename)
            subprocess.check_output([exeName, filepath]).decode('utf-8')

    # Move generated files (.horiz, .ss, .ss2) to the temp folder
    for ext in ['.horiz', '.ss', '.ss2']:
        for file in os.listdir('.'):
            if file.endswith(ext):
                shutil.move(file, os.path.join(temp_folder, file))

def sort_horiz_files(i, output_filename):
    """
    Merges and sorts .horiz files by helix content percentage.
    """
    temp_folder = '_external_app'
    output_folder = f'Data_output/{output_filename}_appendix'

    # Create output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    file_names = [f for f in os.listdir(temp_folder) if f.endswith('.horiz')]
    data = []

    for file_name in file_names:
        with open(os.path.join(temp_folder, file_name)) as f:
            content = f.read()

        lines = content.splitlines()
        peptide_num = file_name.replace('.horiz', '').replace('peptide', '')
        aa, pred = lines[4][6:], lines[3][6:]
        H_count = pred.count('H')
        total_count = len(pred)
        H_percent = H_count / total_count
        data.append((int(peptide_num), aa, pred, H_percent))

    df = pd.DataFrame(data, columns=['peptide_num', 'AA', 'Pred', 'H_percent'])
    df_sorted = df.sort_values(by=['H_percent'], ascending=False)
    df_sorted.to_csv(os.path.join(output_folder, f'sorted_helix_{i}.csv'), index=False)

    return df_sorted

def create_helixpool(i, N_for_docking, output_filename):
    """
    Creates a helix pool from the sorted helix data.
    """
    if N_for_docking == 0:
        return []

    output_folder = f'Data_output/{output_filename}_appendix'
    csv_file = os.path.join(output_folder, f'sorted_helix_{i}.csv')
    df = pd.read_csv(csv_file)

    max_H_percent = df['H_percent'].max()
    max_H_percent_rows = df[df['H_percent'] == max_H_percent]
    random_rows = np.random.choice(max_H_percent_rows.index.values, size=N_for_docking, replace=False)
    helixpool = [df.loc[j, 'AA'] for j in random_rows]

    with open(os.path.join(output_folder, f'helixpool{i}.csv'), 'w') as f:
        f.write('helixpool\n')
        for p in helixpool:
            f.write(f"{p}\n")

    return helixpool
