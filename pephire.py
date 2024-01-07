"""
Version 1.0,
written by Lu Peng and Jun Ma, 2023.12.15
"""


from pephire_supply import ladderpath as lp
from pephire_supply import genPeptides as gp
from pephire_supply import psipredHelix as ph
from pephire_supply import hdockScore as hs

import os
import sys
import ast
import glob
import pandas as pd

def read_parameters(file_path):
    """
    Reads parameters from a given file and returns them as a dictionary.

    Args:
    file_path (str): The path to the parameters file.

    Returns:
    dict: Dictionary of parameters.
    """
    parameters = {}
    with open(file_path, 'r') as file:
        for line in file:
            # Ignore comments and empty lines
            if line.startswith('#') or not line.strip():
                continue
            
            # Split the line into key and value
            key, value = line.split('=', 1)
            key = key.strip()
            value = value.strip()
            
            # Attempt to convert the value to a Python object
            try:
                parameters[key] = ast.literal_eval(value)
            except (ValueError, SyntaxError):
                # Fallback to string if conversion fails
                parameters[key] = value

    return parameters

def delete_files(suffix):
    """
    Deletes files with a specific suffix in the 'temp' subdirectory.

    Args:
    suffix (str): The file suffix to delete.
    """
    temp_folder = '_external_app'
    for filename in os.listdir(temp_folder):
        if filename.endswith(suffix):
            os.remove(os.path.join(temp_folder, filename))

def select_data(threshold, output_filename):
    """
    Selects peptides with scores below a given threshold from multiple CSV files.

    Args:
    threshold (float): The threshold score for selecting peptides.

    Returns:
    list: List of selected peptides.
    """
    output_folder = f'Data_output/{output_filename}_appendix'
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_data_path = os.path.join(script_dir, output_folder)

    # List all 'peptide_score' CSV files in the output_data folder
    input_files = [os.path.join(output_data_path, f) for f in os.listdir(output_data_path) if f.startswith('peptide_score') and f.endswith('.csv')]
    output_file = os.path.join(script_dir, 'Data_output', f'{output_filename}.csv')

    # Combine data from all input files
    all_data = pd.concat([pd.read_csv(file) for file in input_files])
    selected_data = all_data[all_data['score'] < threshold]
    selected_data.to_csv(output_file, mode='a', header=True, index=False)

    return selected_data['helixpool'].tolist()

if __name__ == "__main__":
  if len(sys.argv) != 2:
    print("Usage: python pephire.py <output_filename>.csv")
    sys.exit(1)

  output_filename = sys.argv[1].replace('.csv', '')

  # if Data_output/{output_filename}_appendix does not exist, create it
  if not os.path.exists(f'Data_output/{output_filename}_appendix'):
    os.makedirs(f'Data_output/{output_filename}_appendix')

  # Path to the parameters file
  file_path = './parameters.txt'

  # Read parameters
  params = read_parameters(file_path)

  # Use the parameters in your script
  pipPool0 = params['pipPool0']
  pdb_files = params['pdb_files']
  N_newPiptide = params['N_newPiptide']
  N_for_docking = params['N_for_docking']
  N_putBack = params['N_putBack']
  N_iteration = params['N_iteration']
  threshold = params['threshold']
  limitLadderonSize = params['limitSize'] * len(pipPool0[0])

  pipPool = pipPool0

  for i in range(N_iteration):
    # Generate new peptides
    PipPoolBook = gp.getPipPoolBook(pipPool, limitLadderonSize=limitLadderonSize)
    peptides = gp.genNewPips(PipPoolBook, pipPool, N=N_newPiptide, noRepetition=True)

    # Helix prediction using PSIPRED
    ph.run_psipred(peptides, i)
    delete_files('.ss')
    delete_files('.ss2')
    ph.sort_horiz_files(i, output_filename)

    # Docking test
    helixpool = ph.create_helixpool(i, N_for_docking, output_filename)
    delete_files('.horiz')
    delete_files('.fasta')
    ph.run_psipred(helixpool, i)
    delete_files('.horiz')
    delete_files('.ss')

    # Scoring
    hs.docking_score(i, N_for_docking, pdb_files, output_filename)
    hs.get_scores(i, output_filename)
    hs.get_peptide_score(i, output_filename)

    # Update peptide pool
    dockingPool = hs.create_dockingpool(i, N_putBack, output_filename)
    unique_dockingPool = set(dockingPool) - set(pipPool)
    pipPool.extend(list(unique_dockingPool))
    print(f'Iteration {i+1} complete, current pipPool: {pipPool}')

  # Final selection of peptides based on score threshold
  select_data(threshold, output_filename)
