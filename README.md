# PepHiRe: Peptide Hierarchical Reconstructor

## Overview
PepHiRe (Peptide Hierarchical Reconstructor) is an innovative computational tool developed to design peptide inhibitors with high specificity and affinity for target proteins. By harnessing the power of the Ladderpath Theory, PepHiRe offers a systematic approach for the de novo generation and screening of peptide sequences that can potentially disrupt protein-protein interactions, which are pivotal in numerous biological processes and diseases.

## System Requirements and Prerequisites
- **Operating System**: Linux
- **Hardware Requirements**:
  - CPU: Quad-core processor or higher
  - RAM: 8GB or higher recommended
  - Disk Space: At least 2GB of free space for optimal performance

## Getting Started with PepHiRe
To get the most recent version of PepHiRe and stay up to date with the latest developments, we recommend cloning the project repository directly from GitHub. 

### Cloning the Repository
Open your terminal and execute the following command to clone the PepHiRe repository:
```
git clone https://github.com/yuernestliu/pephire.git
```
This command creates a local copy of the PepHiRe repository on your machine. After the cloning process is complete, navigate to the project directory:
```
cd pephire
```
### Python Environment Setup
Create a virtual environment and install the required packages:
```
# Create a virtual environment
python3 -m venv pephire_env

# Activate the virtual environment
source pephire_env/bin/activate

# Install required packages from the requirements file
pip install -r requirements.txt
```

## Setting Up Dependencies
Before using PepHiRe, it is essential to set up specific external dependencies. PepHiRe, being a Python-based program, doesn't require a separate installation process after the necessary Python packages are installed from `requirements.txt`. However, PepHiRe relies on the following software tools for computation. Make sure to install these tools in the `_external_app` directory within the pephire folder for optimal integration with PepHiRe.

### PSIPRED
PSIPRED is used for predicting the secondary structure of peptides.

```
# Download PSIPRED and move the tar.gz file to _external_app directory
wget http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/psipred.4.02.tar.gz

# Unzip and compile PSIPRED within _external_app directory
tar -xzvf psipred.4.02.tar.gz
cd psipred/src
make
make install

# Determine the PSIPRED installation directory
cd .. # back to psipred/
pwd  # Outputs the path, replace <sss> with this path.

# Update the runpsipred_single script with the correct directories
vim runpsipred_single

# In vim, press i to enter insert mode and set the following
set execdir = <sss>/bin
set datadir = <sss>/data
# After editing, press Esc, then type :wq and press Enter to save and exit vim.

# Add PSIPRED to your PATH
export PATH=$PATH:<sss>
source ~/.bashrc

# Test PSIPRED
cd example/
runpsipred_single example.fasta
```

### MODPEP
MODPEP is required for peptide modeling and can be downloaded after filling a form on the official site.(http://huanglab.phys.hust.edu.cn/modpep/)

```
# Unzip and compile package within _external_app directory:
tar -xzvf MODPEP.tar.gz

# Determine the MODPEP installation directory
cd MODPEP_v1.0
pwd  # Note down the outputted path as <sss>.

# Add MODPEP to your PATH
export PATH=$PATH:<sss>
source ~/.bashrc

# Test MODPEP:
runpsipred_single peptide.fasta
modpep peptide.fasta models.pdb -n 1 -L ./ -h peptide.ss2
```

### HDOCK
HDOCK is used for molecular docking and can be obtained after filling a form on the official site.(http://hdock.phys.hust.edu.cn/)

```
# Unzip and compile package within _external_app directory:
tar -xzvf HDOCKlite.tar.gz

# Determine the HDOCK installation directory
cd HDOCKlite-v1.1
pwd  # This will display the directory path as <sss>.

# Add HDOCK to your PATH
export PATH=$PATH:<sss>
source ~/.bashrc

# Test HDOCK:
hdock 1CGI_r_b.pdb 1CGI_l_b.pdb -out Hdock.out
createpl Hdock.out top1.pdb -nmax 1 -complex -models
cat model_1.pdb
```
Once you have successfully set up HDOCK, along with the other dependencies, you're ready to proceed to configuring PepHiRe for your projects.

## Configuration
After cloning the PepHiRe repository and installing all dependencies, you'll find several subfolders and files organized as follows:
- Subfolders:
  - `Data_input/`: Contains essential pdb files for peptide docking, including various conformations of MCL-1 protein.
  - `Data_output/`: Intended for storing the output data. Inside, you'll find `-290.csv`, `new_BH3_peptide.xlsx`, and folders named after your output file, 
  e.g., `<output_filename>.csv` and `<output_filename>_appendix`.
  - `_external_app/`: Designed to house external applications like PSIPRED, MODPEP, and HDOCK, along with specific files needed for MODPEP operations (e.g., `1kv6_C.pdb`, `helix.pdb`, etc.). Detailed installation instructions for these tools are provided earlier in the documentation.
  - `pephire_supply/`: Contains supporting scripts used by the main script `pephire.py`.
- Files:
  - `parameters.txt`: Sets various parameters for peptide generation. Customize the parameters to tailor the peptide generation process.
  - `requirements.txt`: Lists necessary Python libraries. Install these libraries using pip.
  - `pephire.py`: The main script to execute the PepHiRe tool.
  - `README.md`: This file provides detailed information about the project.
  
To configure PepHiRe, modify `parameters.txt` according to your requirements. Execute the `pephire.py` script in the terminal with the specified output file name:
```
python pephire.py <output_filename>.csv
```

## File and Folder Structure Detailed Description
- `Data_input/`: Stores pdb files for peptide docking, including five conformations of the MCL-1 protein.
- `Data_output/`:
  - `-290.csv`: Lists peptides with effective binding to MCL-1 protein, serving as a reference.
  - `new_BH3_peptide.xlsx`: Contains the initial set of 8 BH3 peptides and the final 5 selected peptides for reference.
  - `<output_filename>.csv`: Generated by running `pephire.py`, this file lists new peptides with their docking scores.
  - `<output_filename>_appendix/`: A subfolder with auxiliary results for user reference, including sorted helix files, model pdb files, HDOCK output, and score files.
- `_external_app/`: Contains auxiliary files required for the operation of PSIPRED, MODPEP, and HDOCK. This folder should include the following additional files necessary for running MODPEP and obtaining secondary structures:
`1kv6_C.pdb`, `helix.pdb`, `rotamer.pdb`, `torsionC.pdb`, `torsionN.pdb`.
Users are expected to install these external applications in this directory. Detailed installation instructions are provided in the previous sections. This directory does not include installation packages due to licensing and distribution constraints of these external tools.
- `pephire_supply/`: Contains crucial supporting scripts, such as `genPeptides.py`, `ladderpath.py`, `hdockScore.py`, and `psipredHelix.py`, which are utilized by the main script `pephire.py` for various tasks in the peptide generation process.

## Usage Instructions
To generate new peptide sequences using PepHiRe, execute `pephire.py` with the following command:
```
# Run the main script with the output filename as an argument
python pephire.py <output_filename>.csv
```
This will engage the algorithm as described by the Ladderpath Theory, iterating over the provided sequences.The output is saved as `<output_filename>.csv` in the `Data_output` folder, listing selected peptides with their docking scores.

## Additional Note
When running the `pephire.py` script using the command `python pephire.py <output_filename>.csv`, temporary files may be generated in the current directory due to the operational requirements of PSIPRED, MODPEP, and HDOCK. These temporary files are automatically moved to the `_external_app` folder, and important process files are saved in the `Data_output/<output_filename>_appendix` subfolder. Users should not be alarmed by the temporary appearance and disappearance of these files in the current directory.

## Algorithm Logic and Result Interpretation
The `pephire.py` script implements the algorithmic logic of the Ladderpath Theory to produce new peptide candidates. This method involves decomposing complex sequences into simpler recurring elements (ladderons) and recombining them based on a probabilistic model to generate new sequences.

The resulting `<output_filename>.csv` file contains peptide sequences, names of docking targets, and their average scores. Scores represent the average docking affinity across different protein conformations, aiding in the identification of effective peptide inhibitors. Should the file be empty, users are advised to execute the script multiple times or adjust the `parameters.txt` to ensure diverse generation of peptide sequences.

## License and Copyright
PepHiRe is licensed under the Apache License 2.0. This license grants you permission to use, modify, and distribute this software, but with certain conditions to protect the original work and its contributors. For detailed information, please refer to the LICENSE file found in the root directory of the project.
