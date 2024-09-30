# TSR Package

**TSR Package** is a Python tool for retrieving Protein Data Bank (PDB) files and generating key/triplet files for protein structure analysis. This package allows for batch downloading of PDB files and processing them to extract key features.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
  - [Retrieve PDB Files](#retrieve-pdb-files)
  - [Generate Keys and Triplets](#generate-keys-and-triplets)
  - [Using with CSV Input](#using-with-csv-input)
- [Arguments](#arguments)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)

## Installation

### Cloning the Repository

To get started with the `tsr_package`, clone the repository from GitHub:

```bash
git clone https://github.com/pooryakhajouie/TSR-Package.git
cd TSR-Package
```

### Installing the Package
1. It's recommended to create a virtual environment:

```bash
python3 -m venv tsrenv
source tsrenv/bin/activate  # Mac/Linux
tsrenv\Scripts\activate  # Windows
```

2. Install the package using pip:

```bash
pip install .
```

3. Alternatively, you can install the package from the built wheel:

```bash
pip install dist/tsr_package-0.1.0-py3-none-any.whl
```

4. Install the necessary dependencies:

```bash
pip install -r requirements.txt
```

## Usage
### Retrieve PDB Files
To retrieve PDB files using the `retrieve_pdb_files` function:

```python
from tsr_package.tsr.retrieve_pdb_files import retrieve_pdb_files

# Retrieve PDB files for the specified PDB IDs
pdb_ids = ['1GTA', '1GTB', '1LBE']
retrieve_pdb_files('Dataset/', pdb_ids)
```
This will download the PDB files into the specified `Dataset/` directory.
Protein IDs are not case-sensitive, so you may use lowercase and uppercase letters to address the proteins.

### Generate Keys and Triplets
To generate keys or triplet files for the proteins:

```python
from tsr_package.tsr.generate_keys_and_triplets import process_protein_data

# Define the directory where PDB files are stored
data_dir = "Dataset/"
# Define the list of PDB files and corresponding chains
input_files = ["1GTA", "1GTB", "1LBE"]
chain = ["A", "A", "A"]  # specify chains for each PDB file
output_option = "keys"  # choose 'keys', 'triplets', or 'both'. If none, the function will generate both.

# Process protein data to generate key files
process_protein_data(data_dir, input_files, chain=chain, output_option=output_option)
```

### Using a CSV file as Input
You can pass a CSV file as input to process multiple PDB files with chain information. The CSV file should have the following format:

|protein         |chain        |
|----------------|-------------|
|1GT             |A            |
|1GTB            |A            |
|1LBE            |A            |

To process the CSV file:

```python
from tsr_package.tsr.generate_keys_and_triplets import process_protein_data

# Define the directory and CSV file path
data_dir = "Dataset/"
csv_file = "sample_details.csv"

# Process the CSV input
process_protein_data(data_dir, csv_file, output_option="keys")
```

## Arguments
- `data_dir`: Directory where the PDB files are located or where they will be downloaded.
- `input_files`: A list of PDB IDs or the path to a CSV file containing protein IDs and chains.
- `chain`: A list of chains corresponding to each PDB file (optional if using a CSV file).
- `output_option`: Either "keys" to generate key files or "triplets" to generate triplet files.

## Examples
### Example 1: Retrieving PDB Files and Generating Keys

```python
from tsr_package.tsr.retrieve_pdb_files import retrieve_pdb_files
from tsr_package.tsr.generate_keys_and_triplets import process_protein_data

# Step 1: Retrieve PDB files
data_dir = "Dataset/" # It is also the default directory if not declared
pdb_ids = ["1GTA", "1gtb", "1lbe"] # Not case-sensitive
chain = ["A", "A", "A"] # Case-sensitive
retrieve_pdb_files(pdb_ids, data_dir)

# Step 2: Generate key files for the proteins
process_protein_data(data_dir, pdb_ids, chain=chain, output_option="keys") # Modify the output option as desired
```

### Example 2: Using CSV File for Input

```python
from tsr_package.tsr.generate_keys_and_triplets import process_protein_data

# Use CSV input for batch processing
data_dir = "Dataset/"
csv_file = "sample_details.csv"
process_protein_data(data_dir, csv_file, output_option="triplets")
```

## Contributing
Contributions are welcome! If you'd like to improve this package, feel free to fork the repository and submit a pull request.

1. Fork the repository
2. Create a new branch (`git checkout -b feature-branch`)
3. Commit your changes (`git commit -am 'Add new feature'`)
4. Push the branch (`git push origin feature-branch`)
5. Open a pull request

## License
This project is licensed under the MIT License - see the LICENSE file for details.
