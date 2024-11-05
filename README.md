# TSR Package

**TSR Package** is a Python tool for retrieving Protein Data Bank (PDB) files and generating key/triplet files for protein structure analysis.

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
To retrieve PDB files using the `PDB_DL` function:

```python
from tsr_package.tsr.PDB_DL import PDB_DL

# Retrieve PDB files for the specified PDB IDs
pdb_ids = ["1GTA", "1GTB", "1lbe"]
PDB_DL(pdb_ids, 'Dataset/')
```
Or you can use a CSV file to download the PDB files:
```python
from tsr_package.tsr.PDB_DL import PDB_DL

data_dir = "Dataset"
csv_file = "sample_details.csv"
PDB_DL(csv_file, data_dir)
```

This will download the PDB files into the specified `Dataset/` directory. If the directory is not specified, the default directory for storing the PDB files would also be `Dataset/`.
Protein IDs are not case-sensitive, so you may use lowercase and uppercase letters to address the proteins.

### Generate Keys and Triplets
To generate keys or triplet files for the proteins:

```python
from tsr_package.tsr.TSR import TSR

# Define the directory where PDB files are stored
data_dir = "Dataset"
# Define the list of PDB files and corresponding chains
input_files = ["1GTA", "1GTB", "1LBE"]
chain = ["A", "A", "A"]  # specify chains for each PDB file
output_option = "keys"  # choose 'keys', 'triplets', or 'both'. If none, the function will generate both.

# Process protein data to generate key files
TSR(data_dir, input_files, chain=chain, output_option=output_option)
```
Protein's chains are case-sensitive since there are chains with both lower and uppercase letters.

### Additional TSR Functions

Besides the basic TSR function, the package provides specialized TSR functions for various analyses:

- **SSETSR**: Extracts secondary structure elements (SSEs) for each protein in a dataset. Ideal for focused structural analysis of SSEs such as helices and sheets.

- **CrossTSR**: Designed to analyze interactions between pairs of chains across different proteins. It can take single PDB files, lists, or CSV inputs with specified chain pairs to analyze cross-chain interactions.

- **DrugTSR**: Specifically built to analyze interactions between proteins and drugs. For each protein-drug pair, it generates key and triplet files, facilitating drug-binding site analysis.

Each of these functions can be customized using the same directory structure and options as the core TSR function.


### Using a CSV file as Input
You can pass a CSV file as input to process multiple PDB files with chain information. The CSV file should have the following format:

|protein         |chain        |
|----------------|-------------|
|1GTA            |A            |
|1GTB            |A            |
|1LBE            |A            |

To process the CSV file:

```python
from tsr_package.tsr.TSR import TSR

# Define the directory and CSV file path
data_dir = "Dataset"
csv_file = "sample_details.csv"

# Process the CSV input
TSR(data_dir, csv_file, output_option="keys")
```

## Arguments
- `data_dir`: Directory where the PDB files are located or where they will be downloaded.
- `input_files`: A list of PDB IDs or the path to a CSV file containing protein IDs and chains.
- `chain`: A list of chains corresponding to each PDB file (optional if using a CSV file).
- `drug_name` and `drug_id`: (For DrugTSR) Specify the drug name and ID for each protein-drug pair.
- `output_option`: Choose "keys", "triplets", or "both". Defaults to "both" if not specified.
- `aa_grouping`: Optional. Set to True if you want to use amino acid grouping labels instead of individual labels.
- `mirror_image`: Optional. Set to True if you want TSR to consider mirror image triangles.
- `size_filter`: Optional. Integer value to keep keys with `mxDist` below the specified threshold.

Additional arguments specific to each function may be found in the detailed documentation on the [TSR website].

## Examples
### Example 1: Retrieving PDB Files and Generating Keys

```python
from tsr_package.tsr.PDB_DL import PDB_DL
from tsr_package.tsr.TSR import TSR

# Step 1: Retrieve PDB files
data_dir = "Dataset" # It is also the default directory if not declared
pdb_ids = ["1GTA", "1gtb", "1lbe"] # Not case-sensitive
chain = ["A", "A", "A"] # Case-sensitive
PDB_DL(pdb_ids, data_dir)

# Step 2: Generate key files for the proteins
TSR(data_dir, pdb_ids, chain=chain, output_option="keys") # Modify the output option as desired
```

### Example 2: Using CSV File for Input

```python
from tsr_package.tsr.PDB_DL import PDB_DL
from tsr_package.tsr.TSR import TSR

# Use CSV input for batch processing
data_dir = "Dataset"
csv_file = "sample_details.csv"
PDB_DL(csv_file, data_dir)
TSR(data_dir, csv_file, output_option="triplets")
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
