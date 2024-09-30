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
