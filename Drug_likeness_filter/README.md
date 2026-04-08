## Drug Likeness Filter (RDKit)

A simple Python tool to evaluate the drug-likeness of molecules based on Lipinski’s Rule of Five using RDKit.
This script takes SMILES as input (directly, from a file, or manually), calculates key molecular properties, and determines whether each compound is considered drug-like (True) or not (False).

## Features 

- Accepts input in 3 different ways:
    - Direct SMILES from command line
    - CSV or Excel file
    - Manual input
- Calculates key descriptors:
    - Molecular Weight (MW)
    - LogP
    - H-Bond Donors (HBD)
    - H-Bond Acceptors (HBA)
    - Total Polar SUrface Area (TPSA)
    - Rotable bonds
- Applies Lipinski Rule of Five (Allows up to 1 violation)
- Applies Veber, PAINS and Brenk filters
- Exports results to CSV and xlsx

## Requirements

Python 3.9+ and the following libraries:
rdkit-pypi
pandas
openpyxl

## Usage

- Using SMILES directly:
    python druglikeness.py --smiles "CCO,CC(=O)O,Nc1ccccc1"
    Calculates properties for ethanol, acetic acid, and aniline.
    Saves results to properties_output.xlsx by default.
- Using a CSV or excel file:
    Example file molecules.csv
    python druglikeness.py --input_file Molecules.xlsx --output_file results
    The script will read SMILES from the file and save the results to results.xlsx.
- Without arguments
    python druglikeness.py
    The script will prompt you to enter a SMILES string manually.
    Example input: CCO

## Why this matters

The Drug-Likeness Filter matters because it quickly and objectively identifies molecules with a higher likelihood of being viable drugs, saving time and resources while guiding synthesis and assay decisions.