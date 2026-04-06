<<<<<<< HEAD
# Molecular Property Calculator

This project calculates key physicochemical properties of molecules from SMILES.

## Features
- Calculates molecular properties:
  - Molecular Weight
  - LogP
  - H-Bond Donors
  - H-Bond Acceptors
  - TPSA (Topological Polar Surface Area)
  - Rotatable Bonds
- Accepts SMILES input from:
  - Command line
  - CSV or Excel file with a `SMILES` column
- Generates an Excel file with the calculated properties.

## Requirements

Python 3.9+ and the following libraries:
rdkit-pypi
pandas
openpyxl

## Usage

- Using SMILES directly:
    python calc_properties.py --smiles "CCO,CC(=O)O,Nc1ccccc1"
    Calculates properties for ethanol, acetic acid, and aniline.
    Saves results to properties_output.xlsx by default.
- Using a CSV or excel file:
    Example file molecules.csv
    python calc_properties.py --input_file molecules.csv --output_file results.xlsx
    The script will read SMILES from the file and save the results to results.xlsx.
- Without arguments
    python calc_properties.py
    The script will prompt you to enter a SMILES string manually.
    Example input: CCO

## Why this matters

These properties are fundamental in drug discovery for evaluating drug-likeness and ADME behavior.


