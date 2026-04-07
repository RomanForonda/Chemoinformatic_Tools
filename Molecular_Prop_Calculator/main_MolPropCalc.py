from rdkit import Chem      #Importing the Chem module from RDKit to handle chemical structures
from rdkit.Chem import Draw      #Importing the Draw module to visualize chemical structures
from rdkit.Chem import Descriptors      #Importing the Descriptors module to calculate molecular properties
import pandas as pd     #Importing pandas for data manipulation and analysis
import argparse    #Importing argparse to handle command-line arguments
import os   #Importing os to handle file paths and operations

def calculate_properties(smiles_list):      
    """
    Calculates basic molecular properties and structure representations for a list of SMILES.
    Returns a DataFrame with the results.
    """
    results = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)    #Converting the SMILES string to a molecule object using RDKit
        if mol is None:
            props = {"SMILES": smiles, "Error": "Invalid SMILES"}
        else:
            props = {
                "SMILES": smiles,
                "Molecular Weight": Descriptors.MolWt(mol),
                "LogP": Descriptors.MolLogP(mol),
                "H-Bond Donors": Descriptors.NumHDonors(mol),
                "H-Bond Acceptors": Descriptors.NumHAcceptors(mol),
                "TPSA": Descriptors.TPSA(mol),
                "Rotatable Bonds": Descriptors.NumRotatableBonds(mol)
            }
        results.append(props)
    return pd.DataFrame(results)

def read_smiles_from_file(file_path):
    """
    Reads SMILES from an Excel or CSV file.
    Expects a column named 'SMILES'.
    """
    ext = os.path.splitext(file_path)[1].lower()
    if ext == ".xlsx":
        df = pd.read_excel(file_path)
    elif ext == ".csv":
        df = pd.read_csv(file_path)
    else:
        raise ValueError("Unsupported file format. Use .xlsx or .csv")
    
    if 'SMILES' not in df.columns:
        raise ValueError("File must contain a column named 'SMILES'")
    
    return df['SMILES'].tolist()

def main():     #Main function to handle command-line arguments and execute the property calculation
    parser = argparse.ArgumentParser(description="Calculate molecular properties from SMILES")
    parser.add_argument("--smiles", type=str, help="Single SMILES or comma-separated list of SMILES")
    parser.add_argument("--input_file", type=str, help="Excel/CSV file with 'SMILES' column")
    parser.add_argument("--output_file", type=str, default="properties_output.xlsx",
                        help="Output Excel file name")
    
    args = parser.parse_args()
    
    smiles_list = []
    
    if args.smiles:
        smiles_list = [s.strip() for s in args.smiles.split(",")]
    elif args.input_file:
        smiles_list = read_smiles_from_file(args.input_file)
    else:
        smiles = input("Enter SMILES: ")
        smiles_list = [smiles.strip()]
    
    df_props = calculate_properties(smiles_list)
    df_props.to_excel(args.output_file, index=False)
    print(f"Properties calculated and saved to '{args.output_file}'")
    print(df_props)

if __name__ == "__main__":
    main()