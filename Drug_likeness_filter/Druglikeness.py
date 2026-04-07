from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd
import argparse
import os


def calculate_properties(smiles):
    """Calculate molecular properties from a SMILES string"""
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return None
    
    return {
        "SMILES": smiles,
        "MW": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol)
    }


def lipinski_filter(props):
    """Apply Lipinski Rule of 5 (allow max 1 violation)"""
    violations = 0
    
    if props["MW"] > 500:
        violations += 1
    if props["LogP"] > 5:
        violations += 1
    if props["HBD"] > 5:
        violations += 1
    if props["HBA"] > 10:
        violations += 1
    
    return violations <= 1


def read_smiles_from_file(file_path):
    """Read SMILES from CSV or Excel file"""
    ext = os.path.splitext(file_path)[1].lower()
    
    if ext == ".csv":
        df = pd.read_csv(file_path)
    elif ext == ".xlsx":
        df = pd.read_excel(file_path)
    else:
        raise ValueError("Unsupported file format. Use .csv or .xlsx")
    
    if "SMILES" not in df.columns:
        raise ValueError("File must contain a column named 'SMILES'")
    
    return df["SMILES"].tolist()


def main():
    parser = argparse.ArgumentParser(description="Drug-likeness filter using Lipinski rules")
    
    parser.add_argument("--smiles", type=str,
                        help="Single SMILES or comma-separated SMILES list")
    
    parser.add_argument("--input_file", type=str,
                        help="CSV or Excel file with a 'SMILES' column")
    
    parser.add_argument("--output_file", type=str,
                        default="drug_likeness_results.csv",
                        help="Output file name (CSV)")
    
    args = parser.parse_args()
    
    # Decide input method
    if args.smiles:
        smiles_list = [s.strip() for s in args.smiles.split(",")]
    
    elif args.input_file:
        smiles_list = read_smiles_from_file(args.input_file)
    
    else:
        smiles = input("Enter SMILES: ")
        smiles_list = [smiles.strip()]
    
    results = []
    
    for smi in smiles_list:
        props = calculate_properties(smi)
        
        if props is None:
            results.append({"SMILES": smi, "Status": "Invalid"})
            continue
        
        is_drug_like = lipinski_filter(props)
        
        props["Drug-like"] = is_drug_like
        results.append(props)
    
    # Convert to DataFrame
    df = pd.DataFrame(results)
    
    # Save results
    df.to_csv(args.output_file, index=False)
    
    print(f"\nResults saved to: {args.output_file}\n")
    print(df)


if __name__ == "__main__":
    main()