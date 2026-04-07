from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd
import argparse


def calculate_properties(smiles):
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


def main():
    parser = argparse.ArgumentParser(description="Drug-likeness filter")
    parser.add_argument("--smiles", type=str, help="Comma-separated SMILES")
    args = parser.parse_args()
    
    smiles_list = [s.strip() for s in args.smiles.split(",")]
    
    results = []
    
    for smi in smiles_list:
        props = calculate_properties(smi)
        
        if props is None:
            results.append({"SMILES": smi, "Status": "Invalid"})
            continue
        
        is_drug_like = lipinski_filter(props)
        
        props["Drug-like"] = is_drug_like
        results.append(props)
    
    df = pd.DataFrame(results)
    
    print(df)
    df.to_csv("drug_likeness_results.csv", index=False)


if __name__ == "__main__":
    main()