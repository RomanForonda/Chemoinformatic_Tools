from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams
import pandas as pd
import argparse
import os

# ----------------------------------------
# Funciones de cálculo de propiedades
# ----------------------------------------
def calculate_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "SMILES": smiles,
        "MW": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
        "TPSA": Descriptors.TPSA(mol),
        "RotatableBonds": Descriptors.NumRotatableBonds(mol)
    }

# ----------------------------------------
# Filtros Lipinski y Veber
# ----------------------------------------
def lipinski_filter(props):
    violations = 0
    if props["MW"] > 500: violations += 1
    if props["LogP"] > 5: violations += 1
    if props["HBD"] > 5: violations += 1
    if props["HBA"] > 10: violations += 1
    return violations <= 1

def veber_filter(props):
    return props["RotatableBonds"] <= 10 and props["TPSA"] <= 140

# ----------------------------------------
# Filtros PAINS y Brenk usando RDKit FilterCatalog
# ----------------------------------------
def create_filter_catalog():
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    return FilterCatalog.FilterCatalog(params)

def check_pains_brenk(mol, catalog):
    if mol is None:
        return {"PAINS": None, "Brenk": None}
    matches = catalog.GetMatches(mol)
    pains_flag = any(f.GetCategory() in ["PAINS A","PAINS B","PAINS C"] for f in matches)
    brenk_flag = any(f.GetCategory() == "Brenk" for f in matches)
    return {"PAINS": pains_flag, "Brenk": brenk_flag}

# ----------------------------------------
# Leer SMILES desde archivo
# ----------------------------------------
def read_smiles_from_file(file_path):
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

# ----------------------------------------
# Función principal
# ----------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Drug-likeness filter with Lipinski, Veber, PAINS, and Brenk")
    parser.add_argument("--smiles", type=str, help="Single SMILES or comma-separated list")
    parser.add_argument("--input_file", type=str, help="CSV or Excel file with a 'SMILES' column")
    parser.add_argument("--output_file", type=str, default="drug_likeness_results.csv", help="Output CSV file")
    args = parser.parse_args()

    if args.smiles:
        smiles_list = [s.strip() for s in args.smiles.split(",")]
    elif args.input_file:
        smiles_list = read_smiles_from_file(args.input_file)
    else:
        smiles = input("Enter SMILES: ")
        smiles_list = [smiles.strip()]

    catalog = create_filter_catalog()
    results = []

    for smi in smiles_list:
        props = calculate_properties(smi)
        if props is None:
            results.append({"SMILES": smi, "Status": "Invalid"})
            continue

        props["Lipinski"] = lipinski_filter(props)
        props["Veber"] = veber_filter(props)

        mol = Chem.MolFromSmiles(smi)
        flags = check_pains_brenk(mol, catalog)
        props.update(flags)

        results.append(props)

    df = pd.DataFrame(results)
    df.to_csv(args.output_file, index=False)

    print(f"\nResults saved to: {args.output_file}\n")
    print(df)

if __name__ == "__main__":
    main()