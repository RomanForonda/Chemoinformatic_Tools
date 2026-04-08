from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
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
# Filtros PAINS y Brenk
# ----------------------------------------
def create_catalogs():
    params_pains_a = FilterCatalogParams()
    params_pains_a.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    pains_a_catalog = FilterCatalog(params_pains_a)

    params_pains_b = FilterCatalogParams()
    params_pains_b.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    pains_b_catalog = FilterCatalog(params_pains_b)

    params_pains_c = FilterCatalogParams()
    params_pains_c.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
    pains_c_catalog = FilterCatalog(params_pains_c)

    params_brenk = FilterCatalogParams()
    params_brenk.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    brenk_catalog = FilterCatalog(params_brenk)

    return pains_a_catalog, pains_b_catalog, pains_c_catalog, brenk_catalog

def check_pains_brenk(mol, pains_a, pains_b, pains_c, brenk):
    if mol is None:
        return {"PAINS": None, "Brenk": None}
    pains_flag = (
        any(pains_a.GetMatches(mol)) or
        any(pains_b.GetMatches(mol)) or
        any(pains_c.GetMatches(mol))
    )
    brenk_flag = any(brenk.GetMatches(mol))
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
    parser.add_argument("--output_file", type=str, default="drug_likeness_results", help="Base name for output files")
    args = parser.parse_args()

    # Input
    if args.smiles:
        smiles_list = [s.strip() for s in args.smiles.split(",")]
    elif args.input_file:
        smiles_list = read_smiles_from_file(args.input_file)
    else:
        smiles = input("Enter SMILES: ")
        smiles_list = [smiles.strip()]

    pains_a, pains_b, pains_c, brenk = create_catalogs()
    results = []

    # Procesamiento
    for smi in smiles_list:
        props = calculate_properties(smi)
        if props is None:
            results.append({"SMILES": smi, "Status": "Invalid"})
            continue

        props["Lipinski"] = lipinski_filter(props)
        props["Veber"] = veber_filter(props)

        mol = Chem.MolFromSmiles(smi)
        flags = check_pains_brenk(mol, pains_a, pains_b, pains_c, brenk)
        props.update(flags)

        results.append(props)

    df = pd.DataFrame(results)

    # ----------------------------------------
    # Guardar en CSV + Excel
    # ----------------------------------------
    base_name = os.path.splitext(args.output_file)[0]

    csv_file = base_name + ".csv"
    excel_file = base_name + ".xlsx"

    df.to_csv(csv_file, index=False, sep=";", encoding="utf-8-sig")
    df.to_excel(excel_file, index=False)

    print(f"\nResults saved to:")
    print(f"CSV: {csv_file}")
    print(f"Excel: {excel_file}\n")

    print(df)

# ----------------------------------------
if __name__ == "__main__":
    main()