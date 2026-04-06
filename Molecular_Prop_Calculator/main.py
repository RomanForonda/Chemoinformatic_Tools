from rdkit import Chem
from rdkit.Chem import Descriptors

def calculate_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return None
    
    properties = {
        "Molecular Weight": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "H-Bond Donors": Descriptors.NumHDonors(mol),
        "H-Bond Acceptors": Descriptors.NumHAcceptors(mol),
        "TPSA": Descriptors.TPSA(mol),
        "Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
        "Aromatic Rings": Descriptors.NumAromaticRings(mol),
    }
    
    return properties


if __name__ == "__main__":
    smiles = input("Enter SMILES: ")
    
    props = calculate_properties(smiles)
    
    if props is None:
        print("Invalid SMILES")
    else:
        for key, value in props.items():
            print(f"{key}: {value:.2f}")