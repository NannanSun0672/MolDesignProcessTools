"""
Filter ligands using rules
(1)Ring counts >=3
(2)Heavy atoms <=60
(3)氢键供体和氢键受体总数>=5
(4)tPSA<=100拓扑极表面积
(5)1<=logP<=5
(6)Rotable bonds <=6
(7) QED >0.3
"""
import os
import sys
import rdkit
import rdkit.Chem as Chem
from rdkit.Chem import Lipinski,Descriptors
from rdkit.Chem import QED
from rdkit.Chem.Crippen import MolLogP
import pandas as pd

def filter_rules(idx,mol):

    RingInfos = mol.GetRingInfo()
    RingNums = RingInfos.NumRings()
    if RingNums >=2 and RingNums <=5:

       heavy_atoms = mol.GetNumHeavyAtoms()
       if heavy_atoms <=60:

            Hd = Lipinski.NumHDonors(mol)
            Ha = Lipinski.NumHAcceptors(mol)
            nums = Hd+Ha
            if nums>=5:
                if Descriptors.TPSA(mol) <=120:
                    if Descriptors.MolLogP(mol):
                        if Lipinski.NumRotatableBonds(mol)<=6:
                            logp = MolLogP(mol)
                            if logp >=1 and logp<=5:
                                return idx,Chem.MolToSmiles(mol)
                            else:
                                return None,None
                        else:
                            return None,None
                    else:
                        return None,None
                else:
                    return None,None
            else:
                return None,None
       else:
           return None,None
    else:
        return None,None

def main(ligands_dir):
    assert ligands_dir
    ligands = set()
    for filename in os.listdir(ligands_dir):
        DataPath = os.path.join(ligands_dir,filename)
        with open(DataPath,"r")as fr:
            for line in  fr.read().splitlines():
                ligands.add(line)
    print(len(ligands))
    mols = [Chem.MolFromSmiles(lig)for lig in ligands]
    smiles = set()
    for idx,mol in enumerate(mols):
        res = filter_rules(idx,mol)
        #print(mol,res)
        if res:
            smiles.add(res)
    print(len(smiles))
    with open("./Data/","w") as fw:
        for line in smiles:
            fw.write(line+"\n")
if __name__ == "__main__":
    ligands_path = ".txt"
    #main(ligands_dir)
    smiles = []
    affinity = []
    with open(ligands_path,"r")as fr:
        for line in fr.read().splitlines():
            smiles.append(line.split("\t")[0])
            affinity.append(line.split("\t")[-1])
    infos = {"smiles":[],"affinity":[]}

    for idx,lig in enumerate(smiles):
        #print(idx,lig)
        mol = Chem.MolFromSmiles(lig)
        #print(idx,mol)
        idx, ligs = filter_rules(idx, mol)
        if idx and ligs is not None:
            #print(idx,ligs)
            infos["smiles"].append(ligs)
            infos["affinity"].append(affinity[idx])
    Dataframe = pd.DataFrame(infos)
    print(Dataframe)
    Dataframe.to_csv("./Data/filtered.csv")
    saved_top1000 = "./Data/top1000.csv"
    with open(saved_top1000,"w")as fw:
        for line in Dataframe["smiles"].values[:1000]:
            fw.write(line+"\n")




