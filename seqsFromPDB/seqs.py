'''
Description: 
Author: Jiangjing
Date: 2021-05-10 02:50:34
LastEditTime: 2022-08-07 16:04:24
LastEditors: luisa_jiangj luisajiangj@outlook.com
'''

import os
import sys
import time
import xlrd
loca = time.strftime('%m%d%H%M')
from utils import *
from tqdm import tqdm


def postprocess_for_MONN():
  '''
  对MONN预测出的多个结果做整合排序
  '''
  result_file = "./MONN_" + 'CXXC-1-long-1' + '_rank.txt' # 排序文件
  results = {}
  intermediate = './test/filter_ligands/MONN_05100910'
  for file in os.listdir(intermediate):
    if file.endswith('.txt'):
      filepath = os.path.join(intermediate, file)
      with open(filepath, 'r') as f:
        for row in f.readlines():
          row = row.strip().split()
          if row[0] not in results.keys():
            results[row[0]] = row[1]
  sorted_res = sorted(results.items(), key=lambda kv:(kv[1], kv[0]), reverse=True)
  with open(result_file, 'w+') as f:
    for r in sorted_res[:1000]:
      f.write(r[0] + '\t' + r[-1] + '\n')


def smiles_to_graph():
  '''
  根据输入列表，将smiles转换为结构图 
  '''
  from rdkit import Chem
  from rdkit.Chem import Draw
  file = './ttttt.txt'
  with open(file, 'r') as fs:
      rows = fs.readlines()
      img_save_dir = './structures'
      if not os.path.exists(img_save_dir):
        os.makedirs(img_save_dir)
      for i, row in enumerate(rows):
        smiles = row.strip().split()[0]
        mol = Chem.MolFromSmiles(smiles)
        impath = os.path.join(img_save_dir, str(i) + '.png')
        img = Draw.MolToImageFile(mol, impath) # 画出分子结构


def get_protein_seq():
  file = './example/6x1a/6x1a_seq.txt'
  fasta = './example/6x1a/6x1a.fasta'
  with open(file, 'r') as f1, open(fasta, 'w') as f2:
    line = f1.readline().strip().split()
    print(line)
    new_line = ''.join(line)
    f2.write(new_line + '\n')
    print(new_line)


def get_protein_seq_v1():
  ''' 
  获取蛋白序列 
  '''
  protein_path = "./project_yupeilin/user_103f1722-4bb9-43d8-84d7-07566cdec379/user_103f1722-4bb9-43d8-84d7-07566cdec379_pocket1.pdb"
  mapping_codes = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E",
    "PHE": "F", "GLY": "G", "HIS": "H", "LYS": "K",
    "ILE": "I", "LEU": "L", "MET": "M", "ASN": "N",
    "PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S",
    "THR": "T", "VAL": "V", "TYR": "Y", "TRP": "W"
  }

  seq = ''
  amino = []
  with open(protein_path, 'r') as f:
    rows = f.readlines()
    for row in rows:
      row = row.split()
      if row[0] == 'ATOM':
        if row[3] in mapping_codes.keys():
          pair = [mapping_codes[row[3]], row[5]]
          if pair not in amino:
            amino.append(pair)
          # seq.append(mapping_codes[row[3]])
        else:
          # seq.append('X')
          amino.append(['X', row[5]])
  # seq = ''.join(seq)
  for ele in amino:
    seq = seq + ele[0]
  print(seq)


def cal_ligand_similarity():
  ''' 
  计算分子拓扑相似性 
  '''
  smi1 = 'O=C1CCC2CN(C(=O)C3CCc4[nH]ncc4C3)CCC2N1'
  smi2 = 'CN(CCC1CCCCN1C)C(=O)C1CCc2[nH]ncc2C1'
  from rdkit import Chem, DataStructs
  mol1 = Chem.MolFromSmiles(smi1)
  mol2 = Chem.MolFromSmiles(smi2)
  fp1 = Chem.RDKFingerprint(mol1)
  fp2 = Chem.RDKFingerprint(mol2)
  sim = DataStructs.FingerprintSimilarity(fp1, fp2)
  print(sim)


def convert():
  smi_file = './example/p2x4/ligands.txt'
  with open(smi_file, 'r') as f:
    smiles = [smi.strip() for smi in f.readlines()]
  
  save_dir = 'convert'
  if not os.path.isdir(save_dir):
    os.makedirs(save_dir)
    
  for i, smi in enumerate(smiles):
    ligand_smi_path = os.path.join(save_dir, str(i+1) + '.smi')
    cmd = "obabel -:\"%s\" -ocan > %s" % (smi, ligand_smi_path)
    t = os.system(cmd)

    # transform .smi file to 3D .sdf file
    ligand_sdf_path = os.path.join(save_dir, str(i+1) + '.sdf')
    cmd = "obabel %s -O %s --gen3D --minimize --steps 200 --sd --ff MMFF94" % \
          (ligand_smi_path, ligand_sdf_path)
    t = os.system(cmd)

    ligand_pdb_path = ligand_sdf_path.replace('sdf', 'pdb')
    cmd = "obabel %s -O %s" % (ligand_sdf_path, ligand_pdb_path)
    t = os.system(cmd)

    os.remove(ligand_smi_path)
    os.remove(ligand_sdf_path)


def move_ligand():
  ''' 
  读入分子smiles并生成移动到对应受体位子的小分子sdf文件 
  '''
  ligand_base_dir = './move_ligands'
  if not os.path.isdir(ligand_base_dir):
    os.makedirs(ligand_base_dir)

  file = './example/yy_TRPC6/ligands1_1.txt'
  smiles_list = []
  with open(file, 'r') as f:
    for line in f.readlines():
      line = line.strip()
      if not line:
        continue
      if line not in smiles_list:
        smiles_list.append(line)
  
  pocket_path = './project_trpc6/found_pockets/Pocket_hTRPC6.pdb'
  pocket_center = cal_pocket_mean(pocket_path)
  for i, smiles in enumerate(smiles_list):
    print('smiles {}:{}'.format(i+1, smiles))
    ### using rdkit to generate 3D data
    try:
      m1 = Chem.MolFromSmiles(smiles)
      m2 = Chem.AddHs(m1) # 加氢原子
      AllChem.EmbedMolecule(m2) # 2D->3D化
      AllChem.MMFFOptimizeMolecule(m2)   # 使用MMFF94最小化RDKit生成的构象
      m3 = Chem.RemoveHs(m2) # 删除氢原子
      ligand_sdf_path = os.path.join(ligand_base_dir, str(i+1) + '.sdf')
      w = Chem.SDWriter(ligand_sdf_path)
      w.write(m3)
    except:
      continue
    
    ligand_sdf_path2 = os.path.join(ligand_base_dir, str(i+1) + '_move.sdf')
    ret = moving_ligand(ligand_sdf_path, ligand_sdf_path2, pocket_center)
    if ret == 0:
      continue
    os.remove(ligand_sdf_path)


def convert_sdf2smi():
  file = './example/zy/1_pha.sdf'
  mols = [mol for mol in Chem.SDMolSupplier(file)]
  outname = file.replace('.sdf', '.txt')
  with open(outname, 'w') as f:
    for mol in mols:
      smi = Chem.MolToSmiles(mol)
      f.write(smi + '\n')


def filter_molecules(f1, f2):
  ''' 
  过滤大的分子库 
  '''
  from rdkit import Chem
  from rdkit.Chem import QED
  from rdkit.Chem import Lipinski
  from rdkit.Chem.Crippen import MolLogP
  from rdkit.Chem.Descriptors import ExactMolWt
  from rdkit.Chem.rdMolDescriptors import CalcTPSA

  with open(f1, 'r') as f1, open(f2, 'w') as f2:
    for line in tqdm(f1.readlines()):
      try:
        smi = line.strip().split()[0]
        mol = Chem.MolFromSmiles(smi)
        hba = Lipinski.NumHAcceptors(mol)
        hbd = Lipinski.NumHDonors(mol)
        logp = MolLogP(mol)
        rotbonds = Lipinski.NumRotatableBonds(mol)
        rings = mol.GetRingInfo().NumRings()
        hvyatoms = mol.GetNumHeavyAtoms()
        tpsa = CalcTPSA(mol)
        if (hba + hbd) >= 5 and (logp >= 1 and logp <= 5) and (rings >= 2 and rings <=5) \
          and hvyatoms <= 60 and rotbonds <= 6 and tpsa <= 120:
          f2.write(smi + '\n')
      except:
        continue


if __name__ == '__main__':
  if sys.argv[1] == 'p':
    postprocess_for_MONN()
  elif sys.argv[1] == 'g':
    smiles_to_graph()
  elif sys.argv[1] == 's':
    get_protein_seq_v1()
  elif sys.argv[1] == 's1':
    file = sys.argv[2]
    get_protein_seq_v1(file)
  elif sys.argv[1] == 'ls':
    cal_ligand_similarity()
  elif sys.argv[1] == 'c':
    convert()
  elif sys.argv[1] == 'm':
    move_ligand()
  elif sys.argv[1] == 'csdf':
    convert_sdf2smi()
  elif sys.argv[1] == 'f':
    filter_molecules(sys.argv[2], sys.argv[3])
      
