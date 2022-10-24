"""
Create on 28 Feb,2021

@author:nannan.sun@wowjoy.cn

Function:
分子结构相似性工具

1.比较MOSES metric 的 utiles里的代码 和Tanimoto coefficient between the chemical fingerprints

"""

import rdkit
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import QED

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import xlsxwriter
import os
import sys
import shutil
import time
loca = time.strftime('%m%d%H%M')



class LigandSimility(object):

    def __init__(self, output_path, input_path, img_saved_dir):
        self.book_path = output_path
        self.input_path = input_path
        self.img_saved_dir = img_saved_dir
        ligands = self.read_ligands(self.input_path)
        self.mols = []
        for idx, smi in enumerate(ligands):
            m = Chem.MolFromSmiles(smi)
            self.mols.append(m)

        self.sub_Info_dict = {}
        self.sub_Info_dict['ligand'] = ligands

    def read_ligands(self, input_path):
        ligands = []
        if input_path.endswith('.xlsx') or input_path.endswith('.xls'):
            import xlrd
            book = xlrd.open_workbook(input_path)
            sheet1 = book.sheets()[0]
            nrows = sheet1.nrows
            for i in range(1, nrows):
                line = sheet1.row_values(i)
                smiles = line[1].strip() # 选择smiles所在的列
                if len(smiles) < 2:
                    continue
                ligands.append(smiles)
        else:
            with open(input_path,"r") as fr:
                content = fr.readlines()
            ligands = [line.strip() for line in content]
        # print('input ligands:', ligands)

        return ligands

    def fingerprints_based(self):
        """
        采用fingerprint 和Tanimoto 方法计算分子相似性
        """
        ### Topological path-based fingerprint
        topo_fps = [Chem.RDKFingerprint(x) for x in self.mols]
        # print('topo_fps:', len(topo_fps))
        ### Morgan
        morgan_fps = [AllChem.GetMorganFingerprintAsBitVect(x, 2, 2048) for x in self.mols]
        # print('morgan_fps:', len(morgan_fps))
        morgan = DataStructs.BulkTanimotoSimilarity(morgan_fps[0], morgan_fps[:])
        # print('morgan:', morgan)
        smi_res = {"similarty_index": [], "smi_score": []}
        similarity_matrix = np.zeros((len(topo_fps), len(topo_fps)), np.float)
        # print('similarity_matrix ', similarity_matrix)
        sim_index_list = list()
        for i in range(len(morgan_fps)):
            Morgan_based_similarity = list()
            Topological_based_similarity = list()
            for idx, smi in enumerate(morgan_fps):
                morgan_smilarity = DataStructs.FingerprintSimilarity(morgan_fps[i], smi)
                Morgan_based_similarity.append(morgan_smilarity)

            for idx, smi in enumerate(topo_fps):
                topo_smilarity = DataStructs.FingerprintSimilarity(topo_fps[i], smi)
                Topological_based_similarity.append(topo_smilarity)

            for j, simi_score in enumerate(Morgan_based_similarity):
                if simi_score < 0.8 and Topological_based_similarity[j] < 0.8:
                    Reset_similarity = 0.6*Topological_based_similarity[j] + 0.4*simi_score
                    similarity_matrix[i,j] = Reset_similarity
                elif simi_score >= Topological_based_similarity[j]:
                    similarity_matrix[i,j] = simi_score
                elif simi_score < Topological_based_similarity[j]:
                    similarity_matrix[i,j] = Topological_based_similarity[j]
        
        ####存储数据
        similarity_index = list()
        for i in range(similarity_matrix.shape[0]):
            sim_index = ""
            for j in range(similarity_matrix.shape[1]):
                if similarity_matrix[i,j] > 0.5 and i < j:
                    index = j+1
                    sim_index += str(index)+","
            if len(sim_index) == 0:
                similarity_index.append("None")
            else:
                similarity_index.append(sim_index)
        
        self.sub_Info_dict["similarity_index"] = similarity_index
        qeds = [QED.qed(mol) for mol in self.mols]
        self.sub_Info_dict["qed"] = qeds
        
        columns = list(self.sub_Info_dict.keys())
        # print('scores ', self.sub_Info_dict["scores"])
        self.Insert_Pictures(self.sub_Info_dict, "ligands", columns)

        #self.Plot_HeatMap(similarity_matrix)

    def num_convert(self,x):
        aff_value = -np.log10(float(x) * 1e-9)
        return aff_value

    def saved_results(self):
        """
        存储结果
        """
        pass

    def Plot_HeatMap(self,similarity_matrix):
        """
        绘制HeatMap图
        """
        print(similarity_matrix.shape)
        x_ticks = range(similarity_matrix.shape[0])
        y_ticks = range(similarity_matrix.shape[0])  # 自定义横纵轴
        ax = sns.heatmap(similarity_matrix, xticklabels=x_ticks, yticklabels=y_ticks,cmap="Reds",linecolor="black")
        ax.set_title('Similarity matrix')  # 图标题
        ax.set_xlabel('x label')  # x轴标题
        ax.set_ylabel('y label')
        plt.show()
        figure = ax.get_figure()
        figure.savefig('{}/16_similarity_matrix.jpg'.format(self.img_saved_dir))  # 保存图片

    def Insert_Pictures(self, info_dict, sheet_name, column_name):
        book = xlsxwriter.Workbook(self.book_path)
        column_name.append("Images")
        sheet = book.add_worksheet(sheet_name)
        
        img_saved_path = os.path.join(self.img_saved_dir, sheet_name)
        if not os.path.exists(img_saved_path):
            os.makedirs(img_saved_path)
        ### set row's heigth
        for i in range(0, 1000):
            # if i > 0:
            sheet.set_row(i, 230)
        sheet.set_column("A:A", 100)
        sheet.set_column("E:E", 100)
        if "scores" in column_name:
            column_name.remove("scores")
        
        # for idx, name in enumerate(column_name):
        #     sheet.write(0, idx, name)
        for idx, ligand in enumerate(info_dict["ligand"]):
            if idx >= 0:
                m = Chem.MolFromSmiles(ligand)
                img = Draw.MolToImageFile(m, img_saved_path + "/mol" + str(idx) + ".png")
                sheet.write(idx, 0, ligand)
                sheet.write(idx, 1, info_dict["similarity_index"][idx])
                sheet.write(idx, 2, info_dict["qed"][idx])
                sheet.insert_image(idx, 3, img_saved_path + "/mol" + str(idx) + ".png")
        book.close()
        shutil.rmtree(img_saved_path)


if __name__ == "__main__":
    input_path = sys.argv[1] # 输入smiles序列，txt或csv文件格式
    img_saved_dir = "./smilarity" # 图片保存路径
    book_path = img_saved_dir + '/smilarity_{}.xlsx'.format(loca) # xlsx结果保存文件
    ligand_Simility = LigandSimility(book_path, input_path, img_saved_dir)
    ligand_Simility.fingerprints_based()
    print('Done!')
