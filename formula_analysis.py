#计算机模型方剂中的扰动信号，殇和Sab分数
#误删了求方剂对应平均靶点的函数，待补充


import pandas as pd
import Data_input as di
import os
import herb_pairs_from_formula as hpff
import Data_output as do
import PPI_analyse as ppi
import numpy as np
import networkx as nx

#计算方剂针对特定疾病(靶点)的横向和纵向扰动分数。
def perturbation_formula(herbs_list,disease_target_symbol_list,herbs_data):
    herbs_score_dict = {}

    herbs_data = herbs_data[herbs_data['Herb_Chinese_Name'].isin(herbs_list)]
    herbs_data = herbs_data[disease_target_symbol_list]

    #按行求和，添加为新列
    herbs_data['row_sum'] = herbs_data.apply(lambda x:x.sum(),axis=1)
    #按列求和，添加为新行
    herbs_data.loc['col_sum'] = herbs_data.apply(lambda x:x.sum())
    herbs_data_row = list(herbs_data['row_sum'])
    herbs_data_col = list(herbs_data.loc['col_sum'])

    sum_row = 0
    sum_col = 0
    for i in range(len(herbs_data_row)-1):
        sum_row = sum_row + np.abs(herbs_data_row[i])
    for j in range(len(herbs_data_col)-1):
        sum_col = sum_col + np.abs(herbs_data_col[j])

    return sum_row,sum_col

def fomula_sab(herbs_list,herb_pair_from_data):

    if len(herbs_list) == 1:
        return 0
    fomula_sab_sum = 0
    for i in range(0,len(herbs_list)-1):
        for j in range(i+1,len(herbs_list)):
            herb1 = herbs_list[i]
            herb2 = herbs_list[j]
            herb1herb2_sum = 0
            if str(herb1)+str(herb2) in herb_pair_from_data:
                herb1herb2_sum = herb_pair_from_data[str(herb1)+str(herb2)]
            if str(herb2)+str(herb1) in herb_pair_from_data:
                herb1herb2_sum = herb_pair_from_data[str(herb2)+str(herb1)]
            fomula_sab_sum = fomula_sab_sum + herb1herb2_sum
    fomula_sab_sum = float(fomula_sab_sum)/len(herbs_list)
    return fomula_sab_sum


def formula_entropy(herbs_list,herb_pair_from_data_shangshi,bins):#计算方剂中的熵
    if len(herbs_list) == 1:
        return 0

    cut_num = {}
    for k2 in bins:
        cut_num[k2] = 0

    total_num = 0
    for i in range(0,len(herbs_list)-1):
        for j in range(i+1,len(herbs_list)):
            herb1 = herbs_list[i]
            herb2 = herbs_list[j]
            herb1herb2_sum = 0
            if str(herb1)+str(herb2) in herb_pair_from_data_shangshi:
                herb1herb2_sum = herb_pair_from_data_shangshi[str(herb1)+str(herb2)]
            if str(herb2)+str(herb1) in herb_pair_from_data_shangshi:
                herb1herb2_sum = herb_pair_from_data_shangshi[str(herb2)+str(herb1)]

            for k in cut_num.keys():
                if herb1herb2_sum in k:
                    cut_num[k] = cut_num[k] + 1
            total_num += 1

    # 求各个区间的熵
    total_num_entr = 0
    for k2 in cut_num.keys():
        if cut_num[k2] ==0:
            entr_cut = 0
        else:
            p = float(cut_num[k2]) / total_num
            entr_cut = -p * np.log10(p)
        total_num_entr += entr_cut
    return total_num_entr



def target_num_disease(herbs_list,tar_nodes_list,h_m_t):
    h_m_t = h_m_t[h_m_t['herb_cn_name'].isin(herbs_list)]
    herb_target = set(h_m_t['TARGET_ID'])
    tar_nodes_list = set(tar_nodes_list)

    inter = len(herb_target.intersection(tar_nodes_list))
    return inter,inter/len(herbs_list),inter/len(herb_target)


def closet_Drug_Target(herbs_list,disease_target_entrezid_list,h_m_t,path_length,gene_symbol_entrezid):

    h_m_t = h_m_t[h_m_t['herb_cn_name'].isin(herbs_list)]
    herb_target = list(h_m_t['TARGET_ID'])

    nodes_pd_target = gene_symbol_entrezid[gene_symbol_entrezid['target'].isin(herb_target)]
    nodes_symbolid = list(set(nodes_pd_target['entrezid']))

    #print(disease_target_entrezid_list)
    #print(nodes_symbolid)
    sum_path = 0
    num = 0
    for t in disease_target_entrezid_list:
        length_t = []
        for m in nodes_symbolid:
            if str(t) in path_length:
                if str(m) in path_length[str(t)]:
                    length_t.append(path_length[str(t)][str(m)])
        if len(length_t)!= 0:
            sum_path = sum_path + min(length_t)
            num = num + 1

    for m in nodes_symbolid:
        length_t = []
        for t in disease_target_entrezid_list:
            if str(m) in path_length:
                if str(t) in path_length[str(m)]:
                    length_t.append(path_length[str(m)][str(t)])
        if len(length_t)!= 0:
            sum_path = sum_path + min(length_t)
            num = num + 1
    return float(sum_path)/num



def formula_analysis_rs(filedir,nodes_pd_target,tar_nodes_list,h_m_t):
    filepath_ppiscore = 'D:/ctm_data/'
    filename_ppiscore = 'herb_PPI_score.csv'
    herbs_data = pd.read_csv(filepath_ppiscore + filename_ppiscore, encoding='ansi')

    filepath_sab = 'D:/formula_result/1算法设计/2药物配伍得分和SAB等指标的关系/'
    filename_sab = '2药物配伍得分和SAB等指标的关系.xlsx'
    herb_pair_from_data = hpff.herb_pair_score_from_Sab(filepath_sab, filename_sab, -1)

    filepath_entr = 'D:/formula_result/1算法设计/2药物配伍得分和SAB等指标的关系/'
    filename_entr = '2药物配伍得分和SAB等指标的关系.xlsx'
    herb_pair_from_data_shangshi = hpff.herb_pair_score_from_shangshi(filepath_entr, filename_entr, -1)
    herb_pair_from_data_value = np.array(list(herb_pair_from_data_shangshi.values()))
    herb_pair_from_data_cut = pd.cut(herb_pair_from_data_value,5)
    bins = herb_pair_from_data_cut.unique()

    G = di.graphFromPPI()
    path_length = dict(nx.all_pairs_shortest_path_length(G))
    gene_symbol_entrezid = di.gene_symbol_entrezid_pd()


    for root, dirs, files in os.walk(filedir, topdown=False):
        for name in files:
            filename = os.path.join(root, name)
            print(filename)
            with open(filename) as fl:
                for line in fl:
                    herbs = line.strip().split(',')
                    fmapsocre = herbs[0]
                    herbs = herbs[1:-1]

                    disease_target_symbol_list = list(set(nodes_pd_target['symbol']))
                    disease_target_entrezid_list = list(set(nodes_pd_target['entrezid']))

                    inter,inter_herbs,inter_target = target_num_disease(herbs, tar_nodes_list, h_m_t)

                    sum_row,sum_col = perturbation_formula(herbs, disease_target_symbol_list,herbs_data)
                    entr = formula_entropy(herbs,herb_pair_from_data_shangshi,bins)
                    fm_sab = fomula_sab(herbs,herb_pair_from_data)
                    closest_herbs = closet_Drug_Target(herbs,disease_target_entrezid_list,h_m_t,path_length,gene_symbol_entrezid)

                    herbs.insert(0,str(fmapsocre))
                    #herbs.insert(0,str(sum_row))
                    herbs.insert(0,str(sum_col))
                    herbs.insert(0,str(entr))
                    herbs.insert(0,str(fm_sab))
                    herbs.insert(0,str(closest_herbs))
                    herbs.insert(0,str(inter))
                    herbs.insert(0,str(inter_herbs))
                    herbs.insert(0,str(inter_target))

                    filers = name + '_rs.csv'
                    do.writedatalisttodata(os.path.join(root, filers),herbs)



if  __name__ == '__main__':
    #herbs_list = ['麻黄','桂枝','甘草','白芍','川芎']
    filepath = 'D:\\ctm_data\\TCMSP-数据\\'
    filename = 'TCMSP_DB_加工.xlsx'

    #获取对应疾病的symbolid靶点
    disease_target = di.disease_target(filepath, filename)
    nodes_list = list(disease_target)
    gene_symbol_entrezid = di.gene_symbol_entrezid_pd()
    nodes_pd_target = gene_symbol_entrezid[gene_symbol_entrezid['target'].isin(list(nodes_list))]
    #nodes_symbolid = nodes_pd_target['symbol']
    #print(nodes_symbolid)

    h_m_t = di.herb_mol_targets(filepath,filename)

    filedir = 'D:\\formula_result\\3方剂分析'
    #formula_analysis_rs(filedir, nodes_symbolid)

    formula_analysis_rs(filedir, nodes_pd_target, nodes_list, h_m_t)