#计算方剂中节点重要性和配伍重要性
import formula_study.Data_input as di
import pandas as pd
import networkx as nx
import numpy as np

def jaccard_gini(targ_mol_herb):#计算jaccard距离,A为对应靶点成分的数量，药物成分的数量
    herb_m = 'v_Herbs_Molecules'
    herb_mol = di.data_from_excel_sheet(filepath + filename, herb_m)  # 计算中药对应的成分
    herb_num = herb_mol.groupby('herb_cn_name')['MOL_ID'].nunique()

    targ_mol,targ_herb = targ_mol_herb_num(targ_mol_herb)

    #疾病对应3001个成分
    disease_mol_num = 3001
    jac_gini_matrix = pd.merge(targ_mol.reset_index(), herb_num.reset_index() ,how= 'left', on= 'herb_cn_name')
    jac_gini_matrix['jaccard_score'] = jac_gini_matrix['MOL_ID'] / (jac_gini_matrix['03_Info_Molecules_MOL_ID'] + disease_mol_num - jac_gini_matrix['MOL_ID'] )
    jac_gini_matrix['gini_score'] = jac_gini_matrix['MOL_ID'] / jac_gini_matrix['03_Info_Molecules_MOL_ID'].apply(lambda x : disease_mol_num if disease_mol_num < x else x)
    rs = jac_gini_matrix.loc[:,['herb_cn_name','jaccard_score','gini_score']]
    #rs.to_csv('jaccard_gimi.csv')

def Graph_from_data():# 将同一疾病的靶点连线，构成图
    filepath = 'D:\\ctm_data\\TCMSP-数据\\'
    filename = 'TCMSP_DB_加工.xlsx'
    sheet_name = 'v_Targets_Diseases'
    tag_id = 'target_ID'
    dis_id = 'disease_ID'
    G = di.data_from_excel_graph(filepath + filename, sheet_name, tag_id, dis_id)  # 将同一疾病的靶点连线，构成图
    return G

def herb_mol_targets(filepath,filename):#计算每种中药对应的成分和靶点
    herb_m = 'v_Herbs_Molecules'
    herb_mol = di.data_from_excel_sheet(filepath + filename, herb_m)  # 计算中药对应的成分
    mols_t = 'v_Molecules_Targets'
    mol_target = di.data_from_excel_sheet(filepath + filename, mols_t)  # 计算中药对应的成分

    herb_mol_target = pd.merge(herb_mol, mol_target,how = 'left',on= 'MOL_ID') #将中药 成分和成分对应的靶点进行关联
    #return herb_mol_target
    return herb_mol_target

def targets_mol_herb(filepath, filename):#生成目标靶点对应的成分和中药的矩阵
    target_m = 'target_mol'
    target_mol = di.data_from_excel_sheet(filepath + filename, target_m)  # 靶点成分对应关系

    mol_d = 'mol_disease'
    mol_disease = di.data_from_excel_sheet(filepath + filename, mol_d)  # 成分和中药的对应关系

    targ_mol_herb = pd.merge(target_mol,mol_disease,how = 'left',on= 'MOL_ID') #将疾病有关的靶点 成分 中药进行关联
    return targ_mol_herb


def pagerank_score(tag_mol_herb):#计算PageRank加权之后的每种中药得分
    G = Graph_from_data()
    pagerank_rs = nx.pagerank(G)
    herb_pagerank = {}
    for herb in tag_mol_herb['herb_cn_name'].unique():
        targets = tag_mol_herb[tag_mol_herb['herb_cn_name'] == str(herb)]['TARGET_ID']
        pagerank_s = 0
        for target in targets.unique():
            if target in pagerank_rs:
                pagerank_s = pagerank_rs[target] + pagerank_s
        herb_pagerank[herb] = pagerank_s
    r = pd.DataFrame.from_dict(herb_pagerank, orient= 'index')
    #r.to_csv('pagerankscore.csv')
    return herb_pagerank

def targ_mol_herb_num(targ_mol_herb):#计算每种药物对应的靶点成分数量和靶点数量
    targ_mol = targ_mol_herb.groupby('herb_cn_name')['MOL_ID'].nunique()#计算每种药物中 能够关联靶点的成分数量
    targ_herb = targ_mol_herb.groupby('herb_cn_name')['TARGET_ID'].nunique()#计算每种药物中 能够关联的靶点数量
    return targ_mol,targ_herb

def shortest_distance(herb_mol_target):
    filewrite = 'distance.csv'
    herbs = list(herb_mol_target['herb_cn_name'].unique())
    G = Graph_from_data()
    for i in range(len(herbs)-1):
        for j in range(i+1 , len(herbs) ):
            #print(herbs[i],herbs[j])
            distance_list = []
            herb1_targets = herb_mol_target[herb_mol_target['herb_cn_name'] == herbs[i]]['TARGET_ID']
            herb2_targets = herb_mol_target[herb_mol_target['herb_cn_name'] == herbs[j]]['TARGET_ID']
            herb1_targets_list = list(set(herb1_targets.dropna()))
            herb2_targets_list = list(set(herb2_targets.dropna()))
            for targ1 in herb1_targets_list:
                for targ2 in herb2_targets_list:
                    if targ1 in G.nodes() and targ2 in G.nodes() and G.has_edge(targ1, targ2):
                        distance_list.append(nx.shortest_path_length(G, targ1, targ2))
            if len(distance_list) !=0:
                with open(filewrite, 'a') as fw:
                    fw.write(str(herbs[i]))
                    fw.write(",")
                    fw.write(str(herbs[j]))
                    fw.write(",")
                    fw.write(str(np.min(distance_list)))
                    fw.write(",")
                    fw.write(str(np.mean(distance_list)))
                    fw.write('\n')
                    fw.flush()


if __name__ == '__main__':
    filepath = 'D:\\ctm_data\\TCMSP-数据\\'
    filename = 'TCMSP_DB_加工.xlsx'

    targ_mol_herb = targets_mol_herb(filepath,filename)#疾病靶点对应的
    herb_mol_target = herb_mol_targets(filepath, filename) # 计算每种中药对应的成分和靶点


    #for herb_name in herb_mol_target['herb_cn_name'].unique():
    #    print(herb_name)
    #    print((herb_mol_target[herb_mol_target['herb_cn_name'] == herb_name]['TARGET_ID']))

    #filewrite = 'targ_herb.csv'
    #pd.DataFrame(targ_herb).to_csv(filewrite)
    #jaccard_gini(targ_mol,herb_num)
    #pagerank_score(targ_mol_herb)



