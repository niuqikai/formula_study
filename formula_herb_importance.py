#算法模块，计算算法过程中的各种指标和变量
import formula_study.Data_input as di
import pandas as pd
import networkx as nx
import numpy as np
import random as rd

def jaccard_gini(targ_mol_herb):#计算jaccard距离,A为对应靶点成分的数量，药物成分的数量
    herb_mol = di.herb_molecules(filepath , filename)  # 中药对应的成分
    herb_num = herb_mol.groupby('herb_cn_name')['MOL_ID'].nunique()
    targ_mol,targ_herb = targ_mol_herb_num(targ_mol_herb)
    #疾病对应3001个成分

    disease_mol_num = targ_mol_herb['MOL_ID'].nunique()
    jac_gini_matrix = pd.merge(targ_mol.reset_index(), herb_num.reset_index() ,how= 'left', on= 'herb_cn_name')
    jac_gini_matrix['jaccard_score'] = jac_gini_matrix['MOL_ID_x'] / (jac_gini_matrix['MOL_ID_y'] + disease_mol_num - jac_gini_matrix['MOL_ID_x'] )
    jac_gini_matrix['gini_score'] = jac_gini_matrix['MOL_ID_x'] / jac_gini_matrix['MOL_ID_y'].apply(lambda x : disease_mol_num if disease_mol_num < x else x)
    rs = jac_gini_matrix.loc[:,['herb_cn_name','jaccard_score','gini_score']]
    #rs.to_csv('jaccard_gini.csv')

def pagerank_score(tag_mol_herb):#计算PageRank加权之后的每种中药得分
    G = di.Graph_from_data()
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

def Sab(G , nodes):
    distance_total = 0.0
    reduce_num = 0
    for source in nodes:
        source_list = []
        for target in nodes:
            if (source != target) and (source in G.nodes()) and (target in G.nodes()) and (G.has_edge(source ,target)):
                source_list.append(nx.shortest_path_length(G ,source ,target))
        if len(source_list) != 0:
            s_distance = np.min(source_list)
            distance_total = distance_total + s_distance
        else :
            red_num = reduce_num -1
    if len(nodes) - reduce_num != 0:
        rs = distance_total/(len(nodes) - reduce_num)
        return rs

def shortest_distance(herb_mol_target):#计算两味中药之前的平均最短路径
    filewrite = 'distance.csv'
    herbs = list(herb_mol_target['herb_cn_name'].unique())
    G = di.Graph_from_data()
    for i in range(len(herbs)-1):
        for j in range(i+1 , len(herbs) ):
            #print(herbs[i],herbs[j])
            distance_list = []
            herb1_targets = herb_mol_target[herb_mol_target['herb_cn_name'] == herbs[i]]['TARGET_ID']
            herb2_targets = herb_mol_target[herb_mol_target['herb_cn_name'] == herbs[j]]['TARGET_ID']
            herb1_targets_list = list(set(herb1_targets.dropna()))
            herb2_targets_list = list(set(herb2_targets.dropna()))
            Sa = Sab(G , herb1_targets_list)
            Sb = Sab(G , herb2_targets_list)
            print(Sa)
            print(Sb)
            herb1_herb2_targets_list = []
            herb1_herb2_targets_list.extend(herb1_targets_list)
            herb1_herb2_targets_list.extend(herb2_targets_list)
            herb1_herb2_targets_list = list(set(herb1_herb2_targets_list))
            S_ab = Sab(G , herb1_herb2_targets_list)

            for targ1 in herb1_targets_list:
                for targ2 in herb2_targets_list:
                    if targ1 in G.nodes() and targ2 in G.nodes() and G.has_edge(targ1, targ2): #同一种疾病下的成为
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
                    fw.write(",")
                    fw.write(str(S_ab - (Sa + Sb)/2.0))
                    fw.write('\n')
                    fw.flush()

def walk_score_algorithm(df_data ,source, target):#计算二分网络随机游走的分数，df_data为二分网矩阵，source target为源和目标,walk_score为随机游走得分，初始化为1
    t_m_group = dict(df_data.groupby(source)[target].nunique())
    df_data['walk_score'] = df_data.apply(lambda x: (1.0/t_m_group[x[source]] * x['walk_score'] if x[source] in t_m_group else 0) ,axis=1)
    new_walk_score = dict(df_data.groupby(target)['walk_score'].sum())
    df_data['walk_score'] = df_data.apply(lambda x: (new_walk_score[x[target]]  if x[target] in new_walk_score else 0) ,axis=1)
    return df_data

'''
def herb_walk_score(target_molecule):#计算随机游走的数据，对每个中药进行打分。
    t_m = target_molecule[['TARGET_ID','MOL_ID']].drop_duplicates()# 提取靶点和成分列，第一列为起点列
    source = 'TARGET_ID'
    target = 'MOL_ID'
    t_m['walk_score'] = t_m[source].apply(lambda x: 1.0)#初始化分数
    mols_score = walk_score_algorithm(t_m, source ,target)  #成分和对应的分数
    mols_score_dict = {key:values for key, values in zip(mols_score['MOL_ID'], mols_score['walk_score'])}#转换为字典结构

    herb_mols_values = herb_mols[['MOL_ID','herb_cn_name']].drop_duplicates()#去重
    herb_mols_values['walk_score'] = herb_mols_values['MOL_ID'].apply(lambda x: mols_score_dict[x] if x in mols_score_dict else 0)
    source = 'MOL_ID'
    target = 'herb_cn_name'
    herb_score = walk_score_algorithm(herb_mols_values, source ,target)  #药物和对应的分数
    #return herb_score

    #加权的分数算法，权重为毎味中药对应的分子数的倒数。
    h_m_v = herb_mols_values.groupby('herb_cn_name')['MOL_ID'].nunique() #计算
    herb_score_weight = pd.merge(herb_score, h_m_v.reset_index() ,how = 'left' ,on = 'herb_cn_name')#加权
    herb_score_weight['score_weight'] =  herb_score_weight.apply(lambda x : (1.0/x['MOL_ID'] * x['walk_score']),axis=1)#加权 除以毎个味中药在数据库里面的成分
    #herb_score_weight.to_csv('herb_score_weight.csv')
    return herb_score_weight
    #herb_score.to_csv('herb_score.csv')
'''

def herb_walk_score_interation(target_molecule):#计算随机游走的数据，对每个中药进行打分。迭代n次
    t_m = target_molecule[['TARGET_ID','MOL_ID']].drop_duplicates()# 提取靶点和成分列，第一列为起点列
    t_m['walk_score'] = t_m['TARGET_ID'].apply(lambda x: 1.0)  # 初始化分数
    mols_score = ''
    source = 'TARGET_ID'
    target = 'MOL_ID'
    for i in range(9):
        tm = walk_score_algorithm(t_m, source ,target)  #成分和对应的分数
        target , source = source , target
    mols_score = tm
    mols_score_dict = {key:values for key, values in zip(mols_score['MOL_ID'], mols_score['walk_score'])}#转换为字典结构

    herb_mols_values = herb_mols[['MOL_ID','herb_cn_name']].drop_duplicates()#去重
    herb_mols_values['walk_score'] = herb_mols_values['MOL_ID'].apply(lambda x: mols_score_dict[x] if x in mols_score_dict else 0)
    herb_score = ''
    source = 'MOL_ID'
    target = 'herb_cn_name'
    for i in range(9):
        herb_score = walk_score_algorithm(herb_mols_values, source ,target)  #药物和对应的分数
        target , source = source , target
    #return herb_score

    #加权的分数算法，权重为毎味中药对应的分子数的倒数。
    h_m_v = herb_mols_values.groupby('herb_cn_name')['MOL_ID'].nunique() #计算
    herb_score_weight = pd.merge(herb_score, h_m_v.reset_index() ,how = 'left' ,on = 'herb_cn_name')#加权
    herb_score_weight['score_weight'] =  herb_score_weight.apply(lambda x : (1.0/x['MOL_ID_y'] * x['walk_score']),axis=1)#加权 除以毎个味中药在数据库里面的成分
    herb_score_weight.to_csv('herb_score_weight.csv')
    return herb_score_weight
    #herb_score.to_csv('herb_score.csv')


if __name__ == '__main__':
    filepath = 'D:\\ctm_data\\TCMSP-数据\\'
    filename = 'TCMSP_DB_加工.xlsx'

    targ_mol_herb = di.targets_mol_herb(filepath,filename)#疾病靶点对应的
    herb_mol_target = di.herb_mol_targets(filepath, filename) # 计算每种中药对应的成分和靶点
    target_molecule = di.target_mol(filepath, filename, tar='0') #和疾病关联的靶点以及成分
    herb_mols =  di.herb_molecules(filepath, filename) #中药对应的成分

    shortest_distance(herb_mol_target)