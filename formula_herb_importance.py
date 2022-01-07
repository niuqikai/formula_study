#算法模块，计算算法过程中的各种指标和变量
import formula_study.Data_input as di
import pandas as pd
import networkx as nx
import numpy as np
import random as rd

def herb_disease_jaccard_gini(targ_mol_herb):#计算中药和疾病之前的jaccard距离
    herb_mol = di.herb_molecules(filepath, filename)  # 中药对应的成分
    herb_num = herb_mol.groupby('herb_cn_name')['MOL_ID'].nunique()
    targ_mol, targ_herb = targ_mol_herb_num(targ_mol_herb)

    # 疾病对应成分数目
    disease_mol_num = targ_mol_herb['MOL_ID'].nunique()
    jac_gini_matrix = pd.merge(targ_mol.reset_index(), herb_num.reset_index() ,how= 'left', on= 'herb_cn_name')

    A = jac_gini_matrix['MOL_ID_x']
    B = jac_gini_matrix['MOL_ID_y']
    C = disease_mol_num

    rs = jaccard_gini(jac_gini_matrix, A, B, C)

def herb_herb_jaccard_gini(herb_mol_target):#计算中药与中药之间成分和靶点之间的
    herb_pairs = all_herb_pairs(herb_mol_target)
    mols_vector = list(herb_mol_target['MOL_ID'])
    tars_vector = list(herb_mol_target['TARGET_ID'])
    dict_herb_mol_vector,dict_herb_tar_vector,dict_herb_mol_vector_normal,dict_herb_tar_vector_normal = get_all_herb_mol_tar_vector(herb_mol_target)
    for (herb1,herb2) in herb_pairs:
        h_h_m = herb_mol_target[herb_mol_target['herb_cn_name'].isin([herb1,herb2])]
        herb_a = herb_mol_target[herb_mol_target['herb_cn_name']==herb1]
        herb_b = herb_mol_target[herb_mol_target['herb_cn_name']==herb2]

        #获取中药的成分靶点向量
        herb_a_mol_vector = dict_herb_mol_vector[herb1]
        herb_b_mol_vector = dict_herb_mol_vector[herb2]
        herb_a_tar_vector = dict_herb_tar_vector[herb1]
        herb_b_tar_vector = dict_herb_tar_vector[herb2]
        herb_a_mol_vector_normal = dict_herb_mol_vector_normal[herb1]
        herb_b_mol_vector_normal = dict_herb_mol_vector_normal[herb2]
        herb_a_mol_tar_normal = dict_herb_tar_vector_normal[herb1]
        herb_b_mol_tar_normal = dict_herb_tar_vector_normal[herb2]

        #根据成分计算各种指标
        A_mol_num = herb_a['MOL_ID'].nunique()
        B_mol_num = herb_b['MOL_ID'].nunique()
        A_B_mol_num = A_mol_num + B_mol_num - h_h_m['MOL_ID'].nunique()#中药A和中药B的交集
        mol_jaccard_herbs = float(A_B_mol_num)/h_h_m['MOL_ID'].nunique()
        mol_gini_herb = float(A_B_mol_num)/min(A_mol_num,B_mol_num)
        if herb_a_mol_vector_normal*herb_b_mol_vector_normal == 0:
            cos_mol_a_b = 0
        else:
            cos_mol_a_b = np.dot(herb_a_mol_vector,herb_b_mol_vector)/(herb_a_mol_vector_normal*herb_b_mol_vector_normal)

        #根据靶点计算各种指标
        A_tar_num = herb_a['TARGET_ID'].nunique()
        B_tar_num = herb_b['TARGET_ID'].nunique()
        A_B_tar_num = A_tar_num + B_tar_num - h_h_m['TARGET_ID'].nunique()  # 中药A和中药B的交集
        if float(A_B_tar_num) == 0:
            tar_jaccard_herbs = 0
            tar_gini_herb = 0
        else:
            tar_jaccard_herbs = float(A_B_tar_num) / h_h_m['TARGET_ID'].nunique()
            tar_gini_herb = float(A_B_tar_num) / min(A_tar_num, B_tar_num)
        if herb_a_mol_tar_normal*herb_b_mol_tar_normal == 0:
            cos_tar_a_b = 0
        else:
            cos_tar_a_b = np.dot(herb_a_tar_vector,herb_b_tar_vector)/(herb_a_mol_tar_normal*herb_b_mol_tar_normal)
        datalist = [herb1,herb2,mol_jaccard_herbs,mol_gini_herb,tar_jaccard_herbs,tar_gini_herb,cos_mol_a_b,cos_tar_a_b]
        filename = 'herb_herb_walkscore_mol_jaccard_gini.csv'
        writelisttodata(filename, datalist)

def get_all_herb_mol_tar_vector(herb_mol_target):#毎味中药的成分和靶点向量
    herbs_name = list(herb_mol_target['herb_cn_name'].unique())#所有的中药名称
    mols_vector = list(herb_mol_target['MOL_ID'])
    tars_vector = list(herb_mol_target['TARGET_ID'])
    dict_herb_mol_vector = {}
    dict_herb_tar_vector = {}
    dict_herb_mol_vector_normal = {}
    dict_herb_tar_vector_normal = {}

    #如果计算中药向量，使用随机随机游走方法则加上这段
    #成分
    t_m = herb_mol_target[['herb_cn_name','MOL_ID']].drop_duplicates()# 提取靶点和成分列，第一列为起点列
    t_m['walk_score'] = t_m['herb_cn_name'].apply(lambda x: 1.0)  # 初始化分数
    df_data = walk_score_algorithm(t_m, 'herb_cn_name', 'MOL_ID')
    h_w_dict = {key:values for key, values in zip(df_data['MOL_ID'], df_data['walk_score'])}#转换为字典结构
    #靶点
    m_t = herb_mol_target[['MOL_ID','TARGET_ID']].drop_duplicates()#去重
    m_t['walk_score'] = m_t['MOL_ID'].apply(lambda x: h_w_dict[x] if x in h_w_dict else 0)
    df_data = walk_score_algorithm(m_t, 'MOL_ID', 'TARGET_ID')
    m_t_dict = {key:values for key, values in zip(df_data['TARGET_ID'], df_data['walk_score'])}
    #如果计算中药向量，使用随机随机游走方法则加上这段

    for i in range(len(herbs_name)):
        herb_mol_vector = [0 for _ in range(len(mols_vector))]
        herb_tar_vector = [0 for __ in range(len(tars_vector))]
        herb_i = herb_mol_target[herb_mol_target['herb_cn_name']==herbs_name[i]]
        herb_i_mols = list(herb_i['MOL_ID'])
        herb_i_tars = list(herb_i['TARGET_ID'])

        for j in range(len(herb_mol_vector)):
            if mols_vector[j] in herb_i_mols:
                herb_mol_vector[j] = h_w_dict[mols_vector[j]] #1 #随机游走分数，否则默认为1
        for k in range(len(herb_tar_vector)):
            if tars_vector[k] in herb_i_tars:
                herb_tar_vector[k] = m_t_dict[tars_vector[k]] #1 #随机游走分数，否则默认为1
        dict_herb_mol_vector[herbs_name[i]] = herb_mol_vector
        dict_herb_tar_vector[herbs_name[i]] = herb_tar_vector
        dict_herb_mol_vector_normal[herbs_name[i]] = np.linalg.norm(herb_mol_vector)
        dict_herb_tar_vector_normal[herbs_name[i]] = np.linalg.norm(herb_tar_vector)
    return dict_herb_mol_vector,dict_herb_tar_vector,dict_herb_mol_vector_normal,dict_herb_tar_vector_normal


def jaccard_gini(jac_gini_matrix,A,B,C):#计算jaccard距离,A为交集，B，C为两个集合
    jac_gini_matrix['jaccard_score'] = A / (B + C - A)
    jac_gini_matrix['gini_score'] = A / B.apply(lambda x : C if C < x else x)
    rs = jac_gini_matrix.loc[:,['herb_cn_name','jaccard_score','gini_score']]
    return  rs
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

def writelisttodata(filename , datalist):#将列表数据写入文本
    with open(filename,'a') as fl:
        for dl in range(len(datalist) - 1):
            fl.write(str(datalist[dl]))
            fl.write(',')
        fl.write(str(datalist[-1]))
        fl.write('\n')
        fl.flush()

def all_herb_pairs(herb_mol_target):
    herb_pairs = []
    herbs = list(herb_mol_target['herb_cn_name'].unique())
    for i in range(len(herbs) - 1):
        for j in range(i + 1, len(herbs)):
            herb_pairs.append((herbs[i],herbs[j]))
    return  herb_pairs

#def write_to_file(filename , valuelist):
def formula_generate_algorithm(herb_score_dict , pair_score_dict):#生成组方的核心算法
    pair_seed = pair_score_dict.keys()

    for (herb1 ,herb2) in pair_seed: #基于一组药对，衍生一张方子
        formula_list = [herb1 ,herb2]
        before_score = herb_score_dict[herb1] + herb_score_dict[herb2] + pair_score_dict[(herb1,herb2)]

    #for h in herb_score_dict.keys():  #基于一味中药衍生一张方子
    #    formula_list = [h]
    #    before_score = herb_score_dict[h]

        max_score = 0 #组成的方剂的分数最大值
        insert_herb = ''
        max_herb_score = -99999 # 下一个插入的中药中，最大的分数
        while(len(formula_list) < 11):
            for herb in herb_score_dict.keys():
                if herb not in formula_list:
                    herb_score_insert = herb_score_dict[herb]#新加入的中药分数
                    pair_score_list = []
                    pair_score = 0
                    for hb in formula_list :
                        if (herb, hb) in pair_score_dict:#药物组合的分数
                            #pair_score_list.append(pair_score_dict[(herb,hb)])
                            pair_score = pair_score + pair_score_dict[(herb,hb)]
                        if (hb, herb) in pair_score_dict:
                            pair_score = pair_score + pair_score_dict[(hb,herb)]
                    #pair_score = pair_score / len(formula_list) #****设定组方分数为新加入的中药与其他中药的均值，若为求和则去掉
                    '''
                    if len(pair_score_list) == 0:
                        pair_score = 0
                    else:
                        pair_score = np.max(pair_score_list)
                    '''
                    if before_score + herb_score_insert - pair_score > max_herb_score :
                        max_herb_score = before_score + herb_score_insert - pair_score
                        insert_herb = herb
            if max_herb_score > before_score:
                formula_list.append(insert_herb)
                before_score = max_herb_score
                max_herb_score = -99999
            else:
                break
        formula_list.append(before_score)
        writelisttodata('_formula_pair_cos_score.csv',formula_list)


def walk_score_algorithm(df_data ,source, target):#计算二分网络随机游走的分数，df_data为二分网矩阵，source target为源和目标,walk_score为随机游走得分，初始化为1
    t_m_group = dict(df_data.groupby(source)[target].nunique())
    df_data['walk_score'] = df_data.apply(lambda x: (1.0/t_m_group[x[source]] * x['walk_score'] if x[source] in t_m_group and t_m_group[x[source]] !=0 else 0) ,axis=1)
    new_walk_score = dict(df_data.groupby(target)['walk_score'].sum())
    df_data['walk_score'] = df_data.apply(lambda x: (new_walk_score[x[target]]  if x[target] in new_walk_score else 0) ,axis=1)
    return df_data


def herb_walk_score_interation(target_molecule,herb_mols):#计算随机游走的数据，对每个中药进行打分。迭代n次
    t_m = target_molecule[['TARGET_ID','MOL_ID']].drop_duplicates()# 提取靶点和成分列，第一列为起点列
    t_m['walk_score'] = t_m['TARGET_ID'].apply(lambda x: 1.0)  # 初始化分数
    mols_score = ''
    source = 'TARGET_ID'
    target = 'MOL_ID'
    for i in range(1):
        tm = walk_score_algorithm(t_m, source ,target)  #成分和对应的分数
        target , source = source , target
    mols_score = tm
    mols_score_dict = {key:values for key, values in zip(mols_score['MOL_ID'], mols_score['walk_score'])}#转换为字典结构

    herb_mols_values = herb_mols[['MOL_ID','herb_cn_name']].drop_duplicates()#去重
    herb_mols_values['walk_score'] = herb_mols_values['MOL_ID'].apply(lambda x: mols_score_dict[x] if x in mols_score_dict else 0)
    herb_score = ''
    source = 'MOL_ID'
    target = 'herb_cn_name'
    for i in range(1):
        herb_score = walk_score_algorithm(herb_mols_values, source ,target)  #药物和对应的分数
        target , source = source , target
    #return herb_score

    #加权的分数算法，权重为毎味中药对应的分子数的倒数。
    h_m_v = herb_mols_values.groupby('herb_cn_name')['MOL_ID'].nunique() #计算
    herb_score_weight = pd.merge(herb_score, h_m_v.reset_index() ,how = 'left' ,on = 'herb_cn_name')#加权
    herb_score_weight['score_weight'] =  herb_score_weight.apply(lambda x : (1.0/x['MOL_ID_y'] * x['walk_score']),axis=1)#加权 除以毎个味中药在数据库里面的成分
    #herb_score_weight.to_csv('herb_score_weight.csv')
    return herb_score_weight
    #herb_score.to_csv('herb_score.csv')


if __name__ == '__main__':
    filepath = 'D:\\ctm_data\\TCMSP-数据\\'
    filename = 'TCMSP_DB_加工.xlsx'

    targ_mol_herb = di.targets_mol_herb(filepath,filename)#疾病靶点对应的成分和中药
    herb_mol_target = di.herb_mol_targets(filepath, filename) # 计算每种中药对应的成分和靶点
    target_molecule = di.target_mol(filepath, filename, tar='0') #和疾病关联的靶点以及成分
    herb_mols =  di.herb_molecules(filepath, filename) #中药对应的成分
    hmtd = di.herb_mol_targets_disease(filepath,filename)

    filepath = 'D:\\network_ctm\\formula_study\\data\\'
    filename = 'result.xlsx'
    herb_score = 'herb_score_9times'
    pair_score = 'herb_herb_mol_jaccard_gini'

    herb_s = di.data_from_excel_sheet(filepath + filename,herb_score)
    h_s = herb_s[['herb_cn_name','walk_score']]
    h_s_dict = {key:values for key, values in zip(h_s['herb_cn_name'], h_s['walk_score'])}#转换为字典结构

    '''
    #根据sab计算药物配伍得分
    pair_s = di.data_from_excel_sheet(filepath + filename, pair_score)
    p_s = pair_s[['herb1','herb2','sab']]
    p_s_dict = {(key1 ,key2):values for key1, key2 ,values in zip(p_s['herb1'], p_s['herb2'], p_s['sab'])}#转换为字典结构
    '''

    #根据余弦夹角计算药物组合得分
    pair_s = di.data_from_excel_sheet(filepath + filename, pair_score)
    p_s = pair_s[['herb1','herb2','cos_mol']]
    p_s_dict = {(key1 ,key2):values for key1, key2 ,values in zip(p_s['herb1'], p_s['herb2'], p_s['cos_mol'])}#转换为字典结构

    #herb_herb_jaccard_gini(herb_mol_target)

