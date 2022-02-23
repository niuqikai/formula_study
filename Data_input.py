#此程序为数据输入接口，处理各种格式的数据并且转换为矩阵或者图格式
#
import pandas as pd
import networkx as nx
import random as rd

filepath = 'D:\\ctm_data\\TCMSP-数据\\'
filename = 'TCMSP_DB_加工.xlsx'
disease_t = 'v_Targets_Diseases'
disease_file = 'D:\\network_ctm\\formula_study\\diseasename\\diseasename_HeartFailure.csv'


def datafromcsv(fileapath):
    df = pd.DataFrame(fileapath)
    return df

def data_from_excel_sheet(filepath, st_name):
    df = pd.read_excel(filepath, sheet_name = st_name)  # 可以通过sheet_name来指定读取的表单
    return df

def disease_target(filepath , filename):#根据疾病名称读取疾病，靶点矩阵,疾病名称放在diseasename.csv里面
    disease_targ = data_from_excel_sheet(filepath + filename, disease_t)  # 计算中药对应的成分

    disease_list = pd.read_csv(disease_file, sep = '#')
    print(disease_list)
    disease_tar = disease_targ[disease_targ['disease_name'].isin(list(disease_list['disease_name']))]
    return disease_tar

def target_mol(filepath , filename, tar = 'all'): #根据指定的靶点找出相对应的成分，all为默认的全量数据
    target_m = 'v_Molecules_Targets'
    target_molecule = data_from_excel_sheet(filepath + filename, target_m)

    if tar == 'all':
        return target_molecule
    else:
        dt = disease_target(filepath ,filename)
        target_molecule = target_molecule[target_molecule['TARGET_ID'].isin(list(dt['TARGET_ID']))]
        return  target_molecule

def targets_mol_herb(filepath, filename):#生成目标靶点对应的成分和中药的矩阵
    target_molecules = target_mol(filepath ,filename, 0)  # 靶点成分对应关系
    mol_herb = herb_molecules(filepath , filename)  # 成分和中药的对应关系

    targ_mol_herb = pd.merge(target_molecules, mol_herb, how = 'inner',on= 'MOL_ID') #将疾病有关的靶点 成分 中药进行关联
    #targ_mol_herb.to_csv('targ_mol_herb_left.csv')
    return targ_mol_herb

def disease_targetname(filepath , filename):#根据疾病确定靶点名称
    dt = disease_target(filepath, filename)
    return dt[['disease_name','target_name','TARGET_ID']]

def herb_molecules(filepath , filename):#计算中药和成分对应的矩阵
    herb_m = 'v_Herbs_Molecules'
    herb_mol = data_from_excel_sheet(filepath + filename, herb_m)  # 计算中药对应的成分
    return herb_mol

def Graph_from_data():# 将同一疾病的靶点连线，构成图
    sheet_name = 'v_Targets_Diseases'
    tag_id = 'TARGET_ID'
    dis_id = 'disease_ID'
    G = data_from_excel_graph(filepath + filename, sheet_name, tag_id, dis_id)  # 将同一疾病的靶点连线，构成图
    return G

def herb_mol_targets(filepath,filename):#计算每种中药对应的成分和靶点
    herb_mol = herb_molecules(filepath , filename)  # 计算中药对应的成分
    mol_target = target_mol(filepath , filename)  # 成分对应的靶点

    herb_mol_target = pd.merge(herb_mol, mol_target,how = 'inner',on= 'MOL_ID') #将中药 成分和成分对应的靶点进行关联
    #herb_mol_target.to_csv('herb_mol_target_inner.csv')
    return herb_mol_target

def targets_disease(filepath,filename):#获取所有靶点对应的疾病
    herb_m = 'v_Targets_Diseases'
    herb_mol = data_from_excel_sheet(filepath + filename, herb_m)  # 计算中药对应的成分
    return herb_mol

#todo靶点对应的疾病
def herb_mol_targets_disease(filepath,filename):#计算每种中药对应的成分和靶点
    herb_mol = herb_molecules(filepath , filename)  # 计算中药对应的成分
    mol_target = target_mol(filepath , filename)  # 成分对应的靶点
    target_disease = targets_disease(filepath,filename)#靶点对应的疾病

    herb_mol_target = pd.merge(herb_mol, mol_target,how = 'left',on= 'MOL_ID') #将中药 成分和成分对应的靶点进行关联
    herb_mol_targets_dis = pd.merge(herb_mol_target,target_disease,how = 'left',on = 'TARGET_ID')
    #return herb_mol_target
    return herb_mol_target

def data_from_excel_graph(filepath, st_name, tag_id ,disease_id):#根据Excel生成图
    #disease_ID
    #TARGET_ID
    df = pd.read_excel(filepath, st_name)
    nodes_list = list(set(df['TARGET_ID']))
    edges_list = []
    G = nx.Graph()
    #G.add_edges_from(edges_list)
    r = rd.random()
    for dis_id in df['disease_ID'].unique():
        tag_s = df[df['disease_ID'] == str(dis_id)]['TARGET_ID']

        if len(tag_s.to_list()) > 1:
            for i in range(len(tag_s.to_list()) - 2):
                for j in range(i + 1,len(tag_s.to_list()) - 1):
                    edge = (tag_s.to_list()[i] , tag_s.to_list()[j])
                    if r > 0.0:
                        edges_list.append(edge)

    G.add_edges_from(edges_list)
    return G
    #largest_cc = max(nx.connected_components(G), key=len) #最大连通子图包含的节点

#if  __name__ == '__main__':
    #herb_mol_targets(filepath, filename)