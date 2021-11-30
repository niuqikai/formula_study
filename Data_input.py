#此程序为数据输入接口，处理各种格式的数据并且转换为矩阵或者图格式
#
import pandas as pd
import networkx as nx

def datafromcsv(fileapath):
    df = pd.DataFrame(fileapath)
    return df

def data_from_excel_sheet(filepath, st_name):
    df = pd.read_excel(filepath, sheet_name = st_name)  # 可以通过sheet_name来指定读取的表单
    return df

def disease_target(filepath , filename):#根据疾病名称读取疾病，靶点矩阵,疾病名称放在diseasename.csv里面
    disease_t = 'v_Targets_Diseases'
    disease_targ = data_from_excel_sheet(filepath + filename, disease_t)  # 计算中药对应的成分
    #'disease_name'
    disease_file = 'diseasename.csv'
    disease_list = pd.read_csv(disease_file, sep = ';')
    disease_tar =disease_targ[disease_targ['disease_name'].isin(list(disease_list['disease_name']))]
    return disease_tar

def target_mol(filepath , filename, tar = 'all'): #根据指定的靶点找出相对应的成分，all为默认的全量数据
    target_m = 'v_Molecules_Targets'
    target_molecule = data_from_excel_sheet(filepath + filename, target_m)

    if tar == 'all':
        return target_molecule
    else:
        dt = disease_target(filepath ,filename)
        target_molecule = target_molecule[target_molecule['TARGET_ID'].isin(list(dt['target_ID']))]
        return  target_molecule

def targets_mol_herb(filepath, filename):#生成目标靶点对应的成分和中药的矩阵
    target_molecules = target_mol(filepath ,filename, 0)  # 靶点成分对应关系
    mol_herb = herb_molecules(filepath , filename)  # 成分和中药的对应关系

    targ_mol_herb = pd.merge(target_molecules, mol_herb, how = 'left',on= 'MOL_ID') #将疾病有关的靶点 成分 中药进行关联
    return targ_mol_herb

def herb_molecules(filepath , filename):#计算中药和成分对应的矩阵
    herb_m = 'v_Herbs_Molecules'
    herb_mol = data_from_excel_sheet(filepath + filename, herb_m)  # 计算中药对应的成分
    return herb_mol

def Graph_from_data():# 将同一疾病的靶点连线，构成图
    filepath = 'D:\\ctm_data\\TCMSP-数据\\'
    filename = 'TCMSP_DB_加工.xlsx'
    sheet_name = 'v_Targets_Diseases'
    tag_id = 'target_ID'
    dis_id = 'disease_ID'
    G = data_from_excel_graph(filepath + filename, sheet_name, tag_id, dis_id)  # 将同一疾病的靶点连线，构成图
    return G

def herb_mol_targets(filepath,filename):#计算每种中药对应的成分和靶点
    herb_mol = herb_molecules(filepath , filename)  # 计算中药对应的成分
    mol_target = target_mol(filepath , filename)  # 成分对应的靶点

    herb_mol_target = pd.merge(herb_mol, mol_target,how = 'left',on= 'MOL_ID') #将中药 成分和成分对应的靶点进行关联
    #return herb_mol_target
    return herb_mol_target

def data_from_excel_graph(filepath, st_name, tag_id ,disease_id):#根据Excel生成图
    #disease_ID
    #TARGET_ID
    df = pd.read_excel(filepath, st_name)
    nodes_list = list(set(df['target_ID']))
    edges_list = []
    G = nx.Graph()
    #G.add_edges_from(edges_list)

    for dis_id in df['disease_ID'].unique():
        tag_s = df[df['disease_ID'] == str(dis_id)]['target_ID']

        if len(tag_s.to_list()) > 1:
            for i in range(len(tag_s.to_list()) - 2):
                for j in range(i + 1,len(tag_s.to_list()) - 1):
                    edge = (tag_s.to_list()[i] , tag_s.to_list()[j])
                    edges_list.append(edge)

    G.add_edges_from(edges_list)
    return G
    #largest_cc = max(nx.connected_components(G), key=len) #最大连通子图包含的节点


