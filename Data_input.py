#此程序为数据输入接口，处理各种格式的数据并且转换为矩阵或者图格式
#
import pandas as pd
import networkx as nx
import random as rd
import copy as cp

filepath = 'D:\\ctm_data\\TCMSP-数据\\'
filename = 'TCMSP_DB_加工.xlsx'
disease_t = 'v_Targets_Diseases'
disease_file = 'diseasename/diseasename_Atherosclerosis.csv'
neighbor_num = 2

def datafromcsv(fileapath):
    df = pd.DataFrame(fileapath)
    return df

def data_from_excel_sheet(filepath, st_name):
    df = pd.read_excel(filepath, sheet_name = st_name)  # 可以通过sheet_name来指定读取的表单
    return df

#根据数据生成PPI网络
def graphFromPPI():
    filepath_PPI = 'D:/ctm_data/PPI蛋白网络数据/'
    filename_PPI = 'PPI_edges.txt'
    G = nx.Graph()
    with open(filepath_PPI + filename_PPI) as fl:
        for line in fl:
            lines = str(line).split('\t')
            G.add_edge(lines[0],lines[1])
    return G

#PPI网络中的节点的一阶邻居节点和二阶邻居节点。neighbor_num表示一阶还是二阶链接。
def neighbor_symbol_PPI(G,nodes_list,neighbor_num):
    neighbor_set = set(nodes_list)
    #节点和节点的一阶链接
    for nodei in nodes_list:
        if str(nodei) in G.nodes():
            node_neighbors = set(G.neighbors(str(nodei)))
            neighbor_set = neighbor_set|node_neighbors
    #节点和节点的一阶链接

    if neighbor_num == 1:
        return  neighbor_set

    # 节点和节点的二阶链接
    if neighbor_num == 2:
        neighbor_one_set = cp.deepcopy(neighbor_set)
        for nodei in neighbor_one_set:
            if str(nodei) in G.nodes():
                node_neighbors = set(G.neighbors(str(nodei)))
                neighbor_set = neighbor_set|node_neighbors
        #节点和节点的二阶链接
        return neighbor_set

#获取目标靶点的一级二级邻居，输入为entrezid，返回target列表
def target_PPI_neighbor(G,nodes_list,neighbor_num):
    gene_symbol_entrezid = gene_symbol_entrezid_pd()
    nodes_pd_target = gene_symbol_entrezid[gene_symbol_entrezid['target'].isin(list(nodes_list))]
    nodes_entrezid = nodes_pd_target['entrezid']

    neighbor_set = neighbor_symbol_PPI(G, nodes_entrezid, neighbor_num)
    nodes_pd_entrezid = gene_symbol_entrezid[gene_symbol_entrezid['entrezid'].isin(list(neighbor_set))]
    nodes_target = nodes_pd_entrezid['target']
    return nodes_target


def disease_target(filepath , filename):#根据疾病名称读取疾病，靶点矩阵,疾病名称放在diseasename.csv里面
    disease_targ = data_from_excel_sheet(filepath + filename, disease_t)  # 计算中药对应的成分

    disease_list = pd.read_csv(disease_file, sep = '#')
    print(disease_list)
    disease_tar = disease_targ[disease_targ['disease_name'].isin(list(disease_list['disease_name']))]
    disease_tar = disease_tar['TARGET_ID']

    #一阶或者二阶链接靶点
    #G = graphFromPPI()
    #disease_tar_PPI = target_PPI_neighbor(G, disease_tar, neighbor_num)
    #disease_tar = set(disease_tar_PPI)|set(disease_tar)
    #一阶或者二阶链接靶点

    return disease_tar

def target_mol(filepath , filename, tar = 'all'): #根据指定的靶点找出相对应的成分，all为默认的全量数据
    target_m = 'v_Molecules_Targets'
    target_molecule = data_from_excel_sheet(filepath + filename, target_m)

    if tar == 'all':
        return target_molecule
    else:
        dt = disease_target(filepath ,filename)
        target_molecule = target_molecule[target_molecule['TARGET_ID'].isin(list(dt))]
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

#
def herb_mol_targets_disease(filepath,filename):#计算每种中药对应的成分和靶点
    herb_mol = herb_molecules(filepath , filename)  # 计算中药对应的成分
    mol_target = target_mol(filepath , filename)  # 成分对应的靶点
    target_disease = targets_disease(filepath,filename)#靶点对应的疾病

    herb_mol_target = pd.merge(herb_mol, mol_target,how = 'left',on= 'MOL_ID') #将中药 成分和成分对应的靶点进行关联
    herb_mol_targets_dis = pd.merge(herb_mol_target,target_disease,how = 'left',on = 'TARGET_ID')
    #return herb_mol_target
    return herb_mol_target


def targetscore():#根据PageRank算法或者度获取PPI网络中的得分
    filep = 'D:/ctm_data/'
    filen_pagerank = 'PPIpagerank.csv'
    filen_degree = 'PPIdegree.csv'
    '''
    pagerank_dict = {}
    degree_dict = {}
    with open(filep + filen_pagerank) as fl_page:
        for line in fl_page:
            lines = line.strip().split(',')
            pagerank_dict[lines[0]] = lines[1]
    with open(filep + filen_degree) as fl_degree:
        for line in fl_degree:
            lines = line.strip().split(',')
            degree_dict[line[0]] = lines[1]
    #return pagerank_dict,degree_dict
    '''
    pd_pagerank = pd.read_csv(filep + filen_pagerank)
    pd_degree = pd.read_csv(filep + filen_degree)
    return pd_pagerank,pd_degree


def targetid_SYMBOL_pd():#数据库中靶点数据和对应的geneid
    filep = 'D:/ctm_data/'
    filen = 'target_gene.csv'

    gene_symbol_dict = {}
    with open(filep + filen) as fl:
        for line in fl:
            lines = line.strip().split(',')
            gene_symbol_dict[lines[0]] = lines[1]#TAR03025,CCNA2
    #return gene_symbol_dict

    pd_gene_symbol = pd.read_csv(filep + filen)
    return pd_gene_symbol

def gene_ENTREZID_pd():
    filep = 'D:/ctm_data/'
    filen = 'SYMPOL.csv'

    gene_entrezid_dict = {} #PGK1': '5230'
    with open(filep + filen) as fl:
        for line in fl:
            lines = line.strip().split('\t')
            gene_entrezid_dict[lines[1]] = lines[2]
    #return gene_entrezid_dict

    pd_gene_entrezid = pd.read_csv(filep + filen,sep='\t')
    return pd_gene_entrezid

def gene_symbol_entrezid_pd():
    gene_entrezid_pd = gene_ENTREZID_pd()
    gene_symbol_pd = targetid_SYMBOL_pd()

    gene_symbol_entrezid = pd.merge(gene_symbol_pd, gene_entrezid_pd, how='inner',on='symbol')  #
    return gene_symbol_entrezid



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
#    print(symbol_sore_from_PPI())
    #herb_mol_targets(filepath, filename)
#    gene_symbol_entrezid_pd()
#    G = graphFromPPI()

