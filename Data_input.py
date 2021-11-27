#此程序为数据输入接口，处理各种格式的数据并且转换为矩阵
#
import pandas as pd
import networkx as nx

def datafromcsv(fileapath):
    df = pd.DataFrame(fileapath)
    return df

def data_from_excel_sheet(filepath, st_name):
    df = pd.read_excel(filepath, sheet_name = st_name)  # 可以通过sheet_name来指定读取的表单
    return df

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


