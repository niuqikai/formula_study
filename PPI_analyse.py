import networkx as nx
import Data_output as do
import node_importance as ni
import Data_input as di
import pandas as pd
import numpy as np
import networkx.algorithms.community as nx_comm

def graphFromPPI(filepath ,filename):
    G = nx.Graph()
    with open(filepath + filename) as fl:
        for line in fl:
            lines = str(line).split('\t')
            G.add_edge(lines[0],lines[1])
    return G


def nodeimportance(G):
    degree, degree_rs, pagerank_rs, eigenvector_rs, closeness_rs, betweenness_rs = ni.getNodeimpotance(G)
    return degree, degree_rs, pagerank_rs, eigenvector_rs, closeness_rs, betweenness_rs

#节点的一度链接，输入参数为entrezid,网络为PPI,输出为entrezid
#def link_nodes_G(nodes):


def symbol_sore_from_PPI():#根据PPI网络获取节点，返回SYMBOL和得分，示例TAR03276：0.00024348708566608384
    filepath = 'data/'
    filename = 'PPI_node_importance.csv'
    importance_file_ppi = pd.read_csv(filepath + filename)
    #Entrez,degree, degree_rs, pagerank_rs, eigenvector_rs, closeness_rs, betweenness_rs

    gene_entrezid_pd = di.gene_ENTREZID_pd()
    gene_symbol_pd = di.targetid_SYMBOL_pd()

    gene_symbol_entrezid = pd.merge(gene_symbol_pd, gene_entrezid_pd, how='inner',on='symbol')  #
    gene_symbol_entrezid_importance = pd.merge(gene_symbol_entrezid, importance_file_ppi, how='inner',on = 'entrezid')

    #gene_symbol_entrezid_importance.to_csv('ppiimportance.csv')
    #return gene_symbol_entrezid_page,gene_symbol_entrezid_degree
    gene_symbol_entrezid_degree_dict = {key:values for key, values in zip(gene_symbol_entrezid_importance['target'], gene_symbol_entrezid_importance['degree'])}
    gene_symbol_entrezid_pagerank_dict = {key:values for key, values in zip(gene_symbol_entrezid_importance['target'], gene_symbol_entrezid_importance['pagerank_rs'])}
    gene_symbol_entrezid_eigenvector_dict = {key:values for key, values in zip(gene_symbol_entrezid_importance['target'], gene_symbol_entrezid_importance['eigenvector_rs'])}
    gene_symbol_entrezid_closeness_dict = {key:values for key, values in zip(gene_symbol_entrezid_importance['target'], gene_symbol_entrezid_importance['closeness_rs'])}
    gene_symbol_entrezid_betweenness_dict = {key:values for key, values in zip(gene_symbol_entrezid_importance['target'], gene_symbol_entrezid_importance['betweenness_rs'])}

    return  gene_symbol_entrezid_degree_dict,gene_symbol_entrezid_pagerank_dict,gene_symbol_entrezid_eigenvector_dict,gene_symbol_entrezid_closeness_dict,gene_symbol_entrezid_betweenness_dict


def symbol_target_sore_from_PPI():#根据PPI网络获取节点，返回SYMBOL和得分，示例CCNA2：0.00024348708566608384
    filepath = 'data/'
    filename = 'PPI_node_importance.csv'
    importance_file_ppi = pd.read_csv(filepath + filename)
    #Entrez,degree, degree_rs, pagerank_rs, eigenvector_rs, closeness_rs, betweenness_rs

    gene_entrezid_pd = di.gene_ENTREZID_pd()
    gene_symbol_pd = di.targetid_SYMBOL_pd()

    gene_symbol_entrezid = pd.merge(gene_symbol_pd, gene_entrezid_pd, how='inner',on='symbol')  #
    gene_symbol_entrezid_importance = pd.merge(gene_symbol_entrezid, importance_file_ppi, how='inner',on = 'entrezid')
    #gene_symbol_entrezid_importance.to_csv('ppiimportance.csv')
    #return gene_symbol_entrezid_page,gene_symbol_entrezid_degree
    gene_symbol_entrezid_degree_dict = {key:values for key, values in zip(gene_symbol_entrezid_importance['symbol'], gene_symbol_entrezid_importance['degree'])}
    gene_symbol_entrezid_pagerank_dict = {key:values for key, values in zip(gene_symbol_entrezid_importance['symbol'], gene_symbol_entrezid_importance['pagerank_rs'])}
    gene_symbol_entrezid_eigenvector_dict = {key:values for key, values in zip(gene_symbol_entrezid_importance['symbol'], gene_symbol_entrezid_importance['eigenvector_rs'])}
    gene_symbol_entrezid_closeness_dict = {key:values for key, values in zip(gene_symbol_entrezid_importance['symbol'], gene_symbol_entrezid_importance['closeness_rs'])}
    gene_symbol_entrezid_betweenness_dict = {key:values for key, values in zip(gene_symbol_entrezid_importance['symbol'], gene_symbol_entrezid_importance['betweenness_rs'])}

    return  gene_symbol_entrezid_degree_dict,gene_symbol_entrezid_pagerank_dict,gene_symbol_entrezid_eigenvector_dict,gene_symbol_entrezid_closeness_dict,gene_symbol_entrezid_betweenness_dict


if  __name__ == '__main__':
    #symbol_target_sore_from_PPI()
    G = di.graphFromPPI()  # PPI网络

    nodesis = ['10499','8648','367','4306','2099','2908']
    nodesjs = ['10499','6531','6532','1803','5644','146','147','1812',
              '1432','152','154','153','801','43','5292','3757','5294',
              '2099','2100','6331','3778','1111','5468','2908','2147','1128',
              '4843','1131','1133','5742','5743','367','4846','2159','6256','2932',
              '890','4988','3326']
    num = 0
    count = 0
    for nodei in nodesis:
        dis_list = []
        for nodej in nodesjs:
            dis_list.append(nx.shortest_path_length(G,nodei,nodej))
        print(dis_list)
        count = count + np.min(dis_list)
        num = num + 1

    for nodei in nodesjs:
        dis_list = []
        for nodej in nodesis:
            dis_list.append(nx.shortest_path_length(G,nodei,nodej))
        print(dis_list)
        count = count + np.min(dis_list)
        num = num + 1
    print(count/num)
    '''
    path_length = dict(nx.all_pairs_shortest_path_length(G))
    all_pairs_path = 'path_length.csv'
    with open(all_pairs_path,'a') as f:
        for k1 in path_length.keys():
            print(k1)
            for k2 in dict(path_length[k1]).keys():
                f.write(str(k1))
                f.write(',')
                f.write(str(k2))
                f.write(',')
                f.write(str(k1)+str(k2))
                f.write(',')
                f.write(str(path_length[k1][k2]))
                f.write('\n')
                f.flush()
'''

