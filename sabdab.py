# -*- coding: utf-8 -*-
"""
Created on Sat Aug  6 16:11:26 2022

@author: Administrator
"""
import networkx as nx
import numpy as np
import pandas as pd
import random as rd

def graphFromPPI():
    filepath_PPI = 'D:/ctm_data/PPI蛋白网络数据/'
    filename_PPI = 'PPI_edges.txt'
    G = nx.Graph()
    with open(filepath_PPI + filename_PPI) as fl:
        for line in fl:
            lines = str(line).split('\t')
            G.add_edge(lines[0], lines[1])
    return G


def Saa(G, nodes, path_length):
    distance_total = 0.0
    reduce_num = 0
    for source in nodes:
        source_list = []
        for target in nodes:
            if (source != target) and (source in G.nodes()) and (target in G.nodes()) and (source in path_length):
                if target in path_length[source]:
                    source_list.append(path_length[source][target])
        if len(source_list) != 0:
            s_distance = np.min(source_list)
            distance_total = distance_total + s_distance
            reduce_num = reduce_num + 1
    if reduce_num != 0:
        rs = float(distance_total) / reduce_num
        return rs


def Sab(G, nodesA, nodesB, path_length):
    distance_total = 0.0
    reduce_num = 0
    for source in nodesA:
        source_list = []
        for target in nodesB:
            if (source in G.nodes()) and (target in G.nodes() and (source in path_length)):
                if target in path_length[source]:
                    source_list.append(path_length[source][target])
        if len(source_list) != 0:
            s_distance = np.min(source_list)
            distance_total = distance_total + s_distance
            reduce_num = reduce_num + 1

    for source in nodesB:
        source_list = []
        for target in nodesA:
            if (source in G.nodes()) and (target in G.nodes() and (source in path_length)):
                if target in path_length[source]:
                    source_list.append(path_length[source][target])
        if len(source_list) != 0:
            s_distance = np.min(source_list)
            distance_total = distance_total + s_distance
            reduce_num = reduce_num + 1

    if reduce_num != 0:
        rs = float(distance_total) / reduce_num
        return rs

def zdct(herb1_targets_list,targets,path_length):
    degrees = nx.degree(G)
    d_dict = dict(degrees)
    tar1_degree_distr = {}#度分布
    tar2_degree_distr = {}#度分布

    network_degrees_distr = {}

    #获得靶点度分布
    for h in herb1_targets_list:
        h_degree = d_dict[h]
        if h_degree in tar1_degree_distr:
            tar1_degree_distr[h_degree] = tar1_degree_distr[h_degree] + 1
        else:
            tar1_degree_distr[h_degree] =  1

    for h in targets:
        h_degree = d_dict[h]
        if h_degree in tar2_degree_distr:
            tar2_degree_distr[h_degree] = tar2_degree_distr[h_degree] + 1
        else:
            tar2_degree_distr[h_degree] =  1

    #获得网络度
    for (k,v) in degrees:
        if v in network_degrees_distr:
            rs = network_degrees_distr[v]
            rs.append(k)
            network_degrees_distr[v] = rs
        else :
            network_degrees_distr[v] = [k]

    rs_list = []
    for i in range(1000):
        tars1_list = []
        tars2_list = []
        for (k,v) in tar1_degree_distr.items():
            if k < 11:
                tars1 = rd.sample(network_degrees_distr[k],v)
                tars1_list.extend(tars1)
            else:
                tmp = []
                for m in range(k - 10,k + 11):
                    if m in network_degrees_distr:
                        tmp.extend(network_degrees_distr[m])
                tars1 = rd.sample(tmp, v)
                tars1_list.extend(tars1)

        for (k,v) in tar2_degree_distr.items():
            if k < 11:
                tars2 = rd.sample(network_degrees_distr[k], v)
                tars2_list.extend(tars2)
            else:
                tmp = []
                for m in range(k - 10, k + 11):
                    if m in network_degrees_distr:
                        tmp.extend(network_degrees_distr[m])
                tars2 = rd.sample(tmp, v)
                tars2_list.extend(tars2)
        rs_list.append(Sab(G, tars1_list, tars2_list, path_length))
    return np.mean(rs_list),np.var(rs_list)

G = graphFromPPI()

'''
degree_dict = 'degree_dict.csv'
path_length = dict(nx.all_pairs_shortest_path_length(G))
with open(degree_dict,'a') as f:
    for (k,v) in path_length.items():
        for (m,n) in v.items():
            f.write(str(k))
            f.write(',')
            f.write(str(m))
            f.write(',')
            f.write(str(n))
            f.write('\n')
            f.flush()

path_length = {}
dict_file = 'degree_dict.csv'
with open(dict_file) as d:
    for line in d:
        lines = line.split(',')
        path_length[lines[0]] = lines[1]
'''
path_length = dict(nx.all_pairs_shortest_path_length(G))
filepath = 'C:\\Users\\Administrator\\Desktop\\'
filegene1 = 'gene1.txt'
filegene2 = 'gene2.txt'
targets = 'targets.txt'

gene1_pd = pd.read_csv(filepath + filegene1)
herb1_targets_list = list(gene1_pd['gene1'].apply(lambda x: str(x)))
gene2_pd = pd.read_csv(filepath + filegene2)
herb2_targets_list = list(gene2_pd['gene2'].apply(lambda x: str(x)))
tars = pd.read_csv(filepath + targets)
targets = list(tars['targets'].apply(lambda x :str(x)))


S_ab = Sab(G, herb1_targets_list, herb2_targets_list, path_length)
Sa = Saa(G, herb1_targets_list, path_length)
Sb = Saa(G, herb2_targets_list, path_length)
print(S_ab - (Sa + Sb) / 2.0)


mean_rs,var_rs = zdct(herb1_targets_list,targets,path_length)
zd = (S_ab - mean_rs)/var_rs
print(zd)
