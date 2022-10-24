# -*- coding: utf-8 -*-
"""
Created on Sat Aug  6 16:11:26 2022

@author: Administrator
"""
import networkx as nx
import numpy as np
import pandas as pd
import random as rd
import os
#计算PPI网络中的sab值和dab值
def graphFromPPI():
    filepath_PPI = 'C:\\Users\\Administrator\\Desktop\\'
    filename_PPI = 'PPI1.txt'
    G = nx.Graph()
    with open(filepath_PPI + filename_PPI) as fl:
        for line in fl:
            lines = str(line).strip().split('\t')
            G.add_edge(lines[0], lines[1])
    H = nx.subgraph(G, max(nx.connected_components(G), key=len))
    return H


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

#找出到geneA和geneB的最短距离基因排序
def shortest_nodes(G, genesA, genesB, path_length):
    nodes_shortest = {}
    for node in G.nodes():
        lengh = 0
        for nodei in genesA:
            if node in path_length and nodei in path_length[node]:
                lengh = lengh + path_length[node][nodei]
            else:
                lengh = lengh + 10000
        for nodej in genesB:
            if node in path_length and nodej in path_length[node]:
                lengh = lengh + path_length[node][nodej]
            else:
                lengh = lengh + 10000
        nodes_shortest[node] = lengh
    nodes_shortest = sorted(nodes_shortest.items(), key = lambda x:x[1], reverse= False)
    return nodes_shortest

#计算fc节点
def fc(G, genesA, genesB):
    node_dict = {}
    for nodei in genesA:
        for nodej in genesB:
            if nodei != nodej :
                splist = [p for p in nx.all_shortest_paths(G,nodei,nodej)]
                nodes_num = {}
                for lst in splist:
                    for n in lst:
                        if n!=nodei and n!=nodej:
                            if n in nodes_num:
                                nodes_num[n] = nodes_num[n] + 1
                            else:
                                nodes_num[n] = 1

                for k in nodes_num.keys():
                    nodes_num[k] = float(nodes_num[k])/len(splist)
                for k in nodes_num.keys():
                    if k in node_dict:
                        node_dict[k] = node_dict[k] + nodes_num[k]
                    else:
                        node_dict[k] = nodes_num[k]
    node_dict = sorted(node_dict.items(), key = lambda x:x[1], reverse= False)
    return node_dict


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
#shortest_path = dict(nx.all_shortest_paths(G))
filepathgene = 'C:\\Users\\Administrator\\Desktop\\gene1\\'
filepathtar = 'C:\\Users\\Administrator\\Desktop\\tar1\\'

#filegene1 = 'gene1.txt'
#filegene2 = 'gene2.txt'
#targets = 'targets.txt'
genefilelist = []
tarfilelist = []
for root, dirs, files in os.walk(filepathgene):

    for f in files:
        filename = os.path.join(root, f)
        genefilelist.append(filename)

for root, dirs, files in os.walk(filepathtar):
    for f in files:
        filename = os.path.join(root, f)
        tarfilelist.append(filename)

#计算sab
for i in range(0,len(genefilelist)-1):
    for j in range(i+1,len(genefilelist)):
        gene1_pd = pd.read_csv(genefilelist[i])
        herb1_targets_list = list(gene1_pd['gene'].apply(lambda x: str(x)))
        gene2_pd = pd.read_csv(genefilelist[j])
        herb2_targets_list = list(gene2_pd['gene'].apply(lambda x: str(x)))
        S_ab = Sab(G, herb1_targets_list, herb2_targets_list, path_length)
        Sa = Saa(G, herb1_targets_list, path_length)
        Sb = Saa(G, herb2_targets_list, path_length)
        print(genefilelist[i],genefilelist[j],'S_ab', S_ab - (Sa + Sb) / 2.0)

#计算zab
for k in genefilelist:
    for m in tarfilelist:
        g_pd = pd.read_csv(k)
        herb_targets_list = list(g_pd['gene'].apply(lambda x: str(x)))

        tars = pd.read_csv(m)
        targets = list(tars['targets'].apply(lambda x :str(x)))
        dab = Sab(G, herb_targets_list, targets, path_length)
        #print(str(k),str(m),'dab:',dab)
        mean_rs,var_rs = zdct(herb_targets_list,targets,path_length)
        zd = (dab - mean_rs)/var_rs
        print(str(k),str(m),'zd:',zd)


#计算关键基因排序
for i in range(0,len(genefilelist)-1):
    for j in range(i+1,len(genefilelist)):
        gene1_pd = pd.read_csv(genefilelist[i])
        herb1_targets_list = list(gene1_pd['gene'].apply(lambda x: str(x)))
        gene2_pd = pd.read_csv(genefilelist[j])
        herb2_targets_list = list(gene2_pd['gene'].apply(lambda x: str(x)))
        rs_s = shortest_nodes(G, herb1_targets_list, herb2_targets_list, path_length)

        filen = str(i) + str(j) + '_shortest.csv'
        with open(filen,'a') as f:
            for (k,v) in rs_s:
                f.write(str(k))
                f.write(',')
                f.write(str(v))
                f.write('\n')
                f.flush()



#计算FC
for i in range(0,len(genefilelist)-1):
    for j in range(i+1,len(genefilelist)):
        gene1_pd = pd.read_csv(genefilelist[i])
        herb1_targets_list = list(gene1_pd['gene'].apply(lambda x: str(x)))
        gene2_pd = pd.read_csv(genefilelist[j])
        herb2_targets_list = list(gene2_pd['gene'].apply(lambda x: str(x)))
        result_rs = fc(G, herb1_targets_list, herb2_targets_list)
        print(i,genefilelist[i])
        print(j,genefilelist[j])
        filen = str(i) + str(j) + '_fc.csv'
        with open(filen, 'a') as f:
            for (k, v) in result_rs:
                f.write(str(k))
                f.write(',')
                f.write(str(v))
                f.write('\n')
                f.flush()
