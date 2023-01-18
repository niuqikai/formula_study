import networkx as nx
import copy
import os

#根据文本中的边连接获取整个网络,返回所有节点和连边
def graphFromfile(filepath ,filename):
    G = nx.Graph()
    #如果第一行为标题
    i = 0
    with open(filepath + filename) as fl:
        for line in fl:
            if i > 0:
                lines = str(line).split('\t')
                G.add_edge(lines[0], lines[1])
            i = i + 1
    return G

#计算网络结构中的各种中心性
def getNodeimpotance(G):
    degree = nx.degree(G)
    print("degree")
    degree_rs = nx.degree_centrality(G)
    print("degree_rs")
    pagerank_rs = nx.pagerank(G)
    print("pagerank_rs")
    eigenvector_rs = nx.eigenvector_centrality(G)
    print("eigenvector_rs")
    closeness_rs = nx.closeness_centrality(G)
    print("closeness_rs")
    betweenness_rs = nx.betweenness_centrality(G)
    print("betweenness_rs")
    return degree,degree_rs,pagerank_rs,eigenvector_rs,closeness_rs,betweenness_rs

def writeNodeinfortofile(filename,r,*param):
    with open(filename, "a") as fw:
        #写入表头
        getNodeimp = ['node', 'degree', 'degree_rs', 'pagerank_rs', 'eigenvector_rs', 'closeness_rs', 'betweenness_rs']
        for nodei in getNodeimp:
            fw.write(str(nodei))
            fw.write(",")
        fw.write("\n")

        if (type(r) == nx.classes.graph.Graph):#如果是图
            for i in r.nodes():
                fw.write(str(i))
                fw.write(",")
                for rs in range(0, len(param)):
                    if rs != len(param) - 1:
                        fw.write(str(param[rs][str(i)]))
                        fw.write(",")
                    else:
                        fw.write(str(param[rs][str(i)]))
                fw.write("\n")
                fw.flush()

        if (type(r) == type({"a":"1"})):#如果是节点,dict类型
            for i in r.keys():
                fw.write(str(i))
                fw.write(str(","))
                value = r[i]
                for rs in range(0,len(value)):
                    if rs != len(value) - 1:
                        fw.write(str(value[rs]))
                        fw.write(",")
                    else:
                        fw.write(str(value[rs]))
                fw.write("\n")
                fw.flush()

def connected_components_infor(nodes, G):
    rsnodes ={}
    for i in nodes:
        Gcopy = copy.deepcopy(G)

        #出去i点 获取最大联通在子集
        G.remove_node(i)
        largest_cc = max(nx.connected_components(G), key=len)

        #获取最大连通子集以外的节点
        allnodes = copy.deepcopy(G.nodes())
        removenodes = set(allnodes) - largest_cc

        #减去不在最大连通子集的节点，即保留最大连通子集
        G.remove_nodes_from(list(removenodes))

        #计算去掉节点i之后的最大连通子集的平均最短路径和聚集系数
        aspl = nx.average_shortest_path_length(G) # 计算平均最短路径长度
        ac = nx.average_clustering(G)
        G = copy.deepcopy(Gcopy)
        rsnodes[i] = [aspl,ac]
    return rsnodes



