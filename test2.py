import networkx as nx
G = nx.Graph()
edgelist = [(0,1),(1,2)]
G.add_edges_from(edgelist)
n = dict(nx.all_pairs_shortest_path_length(G))
for k1 in n.keys():
    for k2 in dict(n[k1]):
        print(k1,k2,n[k1][k2])