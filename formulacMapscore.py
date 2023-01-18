import cmapPy as cmap
import pandas as pd
import Data_input as di
import os
import Data_output as do
import PPI_analyse as ppi
import numpy as np

#这个程序使用cMap,用来计算衡量组方效果
#但是目前来看效果不好
#因此暂时保留 留待修改

#对于指定的方剂列表,基因列表，返回上调和下调基因,herbs_list为方剂，
#up_down_num为上调或者下调基因数量,disease_gene_list为基因列表,当disease_gene_list为-1时 表示为全量基因
#name为文件名，用来判别对应的节点重要性
def up_down_genes_fromPPI(herbs_list,up_down_num,disease_target_list,name):

    degree,pagerank,eigenvector,closeness,betweenness = ppi.symbol_target_sore_from_PPI()
    Nodeimpdict = {'degree': degree, 'pagerank':pagerank, 'eigenvector':eigenvector, 'closeness':closeness, 'betweenness':betweenness}


    filepath = 'data/'
    filename = 'herb_PPI_score.csv'
    herbs_score_dict = {}
    herbs_data = pd.read_csv(filepath + filename,encoding='ansi')
    herbs_data = herbs_data[herbs_data['Herb_Chinese_Name'].isin(herbs_list)]
    herbs_data_sum = herbs_data.sum()
    herbs_data_sum = pd.DataFrame(herbs_data_sum[5:],columns = ['name']).reset_index()
    herbs_data_sum.columns = ['name','score']
    herbs_score_dict = {key:values for key, values in zip(herbs_data_sum['name'], herbs_data_sum['score'])}
    herbs_score_dict_sort = sorted(herbs_score_dict.items(),key = lambda x:x[1], reverse=True)

    #根据PPI网络中的节点重要性排序
    target_importance_dict = {}
    for imp_name in Nodeimpdict.keys():
        if str(imp_name) in name:
            target_importance_dict = Nodeimpdict[imp_name]

    gene_list = []#基因列表
    if up_down_num != -1:
    #根据up_down_num返回列表
        for gn in range(up_down_num):
            (key,value) = herbs_score_dict_sort[gn]
            gene_list.append(key)

        for gn in range(up_down_num):
            (key,value) = herbs_score_dict_sort[-up_down_num:][gn]
            gene_list.append(key)

    else:
        for (key,value) in herbs_score_dict_sort:
            gene_list.append(key)

    target_dict = {}#排序后的靶点列表
    for tar in disease_target_list:
        if tar in gene_list:
            target_dict[tar] = target_importance_dict[tar]
    target_dict = sorted(target_dict.items(), key=lambda x: x[1], reverse=True)
    target_list = []
    for (k,v) in target_dict:
        target_list.append(k)

    return gene_list,target_list

#disease_target表示疾病靶点，genes列表。
def ks_a(disease_target,genes):
    ks_a_rs = []
    for dt in disease_target:
         rs = (disease_target.index(dt) + 1)/len(disease_target) - (genes.index(dt) +1)/len(genes)
         ks_a_rs.append(rs)
    return max(ks_a_rs)

def ks_b(disease_target,genes):
    ks_b_rs = []
    for dt in disease_target:
        rs = (genes.index(dt) + 1) / len(genes) - (disease_target.index(dt)) / len(disease_target)
        ks_b_rs.append(rs)
    return max(ks_b_rs)


def ES(a,b):
    if a>b:
        return a
    else:
        return -b


#计算排在top和down的靶点基因的数量
#gene_list, target_list分别为基因列表和靶点列表，num_top_down表示top或者down数量的genelist
def top_down_gene_num(gene_list, target_list,num_top_down):
    gene_top = gene_list[:num_top_down]
    gene_down = gene_list[-num_top_down:]
    gene_top.extend(gene_down)
    intersection = list(set(gene_top).intersection(set(target_list)))
    return len(intersection)

#计算靶点排在前面的均值和后面的均值排序均值
def target_rank_avg(target_list, gene_list):
    index_sum = 0
    for gn in target_list:
        top_ix =  int(gene_list.index(gn)) + 1
        down_ix = len(gene_list) - int(gene_list.index(gn))
        index_sum += np.min([top_ix,down_ix])
    return float(index_sum)/len(target_list)


#计算文件夹下的文件中方剂的ES分数,和靶点在上调下调基因的数目
def dirES_top_down_num(filedir,nodes_symbolid,up_down_num):
    for root, dirs, files in os.walk(filedir, topdown=False):
        for name in files:
            filename = os.path.join(root, name)
            print(filename)
            with open(filename) as fl:
                for line in fl:
                    herbs = line.strip().split(',')
                    fmapsocre = herbs[0]
                    herbs = herbs[1:-1]
                    gene_list, target_list = up_down_genes_fromPPI(herbs, up_down_num, list(nodes_symbolid),name)
                    a = ks_a(target_list, gene_list)
                    b = ks_b(target_list, gene_list)
                    es = ES(a, b)

                    min_avg_ix = target_rank_avg(target_list, gene_list)
                    interseclist = []
                    for num_top_down in [100,200,300,500,1000,2000]:
                        intersec = top_down_gene_num(gene_list, target_list, num_top_down)
                        interseclist.append(intersec)

                    herbs.insert(0,str(fmapsocre))
                    herbs.insert(0,str(es))
                    herbs.insert(0,str(abs(es)))
                    herbs.insert(0,str(min_avg_ix))
                    herbs.insert(0,str(interseclist[0]))
                    herbs.insert(0,str(interseclist[1]))
                    herbs.insert(0,str(interseclist[2]))
                    herbs.insert(0,str(interseclist[3]))
                    herbs.insert(0,str(interseclist[4]))
                    herbs.insert(0,str(interseclist[5]))

                    filers = name + '_rs.csv'
                    do.writedatalisttodata(os.path.join(root, filers),herbs)




if  __name__ == '__main__':
    #herbs_list = ['麻黄','桂枝','甘草','白芍','川芎']
    filepath = 'data/'
    filename = 'TCMSP_DB_加工.xlsx'

    #获取对应疾病的symbolid靶点
    disease_target = di.disease_target(filepath, filename)
    nodes_list = list(disease_target)
    gene_symbol_entrezid = di.gene_symbol_entrezid_pd()
    nodes_pd_target = gene_symbol_entrezid[gene_symbol_entrezid['target'].isin(list(nodes_list))]
    nodes_symbolid = nodes_pd_target['symbol']

    filedir = 'data/rs/'
    dirES_top_down_num(filedir,nodes_symbolid,-1)
