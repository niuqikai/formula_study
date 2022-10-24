import formula_herb_importance as fhi
import Data_input as di
import Data_output as do
import herb_pairs_from_formula as hpff
import random as rd
import formula_generate as fg
import pandas as pd
import PPI_analyse as ppi

def formula_tareget(herb_mol_target,herbs):#返回方剂中对应的靶点,以entrezid为返回值。
    targets = herb_mol_target[herb_mol_target['herb_cn_name'].isin(herbs)]
    herb_targets = targets['TARGET_ID']
    pd_gene_symbol = di.gene_symbol_entrezid_pd()
    pd_gene_symbol = pd_gene_symbol[pd_gene_symbol['target'].isin(herb_targets)]
    return list(pd_gene_symbol['entrezid'])


def ChinesesNameToEnglish():
    filepath = 'D:\\ctm_data\\'
    filename = 'chinese_english_herb.csv'
    pd_chinese_english = pd.read_csv(filepath + filename,sep = '\t')
    chinese_english_herb_dict = {key:values for key, values in zip(pd_chinese_english['herb_cn_name'], pd_chinese_english['herb_en_name'])}
    return chinese_english_herb_dict


if  __name__ == '__main__':

    gene_entrezid_dict = di.gene_ENTREZID_dict()
    G = di.graphFromPPI()
    filepath_beijing = 'D:\\张思琪\\背景网络整理0728\\背景网络整理0728\\'
    filename_beijing = '背景网络gene.csv'
    beijing_gene_pd = pd.read_csv(filepath_beijing+filename_beijing)
    beijing_gene_list = list(beijing_gene_pd['gene'])
    print(beijing_gene_list)

with open('edges.csv','a') as f:
    #print('from\t','to\t','type\t','weight')
    for (node1, node2) in G.edges():
        if  node1 in gene_entrezid_dict and node2 in gene_entrezid_dict and node1!=node2:
            #if gene_entrezid_dict[node1] in beijing_gene_list and gene_entrezid_dict[node2] in beijing_gene_list:
                f.write(str(gene_entrezid_dict[node1]))
                f.write('\t')
                f.write(str(gene_entrezid_dict[node2]))
                #f.write('\t')
                #f.write('hyperlink')
                #f.write('\t')
                #f.write('1')
                #f.write('\t')
                f.write('\n')
                f.flush()

