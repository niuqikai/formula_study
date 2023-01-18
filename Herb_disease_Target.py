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
    filepath = 'data/'
    filename = 'chinese_english_herb.csv'
    pd_chinese_english = pd.read_csv(filepath + filename,sep = '\t')
    chinese_english_herb_dict = {key:values for key, values in zip(pd_chinese_english['herb_cn_name'], pd_chinese_english['herb_en_name'])}
    return chinese_english_herb_dict


if  __name__ == '__main__':
    chinese_english_herb_dict = ChinesesNameToEnglish()

    filepath = 'data/'
    filename = 'TCMSP_DB_加工.xlsx'

    herb_mols =  di.herb_molecules(filepath, filename) #中药对应的成分
    targets_mol_herb = di.targets_mol_herb(filepath,filename)#中药靶点成分 inner
    herb_mol_target = di.herb_mol_targets(filepath, filename) # 计算每种中药对应的成分和靶点


    #filepath_herbpair = 'D:\\ctm_data\\'
    #filename = '叶天士新.csv'
    #filename = '第一批100首-药物组成.csv'
    #filename_herbpair = '中成药数据库.csv'

    #herb_pair_from_data = hpff.herb_pair_score_from_data(filepath_herbpair,filename_herbpair,herb_mols)
    #formula_score_dict = {}

    #si1 疾病和疾病对应的靶点
    #disease_targ_name = di.disease_targetname(filepath, filename)
    #disease_targ_name.to_csv('astha_name.csv')
    #print(disease_targ_name)

    #1 中药对应的疾病靶点数量、成分数量和得分
    #herb_mols_targs_score(targ_mol_herb, herb_score)

    #2药物配伍得分和jacard分数等
    #fhi.herb_herb_jaccard_gini(herb_mol_target)

    #随机生成和组方数量相同的方剂
    #herbs = list(targets_mol_herb['herb_cn_name'])
    #gene_herbnum_formula(herbs,herb_score_dict,herb_pair_from_data)

    all_target = []
    disease_tar = di.disease_target(filepath,filename)
    pd_disease_symbol = di.gene_symbol_entrezid_pd()
    pd_disease_symbol = pd_disease_symbol[pd_disease_symbol['target'].isin(disease_tar)]

    #herbs = ['丹参','川芎','延胡索','杜仲','甘草','白术','苦杏仁','葛根','麻黄']
    #herbs = ['丹参','党参','山茱萸','枸杞子','甘草']
    herbs = ['丹参','川芎','延胡索','当归','没药','甘草','黄芩']
    entrezid_disease_target = list(set(pd_disease_symbol['entrezid']))
    rs_formula_target = list(set(formula_tareget(herb_mol_target, herbs)))# 返回方剂中对应的靶点,以
    gene_entrezid_dict = di.gene_ENTREZID_dict()

    G = di.graphFromPPI()
    #输入方剂，返回对应的靶点


    print('Node_id\t','type\t','typelabel\t','size')
    for rs in rs_formula_target:
        if str(rs) in G.nodes() and int(rs) not in entrezid_disease_target:
            if str(rs) in gene_entrezid_dict:
                print(rs,end='\t')
                print(gene_entrezid_dict[str(rs)],end='\t')
                print('Herb_Target',end='\t')
                print('1',end='\t')
                print(str(G.degree(str(rs))),end='\t')
                print()



    for i in entrezid_disease_target:
        if str(i) in G.nodes() and int(i) not in rs_formula_target:
            if str(i) in gene_entrezid_dict:
                print(i, end='\t')
                print(gene_entrezid_dict[str(i)], end='\t')
                print('Disease_Target', end='\t')
                print('2', end='\t')
                print(str(G.degree(str(i))),end='\t')
                print()

    for j in entrezid_disease_target:
        if str(j) in G.nodes() and int(j) in rs_formula_target:
            if str(j) in gene_entrezid_dict:
                print(j, end='\t')
                print(gene_entrezid_dict[str(j)], end='\t')
                print('Together_Target', end='\t')
                print('3', end='\t')
                print(str(G.degree(str(j))),end='\t')
                print()


    all_target.extend(rs_formula_target)
    all_target.extend(entrezid_disease_target)
    all_target = list(set(all_target))
    nodes_no_target = []

    '''
    for (node1,node2) in G.edges():
        if int(node1) in all_target and int(node2) not in all_target:
            if node2 in gene_entrezid_dict:
                print(node2, end='\t')
                print(gene_entrezid_dict[node2], end='\t')
                print('Common_Node', end='\t')
                print('4', end='\t')
                print(str(G.degree(str(node2))), end='\t')
                print()
        if int(node2) in all_target and int(node1) not in all_target:
            if node1 in gene_entrezid_dict:
                print(node1, end='\t')
                print(gene_entrezid_dict[node1], end='\t')
                print('Common_Node', end='\t')
                print('4', end='\t')
                print(str(G.degree(str(node1))), end='\t')
                print()
    '''


with open('edges.csv','a') as f:
    #print('from\t','to\t','type\t','weight')
    for (node1, node2) in G.edges():
        if int(node1) in all_target and int(node2)  in all_target and node1 in gene_entrezid_dict and node2 in gene_entrezid_dict and node1!=node2:
            f.write(str(gene_entrezid_dict[node1]))
            f.write('\t')
            f.write(str(gene_entrezid_dict[node2]))
            f.write('\t')
            f.write('hyperlink')
            f.write('\t')
            f.write('1')
            f.write('\t')
            f.write('\n')
            f.flush()

