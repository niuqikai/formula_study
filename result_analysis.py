import formula_herb_importance as fhi
import Data_input as di
import Data_output as do
import herb_pairs_from_formula as hpff
import random as rd
import formula_generate as fg
import pandas as pd
import PPI_analyse as ppi
#分析生成的方剂和各种中药分析
def Sab_pairscore_compare(writefilename,herb_pair_from_data,herbs_Sab):#比较药对得分和Sab分数
    return 0

def write_herb_pair_score(filename,herb_pair_from_data):#将配伍得分写入文件
    do.writedicttodata(filename , herb_pair_from_data)


def write_to_disease_target():
    print()

def herb_mols_targs_score(targ_mol_herb,herb_score):#针对特定疾病的中药对应的靶点、成分和得分
    targ_mol,targ_herb = fhi.targ_mol_herb_num(targ_mol_herb)
    mols_tars_num = pd.merge(targ_mol,targ_herb,how = 'left',on= 'herb_cn_name').reset_index()

    mols_tars_score_num = pd.merge(mols_tars_num,herb_score,how = 'left',on= 'herb_cn_name').reset_index()
    mols_tars_score_num.to_csv('heartdisease_mols_tars_score_num.csv')

def gene_herbnum_formula(herbs, herb_score_dict, pair_num_dict):#生成固定组方数目的随机方剂
    formula_num = 1000
    formula_num_dict = {}
    for i in range(formula_num):
        num = rd.randint(2, 15)  # 从2味药到14味药
        rd.shuffle(herbs)
        formula_list = herbs[0:num]
        print(formula_list)
        score = fg.compute_formula_score(formula_list , herb_score_dict, pair_num_dict)*10000
        formula_num_dict[score] = formula_list
    do.writeformulatodata('lc_rand_formula.csv',formula_num_dict)
    #return formula_num_dict

def ChinesesNameToEnglish():
    filepath = 'data/'
    filename = 'chinese_english_herb.csv'
    pd_chinese_english = pd.read_csv(filepath + filename,sep = '\t')
    chinese_english_herb_dict = {key:values for key, values in zip(pd_chinese_english['herb_cn_name'], pd_chinese_english['herb_en_name'])}
    return chinese_english_herb_dict


def formula_num_herbs_num(chinese_english_herb_dict):#计算Top50的方剂中药物分布，药对分布，方剂数量分布
    #filepath = 'D:\\formula_result\\1算法设计\\5方剂分析\\'
    #filename = 'atherosclerosis.csv'
    filepath = 'D:\\network_ctm\\formula_study\\'
    filename = 'Atherosclerosis_pagerankob30.csv'

    herbs_dict = {} #方剂中中药数量分布
    formula_num_dict = {}#方剂数量分布
    pair_num_dict = {}

    with open(filepath + filename) as f:
        for l in f:
            formula_num = 0
            herbs = l.strip().split(',')
            for herb in herbs:
                if herb!='' :
                    formula_num += 1
                    if herb in herbs_dict:
                        herbs_dict[herb] = herbs_dict[herb] + 1
                    if herb not in herbs_dict:
                        herbs_dict[herb] = 1

                for herbj in herbs:
                    if herb!=herbj and herb!='' and herbj!='':
                        if (herb,herbj) in pair_num_dict:
                            pair_num_dict[herb,herbj] = pair_num_dict[herb,herbj] + 1
                        else:
                            pair_num_dict[herb,herbj] =  1


            if formula_num in formula_num_dict:
                formula_num_dict[formula_num] = formula_num_dict[formula_num] + 1
            else:
                formula_num_dict[formula_num] = 1
    herbs_dict = sorted(herbs_dict.items(),key= lambda x:x[1], reverse=True)
    print(herbs_dict)
    for (k,v) in herbs_dict:
        print(str(chinese_english_herb_dict[k]).strip(),',',v)
    for (k,v) in formula_num_dict.items():
        print(k,v)

    Herbs = ['Herb1','Herb2','Herb3','Herb4','Herb5','Herb6','Herb7','Herb8','Herb9','Herb10']
    for k in range(10):
            (herbk,v) = herbs_dict[k]
            for j in range(10):
                (herbj,t) = herbs_dict[j]
                if k!=j :
                    p = (str(herbk),str(herbj))
                    if p in pair_num_dict:
                        print(Herbs[k],'\t',Herbs[j],'\t',pair_num_dict[p])
                    else:
                        print(Herbs[k],'\t',Herbs[j],'\t',0)
                else:
                    print(Herbs[k],'\t',Herbs[j],'\t',0)


    '''
    print('Herb',end='\t')
    for (k,v) in herbs_dict[:10]:
        print(chinese_english_herb_dict[k],end='\t')
    print()
    for (k,v) in herbs_dict[:10]:
        print(chinese_english_herb_dict[k],end='\t')
        for (j,t) in herbs_dict[:10]:
            if k != j:
                p = (str(k), str(j))
                if p in pair_num_dict:
                    print(pair_num_dict[p],end='\t')
                else:
                    print(0,end='\t')
            else:
                print(0,end='\t')
        print()
    '''

def formula_tareget(herb_mol_target,herbs):#返回方剂中对应的靶点
    targets = herb_mol_target[herb_mol_target['herb_cn_name'].isin(herbs)]
    herb_targets = targets['TARGET_ID']
    pd_gene_symbol = di.targetid_SYMBOL_pd()
    pd_gene_symbol = pd_gene_symbol[pd_gene_symbol['target'].isin(herb_targets)]
    return list(pd_gene_symbol['symbol'])


def diseaseIDinPPI(filepath,filename):
    disease_tar = di.disease_target(filepath,filename)
    gene_symbol_entrezid = di.gene_symbol_entrezid_pd()
    target_entrezid = gene_symbol_entrezid['entrezid']
    return target_entrezid

if  __name__ == '__main__':
    chinese_english_herb_dict = ChinesesNameToEnglish()
    formula_num_herbs_num(chinese_english_herb_dict)

    filepath = 'data/'
    filename = 'TCMSP_DB_加工.xlsx'

    target_molecule = di.target_mol(filepath, filename, tar='0')#获取疾病对应的靶点和成分，tar=0 表示制定疾病的靶点和成分
    herb_mols =  di.herb_molecules(filepath, filename) #中药对应的成分
    targets_mol_herb = di.targets_mol_herb(filepath,filename)#中药靶点成分 inner
    degree,pagerank,eigenvector,closeness,betweenness = ppi.symbol_sore_from_PPI()

    herb_score = fhi.herb_walk_score_interation(targets_mol_herb,pagerank)#计算对应的药物分数
    herb_score = herb_score[['herb_cn_name','walk_score']].drop_duplicates()

    targ_mol_herb = di.targets_mol_herb(filepath,filename)#疾病靶点对应的成分和中药
    herb_mol_target = di.herb_mol_targets(filepath, filename) # 计算每种中药对应的成分和靶点
    hmtd = di.herb_mol_targets_disease(filepath,filename)

    herb_score['walk_score'] = herb_score.apply(lambda x: x['walk_score']/herb_score['walk_score'].sum(),axis=1)#归一化
    herb_score_dict = {key:values for key, values in zip(herb_score['herb_cn_name'], herb_score['walk_score'])}#转换为字典结构

    '''
    pair_score = 'herb_herb_mol_jaccard_gini'
    pair_s = di.data_from_excel_sheet(filepath + filename, pair_score)
    p_s = pair_s[['herb1','herb2','cos_mol']]
    p_s_dict = {(key1 ,key2):values for key1, key2 ,values in zip(p_s['herb1'], p_s['herb2'], p_s['cos_mol'])}#转换为字典结构
    '''

    filepath_herbpair = 'data/'
    #filename = '叶天士新.csv'
    #filename = '第一批100首-药物组成.csv'
    filename_herbpair = '伤寒金匮.csv'

    herb_pair_from_data = hpff.herb_pair_score_from_data(filepath_herbpair,filename_herbpair,herb_mols)
    formula_score_dict = {}

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

    #输入方剂，返回对应的靶点
    #herbs = ['丹参','党参','山茱萸','枸杞子','甘草']
    #herbs = ['丹参','川芎','延胡索','当归','没药','甘草','黄芩']
    herbs = ['丹参','川芎','延胡索','杜仲','甘草','白术','苦杏仁','葛根','麻黄']
    rs_formula_target = formula_tareget(herb_mol_target, herbs) # 返回方剂中对应的靶点,以
    for rs in rs_formula_target:
        print(rs)

    #print(diseaseIDinPPI(filepath,filename))


    '''
    for formula_list in formulas:
        formula_score_dict[fg.compute_formula_score(formula_list, herb_score_dict, herb_pair_from_data)*10000] = formula_list
    while(len(formula_score_dict.keys())!=formula_nums):
        f_h_l = fg.innit_formula_seed(herb_score_dict, 1, [rd.randint(1, 15)])
        f_score = fg.compute_formula_score(f_h_l[0], herb_score_dict, herb_pair_from_data) * 10000
        formula_score_dict[f_score] = f_h_l[0]
    
    max_score = -99999
    for k in range(1000):
        new_score_dict = fg.Genetic_Algorithm(formula_score_dict,herb_score_dict,herb_pair_from_data,formula_nums)
        formula_score_dict = new_score_dict
    '''
    #write_herb_pair_score('叶天士_pair_score.csv',herb_pair_from_data)


