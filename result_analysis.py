import formula_study.formula_herb_importance as fhi
import formula_study.Data_input as di
import formula_study.Data_output as do
import formula_study.herb_pairs_from_formula as hpff
import random as rd
import formula_study.formula_generate as fg
import pandas as pd

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

if  __name__ == '__main__':
    filepath = 'D:\\ctm_data\\TCMSP-数据\\'
    filename = 'TCMSP_DB_加工.xlsx'

    target_molecule = di.target_mol(filepath, filename, tar='0')#获取疾病对应的靶点和成分，tar=0 表示制定疾病的靶点和成分
    herb_mols =  di.herb_molecules(filepath, filename) #中药对应的成分
    targets_mol_herb = di.targets_mol_herb(filepath,filename)#中药靶点成分 inner

    herb_score = fhi.herb_walk_score_interation(targets_mol_herb)#计算对应的药物分数
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

    formula_nums = 1000
    filepath_herbpair = 'D:\\ctm_data\\'
    #filename = '叶天士新.csv'
    #filename = '第一批100首-药物组成.csv'
    filename_herbpair = '中成药数据库.csv'

    rows_list = fg.generate_formula_list(formula_nums)#随机生成方剂中中药数目
    herb_pair_from_data = hpff.herb_pair_score_from_data(filepath_herbpair,filename_herbpair,herb_mols)
    formula_score_dict = {}

    #si1 疾病和疾病对应的靶点
    #disease_targ_name = di.disease_targetname(filepath, filename)
    #disease_targ_name.to_csv('astha_name.csv')
    #print(disease_targ_name)

    #1 中药对应的疾病靶点数量、成分数量和得分
    #herb_mols_targs_score(targ_mol_herb, herb_score)

    #2药物配伍得分和jacard分数等
    fhi.herb_herb_jaccard_gini(herb_mol_target)

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


