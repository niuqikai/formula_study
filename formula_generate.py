#根据算法生成组方和组方分数评价
import formula_study.formula_herb_importance as fhi
import formula_study.Data_input as di
import formula_study.herb_pairs_from_formula as hpff
import random as rd

#初始化方剂
def innit_formula_seed(herb_score_dict, row_num , col_num_list):#herb_score_dict是中药分数字典,row_num为要生成的方剂数目,col_num为方剂中的中药数
    formula_herb_list = [] #方剂列表，包含有可能的row_num个方剂组成
    herb_score_dict = sorted(herb_score_dict.items(),key = lambda x: x[1], reverse=False)

    for i in range(row_num):
        formula_herb_list_seed = []  # 方剂种子，初始生成的方剂
        while(len(formula_herb_list_seed) < col_num_list[i]):
            benchmark_num = 0
            r = rd.random()
            for (k,v) in herb_score_dict:
                if r < benchmark_num + float(v) and k not in formula_herb_list_seed:
                    formula_herb_list_seed.append(k)
                    break
                benchmark_num = benchmark_num + float(v)
        #print(formula_herb_list_seed)
        formula_herb_list_seed = sorted(formula_herb_list_seed)
        formula_herb_list.append(formula_herb_list_seed)
    return formula_herb_list

#计算方剂得分
def compute_formula_score(formula_list , herb_score_dict, pair_num_dict):
    formula_score = 0.0
    for herb_i in range(len(formula_list)-1):
        for herb_j in range(herb_i + 1 , len(formula_list)):
            if (str(formula_list[herb_i])+str(formula_list[herb_j])) in pair_num_dict or (str(formula_list[herb_j])+str(formula_list[herb_i])) in pair_num_dict:
                if str(formula_list[herb_i])+str(formula_list[herb_j]) in pair_num_dict:
                    formula_score += herb_score_dict[formula_list[herb_i]] * herb_score_dict[formula_list[herb_j]] * pair_num_dict[str(formula_list[herb_i])+str(formula_list[herb_j])]
                elif str(formula_list[herb_j])+str(formula_list[herb_i]) in pair_num_dict:
                    formula_score += herb_score_dict[formula_list[herb_i]] * herb_score_dict[formula_list[herb_j]] * pair_num_dict[str(formula_list[herb_j])+str(formula_list[herb_i])]
            else:
                formula_score = formula_score #- herb_score_dict[formula_list[herb_i]] * herb_score_dict[formula_list[herb_j]]
    return formula_score/len(formula_list)#(len(formula_list)*len(formula_list))

def is_same_herb_in_formula(formulalist):#判断方剂里面是否有重复的中药
    if len(formulalist) == len(set(formulalist)):
        return True
    else:
        return False


def variation(herbs_list,p,herb_score_dict):#变异
    herb_score_dict = sorted(herb_score_dict.items(),key = lambda x: x[1], reverse=False)
    herbs_list_new = []
    for herbi in herbs_list:
        r = rd.random()
        if r > p:
            herbs_list_new.append(herbi)
        else:
            benchmark_num = 0.0
            d = rd.random()
            for (k,v) in herb_score_dict:
                if  d < benchmark_num + float(v): #and k not in herbs_list_new:
                    herbs_list_new.append(k)
                    break
                benchmark_num = benchmark_num + float(v)

    return herbs_list_new

def generate_formula_list(rows_num):#生成1-15的方剂列表
    rows_num_list = []
    for i in range(rows_num):
        num = rd.randint(1,15)#从2味药到15味药
        rows_num_list.append(num)
    return rows_num_list

#根据方剂得分,删除后一半得分较低的方剂，然后，交叉、变异，生成新的方剂组
def Genetic_Algorithm(formulas_score_dict,herb_score_dict,herb_pair_from_data,formula_nums):
    formulas_score_dict = sorted(formulas_score_dict.items(),key = lambda x: x[0], reverse=True)#排序
    print(len(formulas_score_dict),formulas_score_dict)
    #保留一半得分高的方剂
    new_formulas_score_list1 = []
    new_formulas_score_list2 = []

    for i in range(int(len(formulas_score_dict)/2)):
        new_formulas_score_list1.append(formulas_score_dict[i])
    formulas_num = len(new_formulas_score_list1)
    rd.shuffle(new_formulas_score_list1)
    for i in range(0,formulas_num,2):
        (scorei,formulai) = new_formulas_score_list1[i]
        (scorej,formulaj) = new_formulas_score_list1[i+1]

        formulak = []
        formulak.extend(formulai)
        formulak.extend(formulaj)
        rd.shuffle(formulak)
        #组成新的方子,交叉变异
        formula1 = formulak[0:int(len(formulak)/2)]
        formula1 = variation(formula1,0.05,herb_score_dict)#变异，变异系数0.05
        formula1 = sorted(formula1)#相同的方剂序列一致
        formula1_list = (compute_formula_score(formula1, herb_score_dict, herb_pair_from_data)* 10000,formula1)

        formula2 = formulak[int(len(formulak)/2):len(formulak)]
        formula2 = variation(formula2,0.05,herb_score_dict)#变异，变异系数0.05
        formula2 = sorted(formula2)
        formula2_list = (compute_formula_score(formula2, herb_score_dict, herb_pair_from_data) *10000,formula2)
        if is_same_herb_in_formula(formula1):
            new_formulas_score_list2.append(formula1_list)
        if is_same_herb_in_formula(formula2):
            new_formulas_score_list2.append(formula2_list)
    new_formulas_score_list1.extend(new_formulas_score_list2)
    new_formulas_score_dict = {key:value for (key,value) in new_formulas_score_list1}
    while (len(new_formulas_score_dict.keys()) != formula_nums):#有时候会有相同的方剂，使得总体数量减少，重新生成，补齐
        f_h_l = innit_formula_seed(herb_score_dict, 1, [rd.randint(1,15)])
        f_score = compute_formula_score(f_h_l[0], herb_score_dict, herb_pair_from_data) * 10000
        new_formulas_score_dict[f_score] = f_h_l[0]
    return  new_formulas_score_dict


if __name__ == '__main__':
    filepath = 'D:\\ctm_data\\TCMSP-数据\\'
    filename = 'TCMSP_DB_加工.xlsx'

    target_molecule = di.target_mol(filepath, filename, tar='0')#获取疾病对应的靶点和成分，tar=0 表示制定疾病的靶点和成分
    herb_mols =  di.herb_molecules(filepath, filename) #中药对应的成分
    herb_score = fhi.herb_walk_score_interation(target_molecule,herb_mols)#计算对应的药物分数
    herb_score = herb_score[['herb_cn_name','walk_score']].drop_duplicates()

    herb_score['walk_score'] = herb_score.apply(lambda x: x['walk_score']/herb_score['walk_score'].sum(),axis=1)#归一化
    herb_score_dict = {key:values for key, values in zip(herb_score['herb_cn_name'], herb_score['walk_score'])}#转换为字典结构

    '''
    pair_score = 'herb_herb_mol_jaccard_gini'
    pair_s = di.data_from_excel_sheet(filepath + filename, pair_score)
    p_s = pair_s[['herb1','herb2','cos_mol']]
    p_s_dict = {(key1 ,key2):values for key1, key2 ,values in zip(p_s['herb1'], p_s['herb2'], p_s['cos_mol'])}#转换为字典结构
    '''

    formula_nums = 1000
    filepath = 'D:\\ctm_data\\'
    #filename = '叶天士新.csv'
    #filename = '第一批100首-药物组成.csv'
    filename = '中成药数据库.csv'

    rows_list = generate_formula_list(formula_nums)#随机生成方剂中中药数目
    herb_pair_from_data = hpff.herb_pair_score_from_data(filepath,filename,herb_mols)
    formula_score_dict = {}
    formulas = innit_formula_seed(herb_score_dict, formula_nums, rows_list)
    #test_formula = ['茯苓','人参','桂枝','当归']
    #print(compute_formula_score(test_formula, herb_score_dict, herb_pair_from_data)*10000)

    for formula_list in formulas:
        formula_score_dict[compute_formula_score(formula_list, herb_score_dict, herb_pair_from_data)*10000] = formula_list
    while(len(formula_score_dict.keys())!=formula_nums):
        f_h_l = innit_formula_seed(herb_score_dict, 1, [rd.randint(1, 15)])
        f_score = compute_formula_score(f_h_l[0], herb_score_dict, herb_pair_from_data) * 10000
        formula_score_dict[f_score] = f_h_l[0]

    max_score = -99999
    for k in range(1000):
        print("第"+str(k)+"迭代")
        new_score_dict = Genetic_Algorithm(formula_score_dict,herb_score_dict,herb_pair_from_data,formula_nums)
        formula_score_dict = new_score_dict

