import pandas as pd
#import formula_study.formula_herb_importance as fhi
#根据方剂生成中药对
import Data_input as di
import Data_output as do


def herb_pair_score_from_data(filepath, filename,herb_mols):
    herb_pair = {}
    #herbs_list = list(herb_mols['herb_cn_name'].unique())
    #print(herbs_list)
    with open(filepath + filename,encoding='utf-8')  as fl:
        for line in fl:
            herbs_list = []
            herbs = str(line).strip().split(',')
            for i in herbs:
                if i != '':
                    #herbs_list.append(str(standedherb(i,herbs_list)))
                    herbs_list.append(str(i))
            for j in range(0,len(herbs_list) - 1):
                for k in range(j + 1, len(herbs_list)):
                    if (str(herbs_list[j])+str(herbs_list[k])) not in herb_pair:
                        herb_pair[str(herbs_list[j])+str(herbs_list[k])] = 1
                    else:
                        herb_pair[str(herbs_list[j])+str(herbs_list[k])] += 1
                    if (str(herbs_list[k])+str(herbs_list[j])) not in herb_pair:
                        herb_pair[str(herbs_list[k])+str(herbs_list[j])] = 1
                    else:
                        herb_pair[str(herbs_list[k])+str(herbs_list[j])] += 1
    return herb_pair


def herb_pair_score_from_Sab(filepath, filename,herb_mols):
    herb_pair = {}
    st_name = 'Sab和pairscore关系'
    #herbs_list = list(herb_mols['herb_cn_name'].unique())
    #print(herbs_list)
    df = pd.read_excel(filepath + filename, sheet_name = st_name)  # 可以通过sheet_name来指定读取的表单
    herb_pair = {key:values for key, values in zip(df['herb1herb2'], df['sab'])}#转换为字典结构
    return herb_pair

def herb_pair_score_from_shangshi(filepath, filename,herb_mols):
    herb_pair = {}
    st_name = 'Sab和pairscore关系'
    #herbs_list = list(herb_mols['herb_cn_name'].unique())
    #print(herbs_list)
    df = pd.read_excel(filepath + filename, sheet_name = st_name)  # 可以通过sheet_name来指定读取的表单
    herb_pair = {key:values for key, values in zip(df['herb1herb2'], df['上市中成药'])}#转换为字典结构
    return herb_pair

#herb_pair = herb_pair_score_from_data(filepath, filename)
#writedicttodata('herb_pair_from_formula2.csv',herb_pair)
#filepath = 'D:\\ctm_data\\TCMSP-数据\\'
#filename = 'TCMSP_DB_加工.xlsx'
#herb_mols =  di.herb_molecules(filepath, filename) #中药对应的成分

#filepath = 'D:\\ctm_data\\'
#filename = '叶天士.csv'
#do.writestandardformula(filepath + filename,herb_mols)


