import pandas as pd
#import formula_study.formula_herb_importance as fhi
#根据方剂生成中药对
filepath = 'D:\\ctm_data\\'
filename = '中成药数据库.csv'
def herb_pair_score_from_data(filepath, filename):
    herb_pair = {}
    with open(filepath + filename,encoding='utf-8')  as fl:
        for line in fl:
            herbs_list = []
            herbs = str(line).strip().split(',')
            for i in herbs:
                if i != '':
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


def writelisttodata(filename , datalist):#将列表数据写入文本
    with open(filename,'a') as fl:
        for dl in range(len(datalist)):
            fl.write(str(datalist[dl]))
            fl.write('\n')

def writedicttodata(filename , datadict):#将字典数据写入文本
    with open(filename,'a') as fl:
        for dl in datadict.keys():
            try:
                fl.write(str(dl))
                fl.write(",")
                fl.write(str(datadict[dl]))
                fl.write('\n')
            except:
                pass
herb_pair = herb_pair_score_from_data(filepath, filename)
#writedicttodata('herb_pair_from_formula2.csv',herb_pair)



