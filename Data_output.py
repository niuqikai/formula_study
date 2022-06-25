
#根据中药名称清洗数据
def standedherb(herb,herb_mols):
    for herbi in herb_mols:
        if str(herb) in str(herbi) or str(herbi) in str(herb):
            return herbi
    return herb

def writestandardformula(filename,herb_mols):
    with open(filename, encoding='utf-8')  as fl:
        for line in fl:
            herbs_list_new = []
            herbs_list = list(herb_mols['herb_cn_name'].unique())
            herbs = str(line).strip().split(',')
            for i in herbs:
                if i != '':
                    herbs_list_new.append(str(standedherb(i,herbs_list)))
            writeformulatodata('叶天士新.csv',herbs_list_new)

'''
def writeformulatodata(filename , datalist):#将重新写新的方剂
    with open(filename,'a') as fl:
        for dl in range(len(datalist)):
            fl.write(str(datalist[dl]))
            fl.write(',')
        fl.write('\n')
'''
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

def writeformulatodata(filename,datadict):#
    with open(filename,'a') as fl:
        for dl in datadict.keys():
            try:
                fl.write(str(dl))
                fl.write(",")
                for value in datadict[dl]:
                    #print(value)
                    fl.write(str(value))
                    fl.write(",")
                fl.write('\n')
            except:
                pass

def writedatalisttodata(filename,datalist):
    with open(filename,'a') as fl:
        for dl in datalist:
            fl.write(str(dl))
            fl.write(",")
        fl.write('\n')
