import pandas as pd
'''
filepath = 'C:\\Users\\Administrator\\Desktop\\'
filename = '工作簿3.xlsx'
sheetname = 'Sheet1'

data1 = pd.read_excel(filepath + filename,sheetname)
data2 = pd.read_excel(filepath + filename,sheetname)
data = pd.merge(data1, data2 ,how = 'left', on = '古代经典名方基本信息')
data_write = data.groupby(['填报机构_x','填报机构_y'])['古代经典名方基本信息'].unique()
data_write.to_csv('fomula2.csv')

print(data_write)
'''
import os
def file_name(file_dir):
    for root, dirs, files in os.walk(file_dir):
        #print(root) #当前目录路径
        #print(dirs) #当前路径下所有子目录
        if len(files) != 0:
            filename = str(root) + '\\' + str(files[0])
            df = pd.read_excel(filename)
            print(df)
            num_index = int(df['序号'].max())#2
            num_count = int(df['Unnamed: 22'].count())#13

            pivot_num = int((num_count - 1)/num_index) #6
            #print(pd.pivot(df, index= '序号',columns=['Unnamed: 17','Unnamed: 18','Unnamed: 19','Unnamed: 20','Unnamed: 21','Unnamed: 22']))
            data_list_list = []
            for i in range(num_index):
                data_list = []
                data_list.append(df.iloc[1+ i*pivot_num,0:16])
                for j in range(0,pivot_num):
                    data_list.append(df.iloc[[2+ i*num_index + pivot_num],[17,18,19,20,21,22]])
                data_list_list.append(data_list)
            print(data_list_list)
file_dir = 'data\classical_formula'
file_name(file_dir)

