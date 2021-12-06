import pandas as pd
filepath = 'C:\\Users\\Administrator\\Desktop\\'
filename = '工作簿3.xlsx'
sheetname = 'Sheet1'

data1 = pd.read_excel(filepath + filename,sheetname)
data2 = pd.read_excel(filepath + filename,sheetname)
data = pd.merge(data1, data2 ,how = 'left', on = '古代经典名方基本信息')
data_write = data.groupby(['填报机构_x','填报机构_y'])['古代经典名方基本信息'].unique()
data_write.to_csv('fomula2.csv')

print(data_write)