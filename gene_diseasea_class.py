from sklearn.model_selection import train_test_split
import pandas as pd
import xgboost as xgb
from sklearn.metrics import roc_auc_score
# import matplotlib.pyplot as plt
import matplotlib.pylab as plt
import seaborn as sns
import formula_study.formula_herb_importance as fhi
# data = data.fillna(-1)
filepath = 'D:\\libing\\基因分类\\'
filelabel = 'trans.csv'
df = pd.read_csv(filepath + filelabel)

filept = 'D:\\libing\\基因分类\\GSE62254-expression-good_new.csv'
test_x = pd.read_csv(filept)
print(test_x.columns)

filelabel = 'train_label.csv'
df_label = pd.read_csv(filepath + filelabel)
df_label_dict = {key:values for key, values in zip(df_label['names'], df_label['label'])}
print(df_label_dict)
gene_list = list(df['gene'])
print(gene_list)
new_label_list = []
for key in gene_list:
    new_label_list.append(int(df_label_dict[key]))

'''
#xfeaturelist = []
xfeaturelist = list(df.columns)
xfeaturelist = xfeaturelist[1:]
newxfeaturelist = []
for xi in xfeaturelist:
    if xi != '15-Sep' and xi in list(test_x.columns):
        newxfeaturelist.append(str(xi).strip())
xfeaturelist = newxfeaturelist
print(xfeaturelist)
'''
xfeaturelist = [
'CAMKV',
'F10',
'C15orf59',
'QPRT',
'MAGEA12',
'GLI3',
'MAGEA6',
'LINC00958',
'CCNE1',
'NRSN2',

#'KMT2D',
#'ANK3',
#'RYR1',
#'KMT2C',
#'MYCBP2',
#'MDN1',
#'TTN',
'DNAH3',
#'HUWE1',
'XYLT2',
#'MYO15A',
'PIK3CA',
#'COL7A1',
'KALRN',
'RGS12',
'NEB',
#'VPS13B',
'CREBBP',
#'HIVEP3',
#'ARID1A',
#'TP53',
#'CDH1',
#'NRXN1',
#'CNBD1',
'DOCK4',
#'KCNT2',
#'OR4N2',
#'APC',
#'DNAH7',
#'LYST',
#'NLRP8',
#'OR8J3',
#'PIK3CG',
#'RGPD4',
#'USH2A',
#'NPAP1',
#'AFF3',
#'BBS9',
#'FBN2',
'ZBTB20',
'TG',
#'WDFY3',
#'NF1',
'SCAF4',
#'PLEC',
'DNAH2',
#'RNF43',
#'NOTCH4',
'CELSR1',
'FAT1'
#'BDP1',
#'BTBD11'
]

df_x = df
x = df_x[xfeaturelist].values
#y = df_x['label'].values
y = new_label_list
# X_train, X_te, y_train, y_te = train_test_split(x, y, test_size = 0.00001, random_state=8) # 为了看模型在没有见过数据集上的表现，随机拿出数据集中30%的部分做测试
X_train = x
y_train = y

X_test = x
y_test = y  #

dtrain = xgb.DMatrix(data=X_train, label=y_train)
dtest = xgb.DMatrix(data=X_test, label=y_test)
param = {'max_depth': 7, 'eta': 0.075, 'silent': 1, 'lambda': 4, 'gamma': 0.05, 'subsample': 0.7,
         'consample_bytree': 0.8,
         'objective':'multi:softmax','nthread': 4, 'min_child_weight': 1, 'seed': 3,'num_class':3}
# evallist  = [(dtest,'eval'), (dtrain,'train')]
evallist = [(dtrain, 'train')]

num_round = 73
bst = xgb.train(param, dtrain, num_round, evallist)
fnames = bst.feature_names
bst.feature_names = xfeaturelist
# bst.feature_names = ['user_gender','user_age','user_marriage','user_education','academic_stage','risk_chnl_branch_org_level','phone_status']

'''
feature_score = bst.get_fscore()
feat_imp = pd.Series(bst.get_fscore()).sort_values(ascending=False)
feat_imp.plot(kind='bar', title='Feature Importances')
plt.ylabel('Feature Importance Score')
plt.show()
'''

from xgboost import plot_importance
plt.figure(figsize=(12, 6))
plot_importance(bst, max_num_features=20)
plt.title("Featurertances")
plt.show()

bst.feature_names = fnames

def writeresulttofile(writefile,gene_list,new_label_list,rs_train):
    with open(writefile,'a') as wf:
        for i in range(len(gene_list)):
            wf.write(str(gene_list[i]))
            wf.write(',')
            if new_label_list !=0:
                wf.write(str(new_label_list[i]))
                wf.write(',')
            wf.write(str(int(rs_train[i])))
            wf.write('\n')

y_pred_train = bst.predict(dtrain)
rs_train = y_pred_train
#writeresulttofile('train_rs.csv',gene_list,new_label_list,rs_train)


'''
filept = 'D:\\libing\\基因分类\\GSE62254-expression-good.csv'
test_x = pd.read_csv(filept)
#转置

test_x = pd.DataFrame(test_x.values.T,index= test_x.columns,columns=test_x.index)
print(test_x)
test_x.to_csv('GSE62254-expression-good_new.csv')
'''

test_x_value = test_x[xfeaturelist].values
test_x_value = xgb.DMatrix(data=test_x_value)
rs_test = bst.predict(test_x_value)
print(rs_test)
writeresulttofile('test_rs.csv',list(test_x['ID_REF']),0,rs_test)

#fhi.writelisttodata('test.csv',rs_test)

