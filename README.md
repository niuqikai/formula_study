#代码使用说明#<br>

##1 如何根据疾病筛选方剂##<br>
formula_generate.py为入口程序，生成的结果保留在当前文件夹下，结果会保留排序前50方剂的方剂和随机生成的5个方剂。<br>
'diseasename/XX.csv'为筛选疾病的对应的靶点如'diseasename/diseasename_Alzheimer.csv'下面表示Alzheimer以及对应的靶点<br>

##2 参数选择##<br>
程序中主要涉及的参数如下所示<br>
low_formula_num = 1 方剂最小中药数目<br>
up_formula_num = 20 方剂最大中药数目<br>
jun_w = 2 君药的权重，下同<br>
chen_w = 1<br>
zuo_w = 1<br>
shi_w = 1<br>

以及"君臣佐使"比例不同：<br>
jun_p = 0.1 君药比例,下同<br>
chen_p =  0.2<br>
zuo_p = 0.3<br>
shi_p =0.4<br>

设定DL和ob阈值<br>
herb_molecules函数可以根据需求设置<br>
herb_mol = herb_mol[(herb_mol['ob']>30) & (herb_mol['drug-likeness']>0.18)]<br>

direct_flag = 0   1表示直接使用靶点,0表示入口为疾病名称<br>

##3 数据经验抽取##<br>
使用者可以根据名医经验、数据库、古代经典名方中抽取对应的方剂集合进行学习，如data文件夹下面的伤寒金匮.csv表示伤寒金匮的方剂。<br>
您可以根据实际情况替换。<br>

##4 data文件夹下的部分数据来自于论文、TCMSP等##<br>
相关引用已经在论文里说明，使用者可以根据需要按照相同的格式增加OMIM，BindingDB等数据库数据。<br>
仅供科学研究使用，请勿用作商业用途。<br>

联系方式：niuqikai@qq.com<br>
