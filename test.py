import pandas as pd
import os
filepath = 'D:/刘思鸿/1995卫生部标准/'

file_rs = 'cintmed_rs.txt'

#namelist = ['标准号','药名','藏方药名','汉语拼音','英文名','藏方药分类','药材基原','性状','炮制','性味','功效','主治','用法用量','贮藏','鉴别','检查']
namelist = []

def getnamelist(filename):
    with open(filename ,encoding='utf-8') as fl:
        for line in fl:
            if '】:' in line:
                lines = line.split('】:')
                #keybefore = key
                key = lines[0].replace('【','')
                if key not in namelist:
                    namelist.append(key)


def writefiletodata(filename,i):
    with open(filename ,encoding='utf-8') as fl:
        for line in fl:
            if '----------' in line :
                if i > 0:
                #写入上一条数据
                    print(dictvalue)
                    with open(file_rs, 'a') as fl:
                        for name in namelist:
                            if name in dictvalue:
                                fl.write(str(dictvalue[name]))
                                fl.write('*')
                            else:
                                fl.write(' ')
                                fl.write('*')
                        fl.write('\n')

                #for namekey in namelist:

            #print(line.replace('----------',''))
                dictvalue = {}

            else:
                i = i +1
                if '】:' in line:
                    lines = line.split('】:')
                    #keybefore = key
                    key = lines[0].replace('【','')
                    value = lines[1].strip()
                    dictvalue[key] = value
                else :
                    value = value + line.strip()
                    dictvalue[key] = value
        with open(file_rs, 'a') as fl:
            for name in namelist:
                if name in dictvalue:
                    fl.write(str(dictvalue[name]))
                    fl.write('*')
                else:
                    fl.write('  ')
                    fl.write('*')
            fl.write('\n')


for root, dirs, files in os.walk(filepath):

    # root 表示当前正在访问的文件夹路径
    # dirs 表示该文件夹下的子目录名list
    # files 表示该文件夹下的文件list

    # 遍历文件
    for f in files:
        filename = os.path.join(root, f)
        getnamelist(filename)
#print(namelist)


for root, dirs, files in os.walk(filepath):

    # root 表示当前正在访问的文件夹路径
    # dirs 表示该文件夹下的子目录名list
    # files 表示该文件夹下的文件list

    # 遍历文件
    for f in files:
        filename = os.path.join(root, f)
        writefiletodata(filename,0)


