pair_score = {'12':0,'13':4,'14':3,'15':3,'16':2,'23':1,'24':0,'25':1,'26':1,'34':3,'35':4,'36':3,'45':2,'46':2,'56':2}
hscore = [7/8,5/8,15/8,0,1,3/8]

f1 = [1,2,3,4,6]
fmapscore = 0
for i in range(len(f1)-1):
    for j in range(i+1,len(f1)):
        if i!=j:
            print(i,j,hscore[i-1] * hscore[j-1] * pair_score[str(i+1)+str(j+1)])
            if str(f1[i])+str(f1[j]) in pair_score:
                fmapscore = fmapscore + hscore[f1[i]-1] * hscore[f1[j]-1] * pair_score[str(f1[i])+str(f1[j])]
print(fmapscore/len(f1))


