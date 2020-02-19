#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 10:33:18 2020

@author: lincx
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 20:56:39 2020

@author: lincx
"""

import pandas as pd
import numpy as np
import random
from collections import Counter

num = 200
print(num)
times=10000
genes_known = pd.read_csv('mat_training.txt',sep='\t',index_col=0,header=0).index
#genes_all = pd.read_csv('prediction.txt',sep=',',index_col=0,header=0).index
genes_all = pd.read_csv('rank.csv',sep=',',index_col=0,header=0).index
ad_gene=eval(open('ad_gene.txt').read())
genelist=[]

#remove known genes to evaluate
for item in genes_all:
    if item not in genes_known:
        genelist.append(item)
genelist = np.array(genelist)

def p_value(random_list,obeserved):
    count=0
    for item in random_list:
        if item > obeserved:
            count = count+1
    if count == 0:
        p_val = 0.0001
    else:
        p_val = count*1.0/times
    return p_val

def write_file(p_val,obeserved,randomx,f):
    line = 'p-value:'+str(p_val)+'\n'
    print(line)
    f.write(line)
    line = 'obeserved:'+str(obeserved)+'\n'
    print(line)
    f.write(line)
    line = 'random:'+str(randomx)+'\n'
    print(line+'\n\n')
    f.write(line)
    f.close()


#1.evaluation of sequence similarity with known genes

def seq_validation(genelist,num,score):
    print('Step1: evaluation of sequence similarity with known genes')
    def std(x,y):
        x = (x-np.min(y))/(np.max(y)-np.min(y))
        return x
    obeserved = []
    for item in genelist[0:num]:
        if item in score:
            obeserved.append(score[item])
    random_value = [0]*times
    for i in range(times):
        random_list = genelist[random.sample(range(len(genelist)),num)]
        for item in random_list:
            if item in score:
                random_value[i] = random_value[i]+score[item]
        random_value[i] = random_value[i]/num 
    Score=list(score.values())
    seq_random = std(random_value,Score)
    obeserved=std(obeserved,Score)
    obeserved = np.round(np.mean(obeserved),4)
    randomx= np.round(np.mean(seq_random),4)
    p_val=p_value(seq_random,obeserved)
    f=open('../validation_result/seq_evaluation.txt','w')
    write_file(p_val,obeserved,randomx,f)
    return seq_random,obeserved,p_val


score=eval(open('/Users/lincx/Documents/study/paper/AD_pred/github/data/sequence_similarity.txt').read())
seq_random,obeserved,p_val=seq_validation(genelist,num,score)

#2.evaluation of co-expression relationship with known disease genes
def coexpr_validation(param,num,genelist,times):
    data=pd.read_csv(param+'-coexpression.txt',header=0,index_col=0,sep='\t')
    newlist=[]
    for item in genelist:
        if item in data.index:
            newlist.append(item)
    genelist=np.array(newlist)
    print(len(data.index))
    data=data.loc[list(genelist),:]
    data=data.abs().mean(axis=1)
    gene=genelist[:num]
    obeserved=data[gene].mean()
    random_score=[]
    for i in range(times):
        gene=genelist[random.sample(range(len(genelist)),num)]
        random_score.append(data[gene].mean())
    randomx=np.median(random_score)
    pval=p_value(random_score,obeserved)
    print('score of coexpression with AD-associated genes of '+param+' samples')
    f=open('../validation_result/coexpr_score_evaluation_'+param+'.txt','w')
    write_file(pval,obeserved,randomx,f)
    return obeserved,random_score

print('Step2: evaluation of co-expression relationship with known disease genes')

obeserved,random_score=coexpr_validation('GSE84422',num,genelist,times)
#obeserved,random_score=coexpr_validation('Control',num,genelist,times)



#Step3: evaluation of ppi interaction with known genes
def ppi_validation(ad_gene,num,genelist,times):
    print('Step3: evaluation of ppi interaction with known genes')
    def stats_ppi(gene):
        stats_ppi=pd.read_csv('stats_ppi.txt',header=0,index_col=0,sep='\t')
        gene=set(gene).intersection(stats_ppi.index)
        stats_ppi_gene=stats_ppi.loc[gene,:]
        return stats_ppi_gene.mean()['score'],stats_ppi_gene.mean()['num']
    gene=genelist[:num]
    obeserved_score,obeserved_num=stats_ppi(gene)
    random_score=[]
    random_num=[]
    for i in range(times):
        gene=genelist[random.sample(range(len(genelist)),num)]
        x1,x2=stats_ppi(gene)
        random_score.append(x1)
        random_num.append(x2)
    pval_score=p_value(random_score,obeserved_score)
    pval_num=p_value(random_num,obeserved_num)
    randomx=np.median(random_score)
    f=open('../validation_result/ppi_score_evaluation.txt','w')
    print('score of ppi interaction with AD-associated genes')
    write_file(pval_score,obeserved_score,randomx,f)
    f=open('../validation_result/ppi_num_evaluation.txt','w')
    print('number of ppi interaction with AD-associated genes')
    randomx=np.median(random_num)
    write_file(pval_num,obeserved_num,randomx,f)


#ppi_validation(ad_gene,num,genelist,times)



#4.evaluation of AD-associated miRNAs

def mirna_validation(ad_gene,num,genelist,gene_mirna,times):
    print('Step4: evaluation of interaction with AD-associated miRNAs')
    def mirnaOfgenes(genes,gene_mirna):
        mirna_gene = []
        for item in genes:
            if item in gene_mirna:
                mirna_gene = mirna_gene+gene_mirna[item]
        countDict = Counter(mirna_gene)
        mirna_gene=set(mirna_gene)
        return mirna_gene,countDict
    mirna_ad,countDict_ad = mirnaOfgenes(ad_gene,gene_mirna)
    x=sorted(countDict_ad.items(), key=lambda item:item[1], reverse=True)
    top3mirna=[x[0][0],x[1][0],x[2][0]]
    #for item in countDict_ad:
    #print(mirna_ad)
    random_num={}
    obeserved={}
    p_val={}
    random_num['all'] = []
    obeserved['all'] = []
    for item in top3mirna:
        random_num[item] = []
        obeserved[item] = []
    for i in range(times):
        random_list = genelist[random.sample(range(len(genelist)),num)]
        mirna_random,countDict_random = mirnaOfgenes(random_list,gene_mirna)
        for item in top3mirna:
            if item in countDict_random:
                random_num[item].append(countDict_random[item])
            else:
                random_num[item].append(0)
        random_numi=len(mirna_random.intersection(mirna_ad))
        random_num['all'].append(random_numi)
    mirna_top,countDict_top=mirnaOfgenes(genelist[0:num],gene_mirna)
    for item in top3mirna:
        obeserved[item] = countDict_top[item]
        p_val[item]=p_value(random_num[item],obeserved[item])
    obeserved['all']=len(mirna_top.intersection(mirna_ad))
    p_val['all']=p_value(random_num['all'],obeserved['all'])
    obeserved=pd.Series(obeserved)
    random_num=pd.DataFrame(random_num)
    randomx=random_num.median()
    f=open('../validation_result/mirna_evaluation.txt','w')
    write_file(p_val,obeserved,randomx,f)
    return p_val,obeserved,random_num

#gene_mirna=eval(open('gene_mirna.txt').read())
#p_val,obeserved,random_num = mirna_validation(ad_gene,num,genelist,gene_mirna,times)      


