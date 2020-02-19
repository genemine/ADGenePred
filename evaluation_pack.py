#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 11:37:12 2020

@author: lincx
"""

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

times=10000
def p_value(random_list,observed):
    count=0
    for item in random_list:
        if item > observed:
            count = count+1
    if count == 0:
        p_val = 0.0001
    else:
        p_val = count*1.0/times
    return p_val

def write_file(p_val,observed,randomx,f):
    line = 'p-value:'+str(p_val)+'\n'
    print(line)
    f.write(line)
    line = 'observed:'+str(observed)+'\n'
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
    observed = []
    for item in genelist[0:num]:
        if item in score:
            observed.append(score[item])
    random_value = [0]*times
    for i in range(times):
        random_list = genelist[random.sample(range(len(genelist)),num)]
        for item in random_list:
            if item in score:
                random_value[i] = random_value[i]+score[item]
        random_value[i] = random_value[i]/num 
    Score=list(score.values())
    seq_random = std(random_value,Score)
    observed=std(observed,Score)
    observed = np.round(np.mean(observed),4)
    randomx= np.round(np.mean(seq_random),4)
    p_val=p_value(seq_random,observed)
    f=open('../validation_result/seq_evaluation.txt','w')
    write_file(p_val,observed,randomx,f)
    result={}
    result['pvalue']=p_val
    result['observed']=observed
    result['random']=list(seq_random)
    random_score=list(seq_random)
    return p_val,observed,random_score


#2.evaluation of co-expression relationship with known disease genes
def coexpr_validation(param,num,genelist,times):
    print('Step2: evaluation of co-expression relationship with known disease genes')
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
    observed=data[gene].mean()
    random_score=[]
    for i in range(times):
        gene=genelist[random.sample(range(len(genelist)),num)]
        random_score.append(data[gene].mean())
    randomx=np.mean(random_score)
    pval=p_value(random_score,observed)
    print('score of coexpression with AD-associated genes of '+param+' samples')
    f=open('../validation_result/coexpr_score_evaluation_'+param+'.txt','w')
    write_file(pval,observed,randomx,f)
    result={}
    result['pvalue']=pval
    result['observed']=observed
    result['random']=random_score
    return pval,observed,random_score


#observed,random_score=coexpr_validation('AD',num,genelist,times)
#observed,random_score=coexpr_validation('GSE84422',num,genelist,times)
#observed,random_score=coexpr_validation('Control',num,genelist,times)



#Step3: evaluation of ppi interaction with known genes
def ppi_validation(ad_gene,num,genelist,times):
    print('Step3: evaluation of ppi interaction with known genes')
    def stats_ppi(gene):
        stats_ppi1=pd.read_csv('stats_ppi.txt',header=0,index_col=0,sep='\t')
        gene=set(gene).intersection(stats_ppi1.index)
        stats_ppi_gene=stats_ppi1.loc[gene,:]
        return stats_ppi_gene.mean()['score'],stats_ppi_gene.mean()['num']
    gene=genelist[:num]
    observed_score,observed_num=stats_ppi(gene)
    random_score=[]
    random_num=[]
    for i in range(times):
        gene=genelist[random.sample(range(len(genelist)),num)]
        x1,x2=stats_ppi(gene)
        random_score.append(x1)
        random_num.append(x2)
    pval_score=p_value(random_score,observed_score)
    pval_num=p_value(random_num,observed_num)
    randomx=np.mean(random_score)
    f=open('../validation_result/ppi_score_evaluation.txt','w')
    print('score of ppi interaction with AD-associated genes')
    write_file(pval_score,observed_score,randomx,f)
    f=open('../validation_result/ppi_num_evaluation.txt','w')
    print('number of ppi interaction with AD-associated genes')
    randomx=np.mean(random_num)
    write_file(pval_num,observed_num,randomx,f)
    result={}
    result['score']={}
    result['num']={}
    for item in result.keys():
        result[item]['pvalue']=eval('pval_'+item)
        result[item]['observed']=eval('observed_'+item)
        result[item]['random']=eval('random_'+item)
    return pval_score,observed_score,random_score,pval_num,observed_num,random_num


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
    observed={}
    p_val={}
    random_num['all'] = []
    observed['all'] = []
    for item in top3mirna:
        random_num[item] = []
        observed[item] = []
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
        observed[item] = countDict_top[item]
        p_val[item]=p_value(random_num[item],observed[item])
    observed['all']=len(mirna_top.intersection(mirna_ad))
    p_val['all']=p_value(random_num['all'],observed['all'])
    observed=pd.Series(observed)
    random_num=pd.DataFrame(random_num)
    randomx=random_num.mean()
    f=open('../validation_result/mirna_evaluation_'+str(num)+'_genes.txt','w')
    write_file(p_val,observed,randomx,f)
    random_num.to_csv('../validation_result/mirna_interaction_number_of_random_'+str(num)+'_genes.txt',sep='\t')
    return p_val,observed,list(random_num['all'])


def pack():
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
    score=eval(open('sequence_similarity.txt').read())
    gene_mirna=eval(open('gene_mirna.txt').read())
    real_v={}
    ran_rcd={}
    mirna={}
    for num in [100,200,500]:
        print(num)
        real_v[num]={}
        ran_rcd[num]={}
        p_val,real_v[num][0],ran_rcd[num][0]=seq_validation(genelist,num,score)
        p_val,real_v[num][1],ran_rcd[num][1]=coexpr_validation('GSE84422',num,genelist,times)
        p_val,real_v[num][2],ran_rcd[num][2],p_val,real_v[num][3],ran_rcd[num][3],=ppi_validation(ad_gene,num,genelist,times)
        p_val,real_v[num][4],ran_rcd[num][4]= mirna_validation(ad_gene,num,genelist,gene_mirna,times)
    fnew=open('observed.txt','w')
    fnew.write(str(real_v))
    fnew.close()
    fnew=open('random.txt','w')
    fnew.write(str(ran_rcd))
    fnew.close()
pack()

