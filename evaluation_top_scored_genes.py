#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 10:33:18 2020

@author: lincx
"""


import pandas as pd
import numpy as np
import random


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
    obeserved=np.round(obeserved,4)
    line = 'obeserved:'+str(obeserved)+'\n'
    print(line)
    f.write(line)
    randomx=np.round(randomx,4)
    line = 'random:'+str(randomx)+'\n'
    print(line+'\n\n')
    f.write(line)
    f.close()


#1.evaluation of sequence similarity with known genes

def seq_validation(genelist,num):
    print('Step1: evaluation of sequence similarity with known genes')
    score=eval(open('data/sequence_similarity.txt').read())
    obeserved = []
    for item in genelist[0:num]:
        if item in score:
            obeserved.append(score[item])
    random_score= [0]*times
    for i in range(times):
        random_list = genelist[random.sample(range(len(genelist)),num)]
        for item in random_list:
            if item in score:
                random_score[i] = random_score[i]+score[item]
        random_score[i] = random_score[i]/num 
    obeserved = np.round(np.mean(obeserved),4)
    randomx= np.round(np.mean(random_score),4)
    p_val=p_value(random_score,obeserved)
    f=open('validation_result/seq_evaluation.txt','w')
    write_file(p_val,obeserved,randomx,f)
    return p_val,obeserved,random_score



#2.evaluation of co-expression relationship with known disease genes
def coexpr_validation(num,genelist,times):
    print('Step2: evaluation of co-expression relationship with known disease genes')
    data=pd.read_csv('data/coexpression.txt',header=0,index_col=0,sep='\t')
    newlist=[]
    for item in genelist:
        if item in data.index:
            newlist.append(item)
    genelist=np.array(newlist)
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
    print('score of coexpression with AD-associated genes')
    f=open('validation_result/coexpr_score_evaluation.txt','w')
    write_file(pval,obeserved,randomx,f)
    return pval,obeserved,random_score




#Step3: evaluation of ppi interaction with known genes
def ppi_validation(ad_gene,num,genelist,times):
    print('Step3: evaluation of ppi interaction with known genes')
    stats=pd.read_csv('data/stats_ppi.txt',header=0,index_col=0,sep='\t')
    def stats_ppi(gene,stats):
        gene=set(gene).intersection(stats.index)
        stats_ppi_gene=stats.loc[gene,:]
        return stats_ppi_gene.mean()['score'],stats_ppi_gene.mean()['num']
    gene=genelist[:num]
    obeserved_score,obeserved_num=stats_ppi(gene,stats)
    random_score=[]
    random_num=[]
    for i in range(times):
        gene=genelist[random.sample(range(len(genelist)),num)]
        x1,x2=stats_ppi(gene,stats)
        random_score.append(x1)
        random_num.append(x2)
    pval_score=p_value(random_score,obeserved_score)
    pval_num=p_value(random_num,obeserved_num)
    randomx=np.median(random_score)
    f=open('validation_result/ppi_score_evaluation.txt','w')
    print('score of ppi interaction with AD-associated genes')
    write_file(pval_score,obeserved_score,randomx,f)
    f=open('validation_result/ppi_num_evaluation.txt','w')
    print('number of ppi interaction with AD-associated genes')
    randomx=np.median(random_num)
    write_file(pval_num,obeserved_num,randomx,f)




#4.evaluation of AD-associated miRNAs

def mirna_validation(ad_gene,num,genelist,times):
    print('Step4: evaluation of interaction with AD-associated miRNAs')
    gene_mirna=eval(open('data/gene_mirna.txt').read())
    def mirnaOfgenes(genes,gene_mirna):
        mirna_gene = []
        for item in genes:
            if item in gene_mirna:
                mirna_gene = mirna_gene+gene_mirna[item]
        mirna_gene=set(mirna_gene)
        return mirna_gene
    mirna_ad= mirnaOfgenes(ad_gene,gene_mirna)
    random_num = []
    for i in range(times):
        random_list = genelist[random.sample(range(len(genelist)),num)]
        mirna_random= mirnaOfgenes(random_list,gene_mirna)
        random_numi=len(mirna_random.intersection(mirna_ad))
        random_num.append(random_numi)
    mirna_top=mirnaOfgenes(genelist[0:num],gene_mirna)
    obeserved=len(mirna_top.intersection(mirna_ad))
    p_val=p_value(random_num,obeserved)
    randomx=np.mean(random_num)
    f=open('validation_result/mirna_evaluation.txt','w')
    write_file(p_val,obeserved,randomx,f)
    return p_val,obeserved,random_num

### main ###
    
num = 200
print(num)
times=10000
genes_known = pd.read_csv('data/mat_training.txt',sep='\t',index_col=0,header=0).index
genes_all = pd.read_csv('data/prediction.txt',sep='\t',index_col=0,header=0).index
ad_gene=eval(open('data/ad_gene.txt').read())
genelist=[]

#remove known genes to evaluate
for item in genes_all:
    if item not in genes_known:
        genelist.append(item)
genelist = np.array(genelist)

seq_pval,seq_obeserved,seq_random_score=seq_validation(genelist,num)

coexp_pval,coexp_obeserved,coexp_random_score=coexpr_validation(num,genelist,times)

ppi_validation(ad_gene,num,genelist,times)

miRNA_pval,miRNA_obeserved,miRNA_random_num=mirna_validation(ad_gene,num,genelist,times)      


