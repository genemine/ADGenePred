#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 23:09:39 2020

@author: lincx
"""


import pandas as pd
import numpy as np
import random

def to_list(inp):
    f=open(inp)
    result=[]
    for line in f:
        t=line.strip()
        result.append(t)
    return result

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
    f.write(line)
    obeserved=np.round(obeserved,4)
    line = 'obeserved:'+str(obeserved)+'\n'
    f.write(line)
    randomx=np.round(randomx,4)
    line = 'random:'+str(randomx)+'\n'
    print('p-value:'+str(p_val)+',obeserved:'+str(obeserved)+','+'random:'+str(randomx))
    f.write(line)
    f.close()


#1.evaluation of sequence similarity with known genes

def seq_validation(genelist1,genelist2,num,times):
    print('Step1: evaluation of sequence similarity with known genes')
    score=eval(open('data/sequence_similarity.txt').read())
    obeserved = []
    for item in genelist1[0:num]:
        if item in score:
            obeserved.append(score[item])
    random_score= [0]*times
    for i in range(times):
        random_list = genelist2[random.sample(range(len(genelist2)),num)]
        n=0
        for item in random_list:
            if item in score:
                random_score[i] = random_score[i]+score[item]
                n=n+1
        random_score[i] =random_score[i]/n
    #obeserved = np.round(np.mean(obeserved),4)
    randomx= np.round(np.mean(random_score),4)
    p_val=p_value(random_score,np.round(np.mean(obeserved),4))
    f=open('validation_result/seq_evaluation.txt','w')
    write_file(p_val,obeserved,randomx,f)
    return p_val,obeserved,random_score



#2.evaluation of co-expression relationship with known disease genes
def coexpr_validation(data,num,genelist1,genelist2,times):
    print('Step2: evaluation of co-expression relationship with known disease genes')
    data=data.abs().mean(axis=1)
    gene1=genelist1[0:num]
    gene=list(set(gene1).intersection(data.index))
    obeserved=data[gene].mean()
    random_score=[]
    for i in range(times):
        gene=genelist2[random.sample(range(len(genelist2)),num)]
        gene=list(set(gene).intersection(data.index))
        random_score.append(data[gene].mean())
    randomx=np.median(random_score)
    pval=p_value(random_score,obeserved)
    print('score of coexpression with AD-associated genes')
    f=open('validation_result/coexpr_score_evaluation.txt','w')
    write_file(pval,obeserved,randomx,f)
    return pval,obeserved,random_score





#Step3: evaluation of ppi interaction with known genes
def ppi_validation1(ad_gene,num,genelist1,genelist2,times):
    print('Step3: evaluation of ppi interaction with known genes')
    stats=eval(open('data/stats_ppi_Bioplex.txt').read())
    def stats_ppi(gene,stats):
        gene=set(gene).intersection(stats.keys())
        interacts=[]
        for genei in gene:
            interacts=interacts+stats[genei]
        interact_num=len(interacts)
        return interact_num
    gene=genelist1[:num]
    obeserved_num=stats_ppi(gene,stats)
    random_num=[]
    for i in range(times):
        gene=genelist2[random.sample(range(len(genelist2)),num)]
        x2=stats_ppi(gene,stats)
        random_num.append(x2)
    pval_num=p_value(random_num,obeserved_num)
    f=open('validation_result/ppi_num_evaluation.txt','w')
    print('number of ppi interaction with AD-associated genes')
    randomx=np.median(random_num)
    write_file(pval_num,obeserved_num,randomx,f)
    return pval_num,obeserved_num,random_num

def ppi_validation2(ad_gene,num,genelist1,genelist2,times):
    print('Step3: evaluation of ppi interaction with known genes')
    stats=eval(open('data/stats_ppi_human_interactome.txt').read())
    def stats_ppi(gene,stats):
        gene=set(gene).intersection(stats.keys())
        interacts=[]
        for genei in gene:
            interacts=interacts+stats[genei]
        interact_num=len(interacts)
        return interact_num
    gene=genelist1[:num]
    obeserved_num=stats_ppi(gene,stats)
    random_num=[]
    for i in range(times):
        gene=genelist2[random.sample(range(len(genelist2)),num)]
        x2=stats_ppi(gene,stats)
        random_num.append(x2)
    pval_num=p_value(random_num,obeserved_num)
    f=open('validation_result/ppi_num_evaluation.txt','w')
    print('number of ppi interaction with AD-associated genes')
    randomx=np.median(random_num)
    write_file(pval_num,obeserved_num,randomx,f)
    return pval_num,obeserved_num,random_num
def ppi_validation3(ad_gene,num,genelist1,genelist2,times):
    print('Step3: evaluation of ppi interaction with known genes')
    stats=eval(open('data/stats_ppi_string.txt').read())
    def stats_ppi(gene,stats):
        gene=set(gene).intersection(stats.keys())
        interacts=[]
        for genei in gene:
            interacts=interacts+stats[genei]
        interact_num=len(set(interacts))
        return interact_num
    gene=genelist1[:num]
    obeserved_num=stats_ppi(gene,stats)
    random_num=[]
    for i in range(times):
        gene=genelist2[random.sample(range(len(genelist2)),num)]
        x2=stats_ppi(gene,stats)
        random_num.append(x2)
    pval_num=p_value(random_num,obeserved_num)
    f=open('validation_result/ppi_num_evaluation.txt','w')
    print('number of ppi interaction with AD-associated genes')
    randomx=np.median(random_num)
    write_file(pval_num,obeserved_num,randomx,f)
    return pval_num,obeserved_num,random_num

#4.evaluation of AD-associated miRNAs

def mirna_validation(ad_gene,num,genelist1,genelist2,times):
    print('Step4: evaluation of interaction with AD-associated miRNAs')
    gene_mirna=eval(open('data/gene2miRNAs.txt').read())
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
        random_list = genelist2[random.sample(range(len(genelist2)),num)]
        mirna_random= mirnaOfgenes(random_list,gene_mirna)
        random_numi=len(mirna_random.intersection(mirna_ad))
        random_num.append(random_numi)
    mirna_top=mirnaOfgenes(genelist1[0:num],gene_mirna)
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
ad_gene=eval(open('data/ad_gene.txt').read())
genelist1=pd.read_csv('data/prediction.txt',sep='\t',index_col=0,header=None).index
inp='data/genes.txt'
genelist2=to_list(inp)


genelist1 = np.array(genelist1)
genelist2 = np.array(genelist2)

def step1():
    seq_pval,seq_obeserved,seq_random_score=seq_validation(genelist1,genelist2,num,times)
    result=[seq_pval,seq_obeserved,seq_random_score]
    fnew=open('validation_result/seq_result.txt','w')
    fnew.write(str(result))
    fnew.close()

def step2():
    data=pd.read_csv('data/coexpression.txt',header=0,index_col=0,sep='\t')
    coexp_pval,coexp_obeserved,coexp_random_score=coexpr_validation(data,num,genelist1,genelist2,times)
    result=[coexp_pval,coexp_obeserved,coexp_random_score]
    fnew=open('validation_result/coexp_result.txt','w')
    fnew.write(str(result))
    fnew.close()


def step3():
    pval_num,obeserved_num,random_num=ppi_validation1(ad_gene,num,genelist1,genelist2,times)
    result=[pval_num,obeserved_num,random_num]
    fnew=open('validation_result/ppi_result_Bioplex.txt','w')
    fnew.write(str(result))
    fnew.close()
    pval_num,obeserved_num,random_num=ppi_validation2(ad_gene,num,genelist1,genelist2,times)
    result=[pval_num,obeserved_num,random_num]
    fnew=open('validation_result/ppi_result_human_interactome.txt','w')
    fnew.write(str(result))
    fnew.close()
    pval_num,obeserved_num,random_num=ppi_validation3(ad_gene,num,genelist1,genelist2,times)
    result=[pval_num,obeserved_num,random_num]
    fnew=open('validation_result/ppi_result_string.txt','w')
    fnew.write(str(result))
    fnew.close()

def step4():
    miRNA_pval,miRNA_obeserved,miRNA_random_num=mirna_validation(ad_gene,num,genelist1,genelist2,times)     
    result=[miRNA_pval,miRNA_obeserved,miRNA_random_num]
    fnew=open('validation_result/mirna_result.txt','w')
    fnew.write(str(result))
    fnew.close()

step1()
step2()
step3()
step4()



