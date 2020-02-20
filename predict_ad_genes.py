#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  15 09:08:55 2020

@author: lincx
"""

import numpy as np
import pandas as pd
import random
from sklearn import svm
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn.model_selection import StratifiedKFold,GridSearchCV
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve,aucs
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve


data=pd.read_csv('mat_training.txt',index_col=0,header=0,sep='\t')
data_pred=pd.read_csv('mat_all_genes.txt',index_col=0,header=0,sep='\t')

#after parameter search: GridSearchCV
models={}
models['ExtraTrees']=ExtraTreesClassifier(n_estimators = 500,class_weight = "balanced_subsample", n_jobs = -1, bootstrap = True)
models['RandomForest']=RandomForestClassifier(n_estimators = 500,class_weight = "balanced_subsample", n_jobs = -1, bootstrap = True)
models['SVM']=svm.SVC(kernel = "linear",probability = True,class_weight = "balanced",C=0.1,gamma=10)
models['LogisticRegression']=LogisticRegression(solver='lbfgs',multi_class='multinomial',class_weight={0:0.5, 1:0.5})
   
#cross validation score
AUROC={}
AUPRC={}
for modeli in models:
    AUROC[modeli]=[]
    AUPRC[modeli]=[]
#score: to store prediction of four models
score=pd.DataFrame(np.zeros((len(data_pred),4)))
score.index=data_pred.index
score.columns=models.keys()


times=100
for i in range(times):
    print('step:'+str(i))
    index=random.sample(list(data.index),len(data.index))
    data_new=data.loc[index,:]
    X=data_new.iloc[:,0:-1]
    y=data_new.iloc[:,-1]
    skf = StratifiedKFold(n_splits=5)
    for modeli in models:
        y_pred=pd.Series([])
        y_new=pd.Series([])
        for train,test in skf.split(X,y):
            X_train, X_test = X.iloc[train,:], X.iloc[test,:]
            y_train, y_test = y.iloc[train], y.iloc[test]
            model=models[modeli]
            model.fit(X_train,y_train)
            pred=model.predict_proba(X_test)[:,1]
            Pred=model.predict_proba(data_pred)[:,1]
            y_predi=pd.Series(pred)
            y_predi.index=y_test.index
            y_new=y_new.append(y_test)
            y_pred=y_pred.append(y_predi)
            score[modeli]=score[modeli]+Pred
        AUROC[modeli].append(round(roc_auc_score(y_new, y_pred),4))
        AUPRC[modeli].append(round(average_precision_score(y_new, y_pred),4))
        print(modeli+':'+str(AUROC[modeli]))
        print(modeli+':'+str(AUPRC[modeli]))

#model performance
print(AUROC)
print(AUPRC)    
result={}
result['AUROC']=AUROC
result['AUPRC']=AUPRC
fnew=open('model_performace.txt','w')
fnew.write(str(result))
fnew.close()

#calculate the average of 100 times prediction
pred=score/(times*5)
pred=pred.sort_values(by='ExtraTrees',ascending=False)
pred.to_csv('prediction.txt',sep='\t')

