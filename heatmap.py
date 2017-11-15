#!/bin/python
##################################
import pandas as pd
import sys
from sys import*
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np
####################################
## Make your heatmap look nice####
sns.set_context('poster')
sns.set_style('white')
cmap = sns.diverging_palette(240,10,s=90,as_cmap=True,center='light')
###########################################################
########### Load the file ################################
Tissue = pd.read_csv(argv[1], sep = '\t', header=0)
Tissue=Tissue.drop_duplicates()
Tissue.index= Tissue[Tissue.columns[0]]
###########################################################
## You will need to reoder the rows in order to make your final heatmap look nice####
Tissue = Tissue[['c2','c3','c4','c5','c6','c7','c9','c10','c11','c12','c13','c14','c16']]
original_cols = Tissue.columns 
Tissue = Tissue.dropna()

## Another way of standardizing the rows Log Expression - Log (Mean Expression) ##
#Tissue['avg'] = Tissue.mean(axis=1)
#for item in original_cols:
#    Tissue[item] = Tissue[item].apply(lambda x: np.log10(x+1)) - Tissue['avg'].apply(lambda x: np.log(x))

## First we plot the zscore standardized figure and cluster the rows
fig = plt.figure()
g=sns.clustermap(Tissue,row_cluster= True, col_cluster=False,linewidth=0, z_score=0, cmap=cmap)
pd.DataFrame(g.data2d).to_csv(argv[2]+'.txt', sep = '\t')
plt.suptitle(argv[3])
plt.savefig(argv[2]+'.png')

## We take the matrix from above and plot the corresponding raw umi counts 
fig = plt.figure() 
r=sns.clustermap(Tissue.ix[g.data2d.index], row_cluster=False, col_cluster=False, linewidth=0, cmap = cmap)
pd.DataFrame(r.data2d).to_csv(argv[2]+'raw.txt', sep = '\t') 
plt.suptitle(argv[3]+' Raw')
plt.savefig(argv[2]+'raw.png')
