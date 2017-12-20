import sys
from sys import*
import pandas as pd
import numpy as np 
import seaborn as sns
import matplotlib as mpl
from matplotlib import pyplot

sns.set_style('white')
sns.set_context('poster')

def merged_data(expression, metadata):
    umi = pd.read_csv(expression, sep = '\t', header=0)
    meta = pd.read_csv(metadata, sep = '\t', header=0)
    merge= pd.concat([umi, meta.T], axis=0)
    merge= merge.T
    merge[merge.columns[-1]]= merge[merge.columns[-1]].apply(lambda x: float(x))
    return merge 

merged= merged_data(argv[1], argv[2])

fig = plt.figure()
ax1 = plt.subplot(241)
ax2 = plt.subplot(242)
ax3 = plt.subplot(243)
ax4 = plt.subplot(244)
ax5 = plt.subplot(245)
ax6 = plt.subplot(246)
ax7 = plt.subplot(247)
ax8 = plt.subplot(248)
ax= [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
cmap= sns.color_palette('husl',9)

for i in range(len(merged.columns)-1):
    sns.violinplot(data=[merged[merged.V1==0][merged.columns[i]],merged[merged.V1==1][merged.columns[i]],\
                         merged[merged.V1==2][merged.columns[i]],merged[merged.V1==3][merged.columns[i]], \
                         merged[merged.V1==4][merged.columns[i]], merged[merged.V1==5][merged.columns[i]],\
                         merged[merged.V1==6][merged.columns[i]], merged[merged.V1==7][merged.columns[i]],\
                         merged[merged.V1==8][merged.columns[i]], merged[merged.V1==9][merged.columns[i]]],\
                   orient='h', cut=0,ax=ax[i], bw=2, scale='area', palette=cmap)
    ax[i].set_xlim(0,10)
    ax[i].set_title(merged.columns[i])
    for i in [1,2,3,5,6,7]:
        plt.setp(ax[i].get_yticklabels(), visible=False)

fig.savefig(argv[3]+'.pdf')

sns.despine()
