import sys
from sys import*
import pandas as pd
import numpy as np 
import seaborn as sns
import matplotlib as mpl
from matplotlib import pyplot

sns.set_style('white')
sns.set_context('poster')

#' merged_data takes in a data frame with rows of genes of your choice (marker genes), columns of cells, and values umi; 
#' another data frame of 'meta data' has rows of cells, and one column of cluster ids 
#' It transposes the expression data frame to (cells x genes), and adds 'cluster id' as an additional column

def merged_data(expression, metadata):
    umi = pd.read_csv(expression, sep = '\t', header=0)
    meta = pd.read_csv(metadata, sep = '\t', header=0)
    merge= pd.concat([umi, meta.T], axis=0)
    merge= merge.T
    merge[merge.columns[-1]]= merge[merge.columns[-1]].apply(lambda x: float(x))
    return merge 

#' This plots the distribution of each marker gene as a violinplot in a subplot of the main figure for 8 genes and 9 cluster ids
#' merged is the merged data frame where rows are cells, columns are genes, and the last column is the cluster id
#' cluster_ids is a list of the cluster_ids ie: [0,1,2,3,4,5,6,7,8,9]
#' figure_name is the name of the pdf you want to save the violinplots as 
def plot_marker_gene_vplots(merged, cluster_ids, figure_name):
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

    #' V1 is the cluster ID column (in this example there are 10 clusters)
    for i in range(len(merged.columns)-1):
        sns.violinplot(data=[merged[merged.V1==cluster_ids[0]][merged.columns[i]],merged[merged.V1==cluster_ids[1]][merged.columns[i]],\
                          merged[merged.V1==cluster_ids[2]][merged.columns[i]],merged[merged.V1==cluster_ids[3]][merged.columns[i]], \
                          merged[merged.V1==cluster_ids[4]][merged.columns[i]], merged[merged.V1==cluster_ids[5]][merged.columns[i]],\
                          merged[merged.V1==cluster_ids[6]][merged.columns[i]], merged[merged.V1==cluster_ids[7]][merged.columns[i]],\
                          merged[merged.V1==cluster_ids[8]][merged.columns[i]], merged[merged.V1==cluster_ids[9]][merged.columns[i]]],\
                       orient='h', cut=0,ax=ax[i], bw=2, scale='area', palette=cmap)
        ax[i].set_xlim(0,10)
        ax[i].set_title(merged.columns[i])
    #'don't want y labels for all the plots (since y labels are cluster IDs) 
        for i in [1,2,3,5,6,7]:
            plt.setp(ax[i].get_yticklabels(), visible=False)
    sns.despine()
    fig.savefig(figure_name +'.pdf')
  
#' Usage for 8 genes and a dataset with 10 clusters 
plot_marker_gene_vplots(merged_data(argv[1], argv[2]),range(10), argv[3])

