#!/bin/python
import pandas as pd
import sys
from sys import*
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np


# Set your context to poster to make your heatmap look nice  
sns.set_context('poster')
sns.set_style('white')

# Red, blue, white color palette 
cmap = sns.diverging_palette(240,10,s=90,as_cmap=True,center='light')

#' read_tissue takes in a txt file and returns a data frame with index of genes, columns of clusters 
#' and values of average UMI in each cluster 
#' if log_norm = True, then the data is normalized for each gene by taking the log of the average UMI in a cluster, 
#' and subtracting the log of the average UMI for that gene across all clusters 

def read_tissue(tissue_df, log_norm=False):
  tissue_df = pd.read_csv(tissue_df,sep = '\t', header=0)
  tissue_df.drop_duplicates(inplace=True)
  tissue_df.dropna(inplace=True)
  tissue_df.index= tissue_df[tissue_df.columns[0]]
  if log_norm = True:
    tissue_df['avg'] = tissue_df.mean(axis=1)
    for item in tissue_df.columns:
      tissue_df[item] = tissue_df[item].apply(lambda x: np.log10(x+1)) - tissue_df['avg'].apply(lambda x: np.log10(x))
  reutrn (tissue_df)

#' plot_save_figs takes in a clustermap 'hmap' ; a title for the heatmap 'title'; and the name to save the png heatmap 'pic_name'
#' plot_save figs returns a png file
  
def plot_save_figs(title, pic_name, hmap):
  fig = plt.figure()
  hmap
  plt.suptitle(title)
  plt.savefig(pic_name+'.png')

#' plot_heatmaps takes in a dataframe 'processed_tissue'; a name for saving your figures 'output_name', a title fo your plot 'plt_title'; 
#' a list of column names 'col_order' (the default col_order is the columns of the data frame)
#' the col_order needs to be changes in order to enhance the aesthetics of the heatmap 
#' plot_heatmaps outputs two txt files (the clustered dataframe with raw umi values; the clustered dataframe wtih row z-scored values)
#' plot_heatmaps outputs two png files (a clustermap of the data frame using raw values; a clustermap of the data frame using 
#' row z scored values 

def plot_heatmaps (processed_tissue, output_name, plt_title, col_order=processed_tissue.columns[1:]):
  col_order= processed_tissue[col_order] 
  scaled=sns.clustermap(processed_tissue,row_cluster= True, col_cluster=False,linewidth=0, z_score=0, cmap=cmap)
  raw=sns.clustermap(processed_tissue.ix[scaled.data2d.index], row_cluster=False, col_cluster=False, linewidth=0, cmap = cmap)
  pd.DataFrame(scaled.data2d).to_csv(output_name+'.txt', sep = '\t')
  pd.DataFrame(raw.data2d).to_csv(output_name+'raw.txt', sep = '\t') 
  plot_save_figs(plt_title, output_name, scaled)
  plot_save_figs(plt_title+'raw', output_name+'raw', raw)
  

#' Usage 
#' argv[1] is the txt file with rows as genes, columns as clusters, and values as average UMI for each gene in a given cluster
#' argv[2] is a name for what you want to save your heatmaps as 
#' argv[3] is a name for what you want to title your heatmap 

plot_heatmaps(read_tissue(argv[1]) , argv[2], argv[3])

