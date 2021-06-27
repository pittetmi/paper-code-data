# Functions 

import os,sys
import scanpy as sc
import pandas as pd
import numpy as np
import scipy as sp
import anndata
import itertools
import functools
from time import time
import glob
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
from collections import Counter
import multiprocessing
from collections import defaultdict 

######## centroids 

# From Rapolas Zilionis: 
def adata_centroids(label,adata,E=None,gene_list=None):
    
    """
    Calculate average gene expression level per cell label (e.g. cluster).
    input:
        - label: name of column that stores the label of interest in adata.obs
        - adata: AnnData object OR a cell x feature pandas dataframe with label as one of the columns
        - E and gene_list: optional and only used when adata is not an AnnData object. In that case
        the cells x genes sparse expression matrix E and the gene_list must be specified
        
    returns:
        pandas dataframe, centroids x genes
        
    """
    
    if isinstance(adata,anndata.AnnData):
        E = adata.X
        gene_list = adata.var_names
        meta = adata.obs
    else:
        meta = adata
        
        
    labs = meta[label].unique()
    centroids = {}
    for lab in labs:
        msk = (meta[label] == lab).values #use "values" to turn pd.Series into row-label-less np.array,
                                             #sometimes row labels mess up the order

        centroids[lab] = np.array(E[msk,:].mean(axis=0))[0]
    centroids=pd.DataFrame(centroids).T
    centroids.columns = gene_list
    return centroids

######## gmask  

def gmask_outlier_genes(adata, min_counts = 5, min_cells = 5):
    
    gmask = np.array((adata.X >= min_counts).sum(axis=0))[0] >= min_cells
    print(str(gmask.sum()), 'genes are expressed in at least', str(min_cells), 'cells with at least more than', str(min_counts), 'counts per cell\n')
  
    return gmask


######## cmask 

def filter_cells_adata(adata, return_filtered = True, **kwargs):
    
    # Currently supports 3 filter columns 
    
    # Define kwargs 
    flist1 = kwargs.get('flist1', list())
    flist2 = kwargs.get('flist2', list())
    flist3 = kwargs.get('flist3', list())
    
    col1 = kwargs.get('col1', None)
    col2 = kwargs.get('col2', None)
    col3 = kwargs.get('col3', None)
    
    isnot1 = kwargs.get('isnot1', False)
    isnot2 = kwargs.get('isnot2', False)
    isnot3 = kwargs.get('isnot3', False)
    
    # Check whether columns exsist 
    if col1:
        if not col1 in adata.obs.columns: 
            raise Exception('{} is not a column in .obs'.format(col1))
    
    if col2:
        if not col2 in adata.obs.columns: 
            raise Exception('{} is not a column in .obs'.format(col2))
        
    if col3:
        if not col3 in adata.obs.columns: 
            raise Exception('{} is not a column in .obs'.format(col3))
        

    # Define the cmask: 
    if flist1 and flist2 and flist3: 
        
        if set(flist1).issubset(set(adata.obs[col1])) == False: 
            raise Exception('flist1 is not a subset of {}'.format(col1))
        
        if set(flist2).issubset(set(adata.obs[col2])) == False: 
            raise Exception('flist2 is not a subset of {}'.format(col2))
        
        if set(flist3).issubset(set(adata.obs[col3])) == False: 
            raise Exception('flist3 is not a subset of {}'.format(col3))
        
        cmasks1 = adata.obs[col1].isin(flist1)
        cmasks2 = adata.obs[col2].isin(flist2)
        cmasks3 = adata.obs[col3].isin(flist3)
        
        if isnot1 and isnot2 and isnot3:
            cmask = ~cmasks1 & ~cmasks2 & ~cmasks3
            print((1-cmasks1).sum(), 'cells will survive {} filtering.'.format(col1))
            print((1-cmasks2).sum(), 'cells will survive {} filtering.'.format(col2))
            print((1-cmasks3).sum(), 'cells will survive {} filtering.'.format(col3))
            print(cmask.sum(), 'cells will survive all filters.\n')
            
        elif isnot1 and isnot2: 
            cmask = ~cmasks1 & ~cmasks2 & cmasks3
            print((1-cmasks1).sum(), 'cells will survive {} filtering.'.format(col1))
            print((1-cmasks2).sum(), 'cells will survive {} filtering.'.format(col2))
            print(cmasks3.sum(), 'cells will survive {} filtering.'.format(col3))
            print(cmask.sum(), 'cells will survive all filters.\n')
        
        elif isnot1 and isnot3: 
            cmask = ~cmasks1 & cmasks2 & ~cmasks3
            print((1-cmasks1).sum(), 'cells will survive {} filtering.'.format(col1))
            print(cmasks2.sum(), 'cells will survive {} filtering.'.format(col2))
            print((1-cmasks3).sum(), 'cells will survive {} filtering.'.format(col3))
            print(cmask.sum(), 'cells will survive all filters.\n')
        
        elif isnot2 and isnot3: 
            cmask = cmasks1 & ~cmasks2 & ~cmasks3
            print(cmasks1.sum(), 'cells will survive {} filtering.'.format(col1))
            print((1-cmasks2).sum(), 'cells will survive {} filtering.'.format(col2))
            print((1-cmasks3).sum(), 'cells will survive {} filtering.'.format(col3))
            print(cmask.sum(), 'cells will survive all filters.\n')
        
        elif isnot1:
            cmask = ~cmasks1 & cmasks2 & cmasks3
            print((1-cmasks1).sum(), 'cells will survive {} filtering.'.format(col1))
            print(cmasks2.sum(), 'cells will survive {} filtering.'.format(col2))
            print(cmasks3.sum(), 'cells will survive {} filtering.'.format(col3))
            print(cmask.sum(), 'cells will survive all filters.\n')
        
        elif isnot2:
            cmask = cmasks1 & ~cmasks2 & cmasks3 
            print(cmasks1.sum(), 'cells will survive {} filtering.'.format(col1))
            print((1-cmasks2).sum(), 'cells will survive {} filtering.'.format(col2))
            print(cmasks3.sum(), 'cells will survive {} filtering.'.format(col3))
            print(cmask.sum(), 'cells will survive all filters.\n')
            
        elif isnot3:
            cmask = cmasks1 & cmasks2 & ~cmasks3
            print(cmasks1.sum(), 'cells will survive {} filtering.'.format(col1))
            print(cmasks2.sum(), 'cells will survive {} filtering.'.format(col2))
            print((1-cmasks3).sum(), 'cells will survive {} filtering.'.format(col3))
            print(cmask.sum(), 'cells will survive all filters.\n')
        
        else:
            cmask = cmasks1 & cmasks2 & cmasks3
            print(cmasks1.sum(), 'cells will survive {} filtering.'.format(col1))
            print(cmasks2.sum(), 'cells will survive {} filtering.'.format(col2))
            print(cmasks3.sum(), 'cells will survive {} filtering.'.format(col3))
            print(cmask.sum(), 'cells will survive all filters.\n')
                     
    elif flist1 and flist2:
        
        if set(flist1).issubset(set(adata.obs[col1])) == False: 
            raise Exception('flist1 is not a subset of {}'.format(col1))
        
        if set(flist2).issubset(set(adata.obs[col2])) == False: 
            raise Exception('flist2 is not a subset of {}'.format(col2))

        cmasks1 = adata.obs[col1].isin(flist1)
        cmasks2 = adata.obs[col2].isin(flist2)
        
        if isnot1 and isnot2:
            cmask = ~cmasks1 & ~cmasks2 
            print((1-cmasks1).sum(), 'cells will survive {} filtering.'.format(col1))
            print((1-cmasks2).sum(), 'cells will survive {} filtering.'.format(col2))
            print(cmask.sum(), 'cells will survive all filters.\n')
    
        elif isnot1:
            cmask = ~cmasks1 & cmasks2 
            print((1-cmasks1).sum(), 'cells will survive {} filtering.'.format(col1))
            print(cmasks2.sum(), 'cells will survive {} filtering.'.format(col2))
            print(cmask.sum(), 'cells will survive all filters.\n')
        
        elif isnot2:
            cmask = cmasks1 & ~cmasks2 
            print(cmasks1.sum(), 'cells will survive {} filtering.'.format(col1))
            print((1-cmasks2).sum(), 'cells will survive {} filtering.'.format(col2))
            print(cmask.sum(), 'cells will survive all filters.\n')
                  
        else:
            cmask = cmasks1 & cmasks2 
            print(cmasks1.sum(), 'cells will survive {} filtering.'.format(col1))
            print(cmasks2.sum(), 'cells will survive {} filtering.'.format(col2))
            print(cmask.sum(), 'cells will survive all filters.\n')
            
    elif flist1:
        
        if set(flist1).issubset(set(adata.obs[col1])) == False: 
            raise Exception('flist1 is not a subset of {}'.format(col1))
            
        if isnot1:
            cmask = ~adata.obs[col1].isin(flist1)
            print((cmask).sum(), 'cells will survive {} filtering.\n'.format(col1))
        
        else:
            cmask = adata.obs[col1].isin(flist1)
            print(cmask.sum(), 'cells will survive {} filtering.\n'.format(col1))

    # Return filtered adata
    if return_filtered: 
        
        print('Shape before filtering:', adata.shape)
        adata = adata[cmask].copy()
        
        if adata.shape[0] == cmask.sum():
            print('Shape after filtering:', adata.shape, '\n')
        else: 
            raise Exception('Filtered adata shape is not the same as cmask')
        
        if col1 and col2 and col3:
            
            if isnot1 and isnot2 and isnot3:
                
                if (set(adata.obs[col1].unique()) != set(flist1)) and (set(adata.obs[col2].unique()) != set(flist2)) and (set(adata.obs[col3].unique()) != set(flist3)):
                    print('{} unique values after filtering:'.format(col1), list(adata.obs[col1].unique()), '\n')
                    print('{} unique values after filtering:'.format(col2), list(adata.obs[col2].unique()), '\n')
                    print('{} unique values after filtering:'.format(col3), list(adata.obs[col3].unique()), '\n')
                else: 
                    raise Exception('Unique values are not the same as the values in the flists')
            
            elif isnot1 and isnot2:
                
                if (set(adata.obs[col1].unique()) != set(flist1)) and (set(adata.obs[col2].unique()) != set(flist2)) and (set(adata.obs[col3].unique()) == set(flist3)):
                    print('{} unique values after filtering:'.format(col1), list(adata.obs[col1].unique()), '\n')
                    print('{} unique values after filtering:'.format(col2), list(adata.obs[col2].unique()), '\n')
                    print('{} unique values after filtering:'.format(col3), list(adata.obs[col3].unique()), '\n')
                else: 
                    raise Exception('Unique values are not the same as the values in the flists')
                
            elif isnot1 and isnot3:
                
                if (set(adata.obs[col1].unique()) != set(flist1)) and (set(adata.obs[col2].unique()) == set(flist2)) and (set(adata.obs[col3].unique()) != set(flist3)):
                    print('{} unique values after filtering:'.format(col1), list(adata.obs[col1].unique()), '\n')
                    print('{} unique values after filtering:'.format(col2), list(adata.obs[col2].unique()), '\n')
                    print('{} unique values after filtering:'.format(col3), list(adata.obs[col3].unique()), '\n')
                else: 
                    raise Exception('Unique values are not the same as the values in the flists')
            
            elif isnot2 and isnot3:
                
                if (set(adata.obs[col1].unique()) == set(flist1)) and (set(adata.obs[col2].unique()) != set(flist2)) and (set(adata.obs[col3].unique()) != set(flist3)):
                    print('{} unique values after filtering:'.format(col1), list(adata.obs[col1].unique()), '\n')
                    print('{} unique values after filtering:'.format(col2), list(adata.obs[col2].unique()), '\n')
                    print('{} unique values after filtering:'.format(col3), list(adata.obs[col3].unique()), '\n')
                else: 
                    raise Exception('Unique values are not the same as the values in the flists')
                
            elif isnot1:
                
                if (set(adata.obs[col1].unique()) != set(flist1)) and (set(adata.obs[col2].unique()) == set(flist2)) and (set(adata.obs[col3].unique()) == set(flist3)):
                    print('{} unique values after filtering:'.format(col1), list(adata.obs[col1].unique()), '\n')
                    print('{} unique values after filtering:'.format(col2), list(adata.obs[col2].unique()), '\n')
                    print('{} unique values after filtering:'.format(col3), list(adata.obs[col3].unique()), '\n')
                else: 
                    raise Exception('Unique values are not the same as the values in the flists')
            
            elif isnot2:
                
                if (set(adata.obs[col1].unique()) == set(flist1)) and (set(adata.obs[col2].unique()) != set(flist2)) and (set(adata.obs[col3].unique()) == set(flist3)):
                    print('{} unique values after filtering:'.format(col1), list(adata.obs[col1].unique()), '\n')
                    print('{} unique values after filtering:'.format(col2), list(adata.obs[col2].unique()), '\n')
                    print('{} unique values after filtering:'.format(col3), list(adata.obs[col3].unique()), '\n')
                else: 
                    raise Exception('Unique values are not the same as the values in the flists')
            
            elif isnot3:
                
                if (set(adata.obs[col1].unique()) == set(flist1)) and (set(adata.obs[col2].unique()) == set(flist2)) and (set(adata.obs[col3].unique()) != set(flist3)):
                    print('{} unique values after filtering:'.format(col1), list(adata.obs[col1].unique()), '\n')
                    print('{} unique values after filtering:'.format(col2), list(adata.obs[col2].unique()), '\n')
                    print('{} unique values after filtering:'.format(col3), list(adata.obs[col3].unique()), '\n')
                else: 
                    raise Exception('Unique values are not the same as the values in the flists')
            
            else:
                
                if (set(adata.obs[col1].unique()) == set(flist1)) and (set(adata.obs[col2].unique()) == set(flist2)) and (set(adata.obs[col3].unique()) == set(flist3)):
                    print('{} unique values after filtering:'.format(col1), list(adata.obs[col1].unique()), '\n')
                    print('{} unique values after filtering:'.format(col2), list(adata.obs[col2].unique()), '\n')
                    print('{} unique values after filtering:'.format(col3), list(adata.obs[col3].unique()), '\n')
                else: 
                    raise Exception('Unique values are not the same as the values in the flists')
                  
        elif col1 and col2:
            
            if isnot1 and isnot2:
                
                if (set(adata.obs[col1].unique()) != set(flist1)) and (set(adata.obs[col2].unique()) != set(flist2)):
                    print('{} unique values after filtering:'.format(col1), list(adata.obs[col1].unique()), '\n')
                    print('{} unique values after filtering:'.format(col2), list(adata.obs[col2].unique()), '\n')
                else:
                    raise Exception('Unique values are not the same as the values in the flists')
        
            elif isnot1:
                
                if (set(adata.obs[col1].unique()) != set(flist1)) and (set(adata.obs[col2].unique()) == set(flist2)):
                    print('{} unique values after filtering:'.format(col1), list(adata.obs[col1].unique()), '\n')
                    print('{} unique values after filtering:'.format(col2), list(adata.obs[col2].unique()), '\n')
                else:
                    raise Exception('Unique values are not the same as the values in the flists')
            
            elif isnot2:
                
                if (set(adata.obs[col1].unique()) == set(flist1)) and (set(adata.obs[col2].unique()) != set(flist2)):
                    print('{} unique values after filtering:'.format(col1), list(adata.obs[col1].unique()), '\n')
                    print('{} unique values after filtering:'.format(col2), list(adata.obs[col2].unique()), '\n')
                else:
                    raise Exception('Unique values are not the same as the values in the flists')
            
            else:
                
                if (set(adata.obs[col1].unique()) == set(flist1)) and (set(adata.obs[col2].unique()) == set(flist2)):
                    print('{} unique values after filtering:'.format(col1), list(adata.obs[col1].unique()), '\n')
                    print('{} unique values after filtering:'.format(col2), list(adata.obs[col2].unique()), '\n')
                else:
                    raise Exception('Unique values are not the same as the values in the flists')
                                                            
        elif col1:
            
            if isnot1: 
                
                if (set(adata.obs[col1].unique()) != set(flist1)):
                    print('{} unique values after filtering:'.format(col1), list(adata.obs[col1].unique()), '\n')
                else: 
                    raise Exception('Unique values are not the same as the values in the flists')
        
            else:
                
                if (set(adata.obs[col1].unique()) == set(flist1)):
                    print('{} unique values after filtering:'.format(col1), list(adata.obs[col1].unique()), '\n')
                else: 
                    raise Exception('Unique values are not the same as the values in the flists')
                   
        return adata
    
    # Return cmask 
    else:
                  
        return cmask  
    
####### highest average

def highest_average(gene, tpm, group, second):
    if second:
        if (tpm.loc[:, gene].sort_values(ascending = False).index[0] == group) or (tpm.loc[:, gene].sort_values(ascending = False).index[1] == group): 
            return gene
        else:
            return 
    else: 
        if (tpm.loc[:, gene].sort_values(ascending = False).index[0] == group):
            return gene
        else:
            return

####### de testing 

def sc_de_testing(adata, groupby = '', reference = 'rest', method = 'wilcoxon', corr_method = 'benjamini-hochberg', key_added = '', 
                  pval_adj = 0.05, upregulated = False, downregulated = False, logFC = 2, high_average = True, second = True, save_excel = False, fname = ''):
    
    print('Started running sc.tl.rank_genes_groups...')

    sc.tl.rank_genes_groups(adata, groupby = groupby, reference='rest', 
                            use_raw = False, method = method, 
                            corr_method = corr_method, 
                            key_added = key_added, n_genes = len(adata.var_names))
    
    print('Finished running sc.tl.rank_genes_groups...assign to deDfDict\n')
    
    
    deDfDict = {}
    
    if high_average:
        
        tpm = adata_centroids(groupby, adata)
    
    for group in adata.obs[groupby].unique(): 
        
        print(groupby + ':', str(group))
        deDfDict[group] = pd.DataFrame([adata.uns[key_added]['scores'][group],
                                       adata.uns[key_added]['names'][group], 
                                       adata.uns[key_added]['logfoldchanges'][group],
                                       adata.uns[key_added]['pvals'][group],
                                       adata.uns[key_added]['pvals_adj'][group]], 
                                       columns = list(adata.uns[key_added]['names'][group]),
                                       index = ['scores', 'names', 'logfoldchanges', 'pvals', 'pvals_adj']
                                      ).T
        
        deDfDict[group] = deDfDict[group][deDfDict[group]['pvals_adj'] < pval_adj]
        
        if upregulated:
            deDfDict[group] = deDfDict[group][deDfDict[group]['logfoldchanges'] > np.log2(logFC)]
        
        if downregulated:
            deDfDict[group] = deDfDict[group][deDfDict[group]['logfoldchanges'] < np.log2(logFC)]
        
        # Both upregulated and downregulated genes
        else:
            deDfDict[group] = deDfDict[group][deDfDict[group]['logfoldchanges'].abs() > np.log2(logFC)]
            
        
        print(deDfDict[group]['logfoldchanges'].min())
        
        print('Number of DE genes:', deDfDict[group].shape[0])
        
        if high_average: 
            
            with multiprocessing.Pool() as p:
                genes = [gene for gene in p.starmap(highest_average, 
                                                    iterable = [*zip([gene for gene in deDfDict[group]['names']], 
                                                    itertools.repeat(tpm, len(deDfDict[group]['names'])), 
                                                    itertools.repeat(group, len(deDfDict[group]['names'])),
                                                    itertools.repeat(second, len(deDfDict[group]['names'])))]) if gene is not None]
            
            # We don't have to call p.join() because starmap blocks until the result is ready
            
            gmask = deDfDict[group]['names'].isin(genes)
            deDfDict[group] = deDfDict[group][gmask]
            
            if second: 
                print('First and second highest average number of DE genes:', deDfDict[group].shape[0], '\n')
            
            else: 
                print('Highest average number of DE genes:', deDfDict[group].shape[0], '\n')
        
        if save_excel:
            fname = fname + '_{}.xlsx'.format(group)
            deDfDict[group].to_excel(fname)
        
    return(deDfDict)

######## clustermap with fixed size 

def fixedWidthClusterMap(dataFrame, row_cluster = False, col_cluster = False, method = 'average', metric = 'euclidean', vmin = None, vmax = None,
                         dendrogram_ratio=(.05, .05), cmap = plt.cm.get_cmap('RdBu_r'), linewidths = .5, cbar_pos = None, cellSizePixels=50):
    # From: https://stackoverflow.com/questions/52806878/seaborn-clustermap-fixed-cell-size
    # Calulate the figure size, this gets us close, but not quite to the right place
    dpi = plt.rcParams['figure.dpi']
    marginWidth = plt.rcParams['figure.subplot.right']-plt.rcParams['figure.subplot.left']
    marginHeight = plt.rcParams['figure.subplot.top']-plt.rcParams['figure.subplot.bottom']
    Ny,Nx = dataFrame.shape
    figWidth = (Nx*cellSizePixels/dpi)/0.8/marginWidth
    figHeigh = (Ny*cellSizePixels/dpi)/0.8/marginHeight

    # do the actual plot
    grid = sns.clustermap(dataFrame, figsize=(figWidth, figHeigh), linewidths = linewidths, yticklabels = 1, xticklabels = 1, cmap = cmap, vmin = vmin, vmax = vmax,
                          row_cluster = row_cluster, col_cluster = col_cluster, method = method, metric = metric, dendrogram_ratio = dendrogram_ratio, cbar_pos = cbar_pos)

    # calculate the size of the heatmap axes
    axWidth = (Nx*cellSizePixels)/(figWidth*dpi)
    axHeight = (Ny*cellSizePixels)/(figHeigh*dpi)

    # resize heatmap
    ax_heatmap_orig_pos = grid.ax_heatmap.get_position()
    grid.ax_heatmap.set_position([ax_heatmap_orig_pos.x0, ax_heatmap_orig_pos.y0, 
                                  axWidth, axHeight])

    # resize dendrograms to match
    ax_row_orig_pos = grid.ax_row_dendrogram.get_position()
    grid.ax_row_dendrogram.set_position([ax_row_orig_pos.x0, ax_row_orig_pos.y0, 
                                         ax_row_orig_pos.width, axHeight])
    ax_col_orig_pos = grid.ax_col_dendrogram.get_position()
    grid.ax_col_dendrogram.set_position([ax_col_orig_pos.x0, ax_heatmap_orig_pos.y0+axHeight,
                                         axWidth, ax_col_orig_pos.height])
    return grid # return ClusterGrid object

######## colormap

def custom_colormap(colors,positions=[],cmap_name = 'my_cmap',register=False):
    """
    example of use:
            my_cmap = custom_colormap(['#000000','#f800f8',''#748c08'],[-100,2,50])
    
    input:
        colors: list of colors as hex codes
        positions: list of floats, same lengths as colors, indicating at what position
                    to place the pure color, if empty list, will space colors evenly.
                    The range can be all real numbers, will rescale to go from 0 to 1.
        
        register: if True, will "register" the colormap name, which can be useful for some applications
            
    output:
        colormap object.
    
    Would be nice to add:
    option to provide for color names, e.g. 'magenta', not just hex codes
    
    More info about making colormaps:
    https://matplotlib.org/examples/pylab_examples/custom_cmap.html
    
    """
    
    # make position range from 0 to 1:
    if len(positions)==len(colors):
        positions = np.array(positions).astype(float)
        positions = (positions-positions.min())
        positions = positions/positions.max()
    else:
        positions = np.linspace(0,1,len(colors))

    rgbs = []

    #turn hex into rgb,scale to max 1, append position
    for h,pos in zip(colors,positions):
        h = h.strip('#')
        rgb = np.array([int(h[i:i+2], 16) for i in (0, 2 ,4)])
        rgb = rgb/255.
        rgbs.append(list(rgb)+[pos])

    reds = []
    greens = []
    blues = []
    ll = []
    
    # prepare the color dictionary as described in
    # https://matplotlib.org/examples/pylab_examples/custom_cmap.html
    for nr,[r,g,b,pos] in enumerate(rgbs):
        for ll,col in zip([reds,greens,blues],[r,g,b]): #ll - list of lists
            #append the left position, the starting red value
            ll.append([pos,col,col])

    cdict = {}
    for key,value in zip(['red','green','blue'],[reds,greens,blues]):
        cdict[key] = tuple([tuple(i) for i in value])

    # make colormap
    cm = LinearSegmentedColormap(cmap_name, cdict)
    
    #register cmap:
    plt.register_cmap(cmap=cm)
    return cm

####### gene expression bar plot that shows replicate and within cluster single-cell dispersion: 
# Idea from Nico 
def bar_condition_rep_se(replicateSeriesDict, states = None, cmap_condition = None, cmap_state = None, ylab = None, title = None, 
                         edgecolor = None, space = 1.5, width = 0.75, h = 8, capsize = 3, se_scaling = 6, save = None): 
    """
    Barplot with bar = condition mean (when len(conditions) > 1) else bar = cell state mean within single condition. 
    In addition, each bar has scatterdots = biological replicate means, and errorbars = standard error biological replicates. 
    
    Parameters
    ----------
    replicateSeriesDict : defaultdict(dict) 
        Nested dictionary of form ``{'condition': {replicate: pd.Series}}``. 
        For each cell in the adata, pd.Series stores gene expression values for each cells. 
        Therefore, pd.Series is of shape (n_cells, ) with index = cell state assignment. 
        pd.Series can be generated as follows:
        
        replicateSeriesDict = defaultdict(lambda:defaultdict(dict))
        for condition in replicateDict.keys():
            for replicate in replicateDict[condition]: 
                cmask = adata.obs[replicate_col] == replicate
                replicateSeriesDict[condition][replicate] = pd.Series(adata[cmask, gene].X.todense().A1, index = adata[cmask, ].obs[state_col])
        
    states : list, default=None
        Ordered list of cell states to plot.
        
    space : int, default=1.5
        Space between states.
        
    width : int, default=0.75 
        Width of every condition within states. 
    
    h :  int,  default=8 
        Height of the axes 
    
    cmap_condition : dict, default=None
        Dictionary of form ``{condition: colorcode}``
    
    cmap_state : dict, default=None
        Dictionary of form ``{state: colorcode}`` (if len(conditions)==1)

    ylab : str, default=None
        y-axis label. 
    
    title : str, default=None
        Title. 
    
    edgecolor : str, default=None
        Edgecolor bars.
        
    save : str, default=None
        path + filename. 
    
    Returns
    -------
    Barplot figure. 
    
    """

    # Initialize statistic variables 
    replicate_state_mean = defaultdict(lambda:defaultdict(dict))
    replicate_state_se = defaultdict(lambda:defaultdict(dict))
    condition_state_mean = defaultdict(dict)

    for condition in replicateSeriesDict.keys():
   
        for replicate in replicateSeriesDict[condition].keys():
        
            for state in states:
                cmask = replicateSeriesDict[condition][replicate].index == state 
            
                # Calculate replicate mean for each state 
                replicate_state_mean[condition][replicate][state] = replicateSeriesDict[condition][replicate][cmask].mean()
            
                # Calculate replicate standard error for each state 
                replicate_state_se[condition][replicate][state] = replicateSeriesDict[condition][replicate][cmask].std()/np.sqrt(len(replicateSeriesDict[condition][replicate][cmask]))

    for condition in replicateSeriesDict.keys():
    
        # Create temporary dataframe to calculate condition mean from replicates 
        tmp_data = functools.reduce(pd.DataFrame.join, 
                                    [pd.DataFrame(replicateSeriesDict[condition][replicate], columns = [replicate]) for replicate in replicateSeriesDict[condition].keys()])
    
        for state in states: 
            cmask = tmp_data.index == state 
            # Calculate condition mean for each state 
            condition_state_mean[condition][state] = tmp_data[cmask].mean(0).mean()
        
    # Define spacing in array x 
    x = [0]
    
    n = len(replicateSeriesDict)
    for i,j in enumerate(itertools.cycle([space*width, n*width])):
        value = x[-1]+j
        x.append(value)
        if i == 2*len(states):
            break
    
    x = np.array(x)
    x = np.delete(x, 0)
    w = x.max()
    x = np.delete(x, -1)
    x = x[::2]

    # Initialize fig: 
    a,fig,gs=startfig(w, h)
    
    if (len(replicateSeriesDict.keys()) == 1):
        
        for i, (condition, value) in enumerate(replicateSeriesDict.items()):
            value = value.keys()
            
            # Plot mean bar for each condition 
            a.bar([pos + width*i for pos in x],
                  [condition_state_mean[condition][state] for state in condition_state_mean[condition].keys()],
                  width = width,
                  color = [cmap_state[state] for state in states], 
                  edgecolor = edgecolor, 
                  zorder = 4)
        
            # For each replicate, plot error bar 
            xoffsets = dict(zip(value, np.array(list(range(-(len(value) - 1)//2, 0)) + list(range(1, (len(value) - 1)//2 + 2))) \
            if len(value) % 2 == 0 else np.array(range(-(len(value) - 1)//2, (len(value) -1)//2+1))))
        
            for replicate in value:
                a.errorbar([pos + width*i for pos in x] + xoffsets[replicate]/se_scaling, 
                           y = [replicate_state_mean[condition][replicate][state] for state in replicate_state_mean[condition][replicate].keys()],
                           xerr = 0, yerr = [replicate_state_se[condition][replicate][state] for state in replicate_state_se[condition][replicate].keys()], 
                           fmt = 'o', ecolor = 'black', color = 'white', markeredgecolor = 'black', capsize = capsize, zorder = 6)
            
    else: 
        
        for i, (condition, value) in enumerate(replicateSeriesDict.items()):
            value = value.keys()
            
            # Plot meanline for each condition 
            a.bar([pos + width*i for pos in x],
                  [condition_state_mean[condition][state] for state in condition_state_mean[condition].keys()],
                  width = width,
                  color = cmap_condition[condition], 
                  edgecolor = edgecolor, 
                  zorder = 4)
        
            # For each replicate, plot error bar
            xoffsets = dict(zip(value, np.array(list(range(-(len(value) - 1)//2, 0)) + list(range(1, (len(value) - 1)//2 + 2))) \
            if len(value) % 2 == 0 else np.array(range(-(len(value) - 1)//2, (len(value) -1)//2+1))))
        
            for replicate in value:
                a.errorbar([pos + width*i for pos in x] + xoffsets[replicate]/se_scaling, 
                           y = [replicate_state_mean[condition][replicate][state] for state in replicate_state_mean[condition][replicate].keys()],
                           xerr = 0, yerr = [replicate_state_se[condition][replicate][state] for state in replicate_state_se[condition][replicate].keys()], 
                           fmt = 'o', ecolor = 'black', color = 'white', markeredgecolor = 'black', capsize = capsize, zorder = 6)
            
    a.set_ylim(bottom = -0.001)
    a.grid(which='major', axis='y', zorder = 0) 
    a.set_xticks([pos + 0.5*width*(n-1) for pos in x])
    
    if states:
        a.set_xticklabels(states, rotation=90)
    
    if title: 
        a.set_title(title)
        
    a.set_ylabel(ylab)
    showspines(a,bottom=True,left=True)
    gs.tight_layout(fig)
    if save: 
        plt.savefig(save)

        
######## barplot with replicate scatterdots for abundance barplots 
def bar_scatter_rep(replicateData, replicateDict, state = True, space = 1.5, width = 0.6, h = 8, dotsize = 18, cmap_condition = None, 
                  cmap_state = None, xTicksLabs = None, ylab = None, title = None, edgecolor = None, color_scatter = 'black', saveName = None):
    """
    Make a bar for every condition and state & a scatter for every replicate and state 
    
    Parameters
    ----------
    replicateData : dict 
        dict with keys = replicates and values = 1D numpy array that contains y in order of the states
        
        dict with keys = replicates and values = one value (states are merged)
        
    replicateDict : dict
        dict with keys = conditions and values = replicates corresponding to conditions 
    
    state : bool 
        whether to plot states or not (default = True)
        
    space : int
        space between states
        
    width : int 
        width of every condition within states 
    
    h :  int
        height of the axes 
    
    dotsize : int
        dotsize of the scatter
    
    cmap_condition : dict 
        dict with keys = condition and values = hexcode colors 
    
    cmap_state : dict 
        dict with keys = state and values = hexcode colors (only if toplot == True and len(replicateDict.keys())==1)
    
    xTicksLabs : list
        list with states in the order of the dfs in ReplicateDataDict (if states == True)
    
    ylab : str
        y-axis label 
    
    title : str
        title 
    
    edgecolor : str
        edgecolor 
        
    color_scatter : str
        color of the scatter points  
    
    saveName : str
        path + filename 
    
    Returns
    -------
    Figure with meanline for every condition and states & scatter for every replicate and states colored as condition
    
    """
    # Calculate the ymean from the replicateDataDict input 
    ymean = {}
    
    for condition, value in replicateDict.items():
        ymean[condition] = np.vstack([replicateData[replicate] for replicate in value]).T.mean(axis = 1)
    
    # Define spacing in array x 
    x = [0]
    
    if state: 
        n = len(replicateDict)
        for i,j in enumerate(itertools.cycle([space*width, n*width])):
            value = x[-1]+j
            x.append(value)
            if i == 2*len(ymean[list(replicateDict.keys())[0]]):
                break
    else: 
        n = 1 
        for i,j in enumerate(itertools.cycle([space*width, n*width])):
            value = x[-1]+j
            x.append(value)
            if i == 2*len(replicateDict):
                break
        
    x = np.array(x)
    x = np.delete(x, 0)
    w = x.max()
    x = np.delete(x, -1)
    x = x[::2]
    
    # Start figure 
    a,fig,gs=startfig(w, h)
    
    if state and (len(replicateDict.keys()) == 1):
        
        for i, (condition, value) in enumerate(replicateDict.items()):
            # Plot meanline for each condition 
            a.bar([pos + width*i for pos in x],
                  ymean[condition],
                  width = width,
                  color = [cmap_state[state] for state in toplot],  
                  zorder = 4)
        
            # Plot Scatter for each replicate 
            for replicate in value: 
                a.scatter([pos + width*i for pos in x], 
                          [float('nan') if pop==0 else pop for pop in replicateData[replicate]],
                           color = color_scatter, 
                           edgecolor = edgecolor, 
                           s = dotsize,
                           zorder = 8) 
                
    elif state:
        
        for i, (condition, value) in enumerate(replicateDict.items()):
            # Plot meanline for each condition 
            a.bar([pos + width*i for pos in x],
                  ymean[condition],
                  width = width,
                  color = cmap_condition[condition], 
                  zorder = 4)
        
            # Plot Scatter for each replicate 
            for replicate in value: 
                a.scatter([pos + width*i for pos in x], 
                          [float('nan') if pop==0 else pop for pop in replicateData[replicate]],
                           color = color_scatter, 
                           edgecolor = edgecolor, 
                           s = dotsize,
                           zorder = 8)   
    else:
        for i, (condition, value) in enumerate(replicateDict.items()):
            # Plot meanline for each condition 
            a.bar([pos + width*0 for pos in x][i],
                  ymean[condition],
                  width = width,
                  color = cmap_condition[condition], 
                  zorder = 4)
        

            
            # Plot Scatter for each replicate 
            for replicate in value: 
                a.scatter([pos + width*0 for pos in x][i], 
                           [float('nan') if pop==0 else pop for pop in replicateData[replicate]],
                           color = color_scatter, 
                           edgecolor = edgecolor, 
                           s = dotsize, 
                           zorder = 8
                        ) 
    
        
    a.set_ylim(bottom = -0.001)
    a.grid(which='major', axis='y', zorder = 0) 
    a.set_xticks([pos + 0.5*width*(n-1) for pos in x])
    
    if state:
        a.set_xticklabels(xTicksLabs, rotation=90)
    else:
        a.set_xticklabels(replicateDict.keys(), rotation=90)
    
    if title: 
        a.set_title(title)
        
    a.set_ylabel(ylab)
    showspines(a,bottom=True,left=True)
    gs.tight_layout(fig)
    if saveName: 
        plt.savefig(saveName)
        
######## fit sklearn classifiers
# Define classifier functions 

def fit_multinomialNb(X_train, y_train, **kwargs):
    """
    Fit multinomialNb model to training data. 
    
    Parameters
    ----------
    X_train
        Training vectors, where n_samples is the number of samples and n_features is the number of features.
        In the case of scRNA-seq X_train is often a pandas DataFrame with elements = the mean expression of genes j for meta groups cells i. 
        
    y_train 
        Target values/class labels. 
        
    Returns 
    -------
    clfMultiNb: multinomialNb model for the training data. 
    
    """
    from sklearn.naive_bayes import MultinomialNB
    clfMultiNb = MultinomialNB(alpha = 0)
    clfMultiNb.fit(X_train, y_train)

    return(clfMultiNb)

def fit_logReg(X_train, y_train, **kwargs):

    """
    Fit Logistic Regression model to training data. 
    
    Parameters
    ----------
    X_train
        Training vectors, where n_samples is the number of samples and n_features is the number of features.
        In the case of scRNA-seq X_train is often a pandas DataFrame with elements = the mean expression of genes j for meta groups cells i. 
        
    y_train 
        Target values/class labels. 
    
    **kwargs 
        c_param: int, default=1
            Regularization parameter. The strength of the regularization is inversely proportional to C. 
            Must be strictly positive.
        
        max_iter: int, default=100
            The maximum number of iterations to be run.
        
    Returns 
    -------
    clflogReg: Logistic Regression model for the training data. 
    """
    c_param = kwargs.get('c_param', 1)
    max_iter = kwargs.get('max_iter', 100)
    
    from sklearn.linear_model import LogisticRegression  
    clfLogReg = LogisticRegression(C = c_param, max_iter = max_iter)
    clfLogReg.fit(X_train, y_train)

    return(clfLogReg)

def fit_linearSVM(X_train, y_train, **kwargs): 
    """
    Fit Linear SVM model to training data. 
    
    Parameters
    ----------
    X_train: {array-like, sparse matrix} of shape(n_samples, n_features)
        Training vectors, where n_samples is the number of samples and n_features is the number of features.
        In the case of scRNA-seq X_train is often a pandas DataFrame with elements = the mean expression of genes j for meta groups cells i. 
        
    y_train 
        Target values/class labels. 
        
    **kwargs 
        c_param: int, default=1
            Regularization parameter. The strength of the regularization is inversely proportional to C. 
            Must be strictly positive.
        
        max_iter: int, default=1000
            The maximum number of iterations to be run.
        
    Returns 
    -------
    clflinearSVM: Linear SVM model for the training data. 
    """
    c_param = kwargs.get('c_param', 1)
    max_iter = kwargs.get('max_iter', 1000)

    from sklearn.svm import LinearSVC
    clflinearSVM = LinearSVC(C = c_param, max_iter = max_iter)
    clflinearSVM.fit(X_train, y_train)

    return(clflinearSVM)

def fit_MLP(X_train, y_train, **kwargs): 
    """
    Fit Multi-layer Perceptron model to training data. 
    
    Parameters
    ----------
    X_train: {array-like, sparse matrix} of shape(n_samples, n_features)
        Training vectors, where n_samples is the number of samples and n_features is the number of features.
        In the case of scRNA-seq X_train is often a pandas DataFrame with elements = the mean expression of genes j for meta groups cells i. 
        
    y_train 
        Target values/class labels. 
    
    **kwargs
        hidden_layer_sizes: tuple, default=(100,)
        alpha = int, default = 0.0001
        
    Returns 
    -------
    **clfMLP* Multi-layer Perceptron model for the training data. 
    """    
    hidden_layer_sizes = kwargs.get('hidden_layer_sizes', (100,))
    alpha = kwargs.get('alpha', 0.0001)
    
    from sklearn.neural_network import MLPClassifier
    clfMLP = MLPClassifier(hidden_layer_sizes = hidden_layer_sizes, alpha = alpha)
    clfMLP.fit(X_train, y_train)

    return(clfMLP) 

def fit_randomForest(X_train, y_train): 
    """
    Fit Random Forest model to training data. 
    
    Parameters
    ----------
    X_train: {array-like, sparse matrix} of shape(n_samples, n_features)
        Training vectors, where n_samples is the number of samples and n_features is the number of features.
        In the case of scRNA-seq X_train is often a pandas DataFrame with elements = the mean expression of genes j for meta groups cells i. 
        
    y_train 
        Target values/class labels. 
        
    Returns 
    -------
    **clfrandomForest* Random Forest model for the training data. 
    """    
    from sklearn.ensemble import RandomForestClassifier
    clfrandomForest = RandomForestClassifier()
    clfrandomForest.fit(X_train, y_train)

    return(clfrandomForest)


###### classify cells adata 

def adata_classifier(classifier, adata, centroids = None, predict = False, group_xtrain = 'state', group_xtest = 'state', log = True, threshold = 0.0, **kwargs): 
    """
    Sklearn anndata.AnnData Classification. 
    If a dictionary of anndata.AnnData matrices is given the matrices are compared in all combinations/
    directions using the cartesian product of all anndata.AnnData. Run this option when you want to calculate reciprocals for each
    dataset-dataset comparison
    
    If a single anndata.AnnData X_test matrix and a centroids X_train DataFrame is given the 
    anndata.AnnData matrix is classified by the single centroids DataFrame.
    
    Parameters
    ----------
    classifier: 
        Sklearn classifier to use. Supported are: 
            fit_LinearSVM : sklearn.svm.LinearSVC
            fit_multinomialNb : from sklearn.naive_bayes.MultinomialNB
            fit_randomForest : sklearn.ensemble.RandomForestClassifier
            fit_MLP : sklearn.neural_network.MLPClassifier
            fit_logReg : sklearn.linear_model.LogisticRegression  
            
    adata: anndata.AnnData OR dict 
        X_test anndata.AnnData that needs to be classified by one X_train dataset that contains the same features/
        genes as the centroids DataFrame. 
        OR
        X_test/X_train dictionary with keys = datasets to be compared in all combinations/directions (also self-self), 
        values = corresponding annotated.AnnData matrix of that dataset. 
    
    centroids: DataFrame, default=None
        X_train that contains the same features/genes as the anndata.AnnData matrix. 

    predict: bool, default=False
        Classify cells in one X_test with one X_train dataset (forced classification).
    
    group_xtrain: str, default='state'
        Calculate the mean expression of given group of cells from X_train if adata is a dictionary.  
        
    group_xtest: str, default='state'
        Calculate the mean probability of cells from X_test being a state from X_train given that they belong to a group of X_test.obs. 
    
    log: bool, default=True
        Whether to log transform the data before classification 
        
    threshold: float, default=0.0 
        Classify all cells in a group of X_test as a state i from X_train if if the mean probability of cells from that group being state i > threshold. 
        Setting threshold at 0 means don't use threshold. 
        
    **kwargs
        Keyword arguments for classifier.fit(): 
        
        LogReg:
        c_param: int, default=1
            Regularization parameter. The strength of the regularization is inversely proportional to C. 
            Must be strictly positive.
        
        max_iter: int, default=100
            The maximum number of iterations to be run.
        
        Linear SVM:
            c_param: int, default=1
                Regularization parameter. The strength of the regularization is inversely proportional to C. 
                Must be strictly positive.
        
            max_iter: int, default=1000
                The maximum number of iterations to be run.
 
        MLP:
            hidden_layer_sizes: tuple, default=(100,)
            alpha = int, default = 0.0001
        
    Returns
    -------
    If predict: 
        prediction: DataFrame shape(n_cells, ); index = index adata.obs. 
            DataFrame that contains predicted class label per cell i. 
    
    Elif group_xtest: 
        mean_pred_proba_group: dict
            Dictionary with keys = 'Xtrain-to-X_test', 
            value = pandas DataFrame with elements = mean probability of cells from X_test being a state j from X_train 
            given that they belong to a group i of X_test. 
    
    Elif group_xtest and treshold: 
        pred_proba_threshold_group: DataFrame shape (n_cells, ); index = index adata.obs.  
            DataFrame that contains predicted class label per cell i, or 'Unknown' if no class was predicted for a group of cells. 
            All cells from a group have the same classification. 
     
    Dependencies functions 
    ----------------------
    
    mm.adata_centroids
    mm.fit_logReg
    mm.fit_multinomialNb
    mm.fit_linearSVM
    mm.fit_MLP
    mm.fit_randomForest

    """ 
    
    if (isinstance(adata, dict) == False):
    
        if not (centroids.columns == adata.var_names).all(): 
            raise Exception('Error: features X_train are not the same as features X_test')
    
        classes = centroids.index.values 
        
        pseudo_xtrain = centroids.sum(1).mean()/1e4
        
        print('Started fitting...')
        if log: 
            # .mean() here because SCTransform data is not total counts normalized 
            pseudo_xtest = adata.X.sum(1).mean()/1e4
    
        # Fit the data
        if log: 
            clf = classifier(np.log2(centroids + pseudo_xtrain), classes,  **kwargs)
        
        else:
            clf = classifier(centroids + pseudo_xtrain, classes,  **kwargs)
        
        print('Started predicting...')
        # Predict 
        if predict:
            if log:
                prediction = pd.DataFrame(clf.predict(np.log2(adata.X.todense() + pseudo_xtest)), 
                                                              index = adata.obs.index, columns = ['prediction'])
            else:
                prediction = pd.DataFrame(clf.predict(adata.X.todense()), index = adata.obs.index, columns = ['prediction'])
                
            # For each cell i return a column with str predicted state names
            return prediction
        
        elif group_xtest and threshold:
            print('Threshold:', str(threshold))
            
            # Predict
            if hasattr(clf, 'predict_proba'):
                if log:
                    pred_proba = pd.DataFrame(clf.predict_proba(np.log2(adata.X.todense() + pseudo_xtest)), 
                                                                          index = adata.obs.index, columns = clf.classes_) 
                else:
                    pred_proba = pd.DataFrame(clf.predict_proba(adata.X), index = adata.obs.index, columns = clf.classes_)
                    
            else:                                             
                if log:                                       
                    pred_proba = pd.DataFrame(clf.decision_function(np.log2(adata.X.todense() + pseudo_xtest)), 
                                                                            index = adata.obs.index, columns = clf.classes_)
                else:
                    pred_proba = pd.DataFrame(clf.decision_function(adata.X.todense()), index = adata.obs.index, columns = clf.classes_)
                                                                                 
                pred_proba = (pred_proba - pred_proba.min()) / (pred_proba.max() - pred_proba.min())
            
            # Add group column from X_test 
            pred_proba = pred_proba.join(adata.obs[group_xtest])
    
            # Calculate E[p(state X_train | group X_test)] and transpose to make rows normalize to 1 
            mean_pred_proba_group = pred_proba.groupby(group_xtest).mean().T
            
            above_thresholdDict = dict(zip(mean_pred_proba_group.idxmax(0)[mean_pred_proba_group.max(0) > threshold].index, 
                                           mean_pred_proba_group.idxmax(0)[mean_pred_proba_group.max(0) > threshold].values))
            
            pred_proba_threshold_group = pd.DataFrame([above_thresholdDict[state] if state in above_thresholdDict.keys() else 'Unknown' for state in adata.obs[group_xtest]],
                                                       index = adata.obs.index, columns = ['prediction'])
            
            return pred_proba_threshold_group
        
        elif group_xtest:
            
            # Predict
            if hasattr(clf, 'predict_proba'):
                if log:
                    pred_proba = pd.DataFrame(clf.predict_proba(np.log2(adata.X.todense() + pseudo_xtest)), 
                                                                          index = adata.obs.index, columns = clf.classes_)    
                else:
                    pred_proba = pd.DataFrame(clf.predict_proba(adata.X), index = adata.obs.index, columns = clf.classes_)
                                                    
            else:                                             
                if log:                                       
                    pred_proba = pd.DataFrame(clf.decision_function(np.log2(adata.X.todense() + pseudo_xtest)), 
                                                                            index = adata.obs.index, columns = clf.classes_)
                else:
                    pred_proba = pd.DataFrame(clf.decision_function(adata.X.todense()), index = adata.obs.index, columns = clf.classes_)
                                                                                 
                pred_proba = (pred_proba - pred_proba.min()) / (pred_proba.max() - pred_proba.min())
            
            # Add group column from X_test 
            pred_proba = pred_proba.join(adata.obs[group_xtest])
    
            # Calculate E[p(state X_train | group X_test)] and transpose to make rows normalize to 1 
            mean_pred_proba_group = pred_proba.groupby(group_xtest).mean().T
                
            # Return probability of cells from X_test being a state j from X_train given that they belong to a group i of X_test
            return mean_pred_proba_group
            
    else:
        # adata = dict of datasets 
        centroids = {}
        pred_proba = {}
        mean_pred_proba_group = {}
        
        # Calculate the training centroids: 
        for dataset in adata.keys():
            centroids[dataset] = adata_centroids(group_xtrain, adata[dataset])
            
        for X_train, X_test in itertools.product(adata.keys(), repeat=2):
            name = X_train + '_to_' + X_test
            print('Started running ' + name)
    
            if not (centroids[X_train].columns == adata[X_test].var_names).all(): 
                raise Exception('Error: features adata X_train are not the same as features adata X_test')
    
            classes = centroids[X_train].index.values 
            
            pseudo_xtrain = adata[X_train].X.sum(1).mean()/1e4
            
            if log: 
                pseudo_xtest = adata[X_test].X.sum(1).mean()/1e4
    
            # Fit the data 
            if log:
                clf = classifier(np.log2(centroids[X_train] + pseudo_xtrain), classes,  **kwargs)
            
            else:
                clf = classifier(centroids[X_train], classes, **kwargs)
                
            # Predict
            if hasattr(clf, 'predict_proba'):
                if log:
                    pred_proba[name] = pd.DataFrame(clf.predict_proba(np.log2(adata[X_test].X.todense() + pseudo_xtest)), 
                                                                                index = adata[X_test].obs.index, columns = clf.classes_)
                else:
                    pred_proba[name] = pd.DataFrame(clf.predict_proba(adata[X_test].X), index = adata[X_test].obs.index, columns = clf.classes_)
            
            else:
                if log:
                    pred_proba[name] = pd.DataFrame(clf.decision_function(np.log2(adata[X_test].X.todense() + pseudo_xtest)), 
                                                                                  index = adata[X_test].obs.index, columns = clf.classes_)
                else:
                    pred_proba[name] = pd.DataFrame(clf.decision_function(adata[X_test].X.todense()), index = adata[X_test].obs.index, columns = clf.classes_)
                    
                pred_proba[name] = (pred_proba[name] - pred_proba[name].min()) / (pred_proba[name].max() - pred_proba[name].min())
                
            # Add group column from X_test 
            pred_proba[name] = pred_proba[name].join(adata[X_test].obs[group_xtest])
    
            # Calculate E[p(state X_train | group X_test)] and transpose to make rows normalize to 1 
            mean_pred_proba_group[name] = pred_proba[name].groupby(group_xtest).mean().T

        # Return probability of cells from X_test being a state j from X_train given that they belong to a group i of X_test
        return mean_pred_proba_group
    
####### calculate reciprocals using output from classify cells adata 

# Define a calculate reciprocals function 
def calculate_reciprocals(mean_predp_Dict):
    """
    Calculate reciprocals from a mean_predp_* Dictionary. 
    
    Parameters
    ----------
    mean_predp_Dict 
        Dictionary with keys = Xtrain-to-X_test, 
        value = pandas DataFrame with matrix elements = 
        mean probability of cells from X_test being a state i from X_train 
        given that they belong to a group j of X_test.  
   
    Returns
    -------
    **reciprocals** Dictionary with keys = datasetOne_datasetTwo comparison, 
                    values = pandas DataFrame with matrix elements reciprocals of state_datasetOne-state_datasetTwo 
    
    """
    
    
    tmpDict = mean_predp_Dict.copy()
    reciprocals = {}

    for classificationOne, classificationTwo in itertools.product(mean_predp_Dict.keys(), repeat = 2):
    
        # check whether the cartesian product of the set of all classification (keys) equals the reverse of a classification: 
        if classificationOne.split('_')[2] == classificationTwo.split('_')[0]:
        
            if classificationOne.split('_')[0] == classificationTwo.split('_')[2]:
                print('Multiplying', classificationOne, 'with', classificationTwo)
            
                # Generate dict keys for the return reciprocal dict: 
                name = str(classificationOne.split('_')[0]) + '_' + str(classificationTwo.split('_')[0])
                nameT = str(classificationTwo.split('_')[0]) + '_' + str(classificationOne.split('_')[0])
            
                # Calculate the self-self classification reciprocals 
                # Check whether the cartesian product of the set of all classification (keys) equals a self-self classification:
                if (classificationOne.split('_')[0] == classificationTwo.split('_')[0]) and (classificationOne.split('_')[2] == classificationTwo.split('_')[2]):
        
                    # Multiply the self-self classifications with it's transpose that has column and index orders fixed to calculate the reciprocal:  
                    reciprocals[name] = tmpDict[classificationOne] * tmpDict[classificationOne].T.loc[tmpDict[classificationOne].index, 
                                                                                                      tmpDict[classificationOne].columns]
        
                else:
                    # Calculate the non-self-self reciprocals 
                
                    # Check whether we should take the transpose of the classificationTwo matrix: 
                    if not set(tmpDict[classificationOne].columns) == set(tmpDict[classificationTwo].columns):
                    
                        tmpDict[classificationTwo] = tmpDict[classificationTwo].T
               
                    # Check that the tranpose worked otherwise exit function:  
                    if not set(tmpDict[classificationOne].columns) == set(tmpDict[classificationTwo].columns):
                        raise Exception('Transpose did not work')
                       
                    # Check whether index and columns orders are of classificationOne matrix are equal to those of the transpose of the classificationTwo matrix:  
                    indexOrder = tmpDict[classificationOne].index == tmpDict[classificationTwo].index 
                    columnOrder = tmpDict[classificationOne].columns == tmpDict[classificationTwo].columns     
                    if all(indexOrder) == False or all(columnOrder) == False:
                    
                        # Assign the classificationTwo matrix with column and index re-ordered to equal those of the classificationOne matrix: 
                        tmpDict[classificationTwo] = tmpDict[classificationTwo].loc[tmpDict[classificationOne].index, 
                                                                                    tmpDict[classificationOne].columns]
               
                    # Verify that the reordering was succesful by re-checking the column and index orders: 
                    indexOrder = tmpDict[classificationOne].index == tmpDict[classificationTwo].index
                    columnOrder = tmpDict[classificationOne].columns == tmpDict[classificationTwo].columns
                    
                    # Check that re-ordering was succesful otherwise exit function: 
                    if all(indexOrder) == False or all(columnOrder) == False:
                        raise Exception('Reordering did not work')
                 
                    # Multiply the non-self-self classifications to calculate the reciprocal: 
                    reciprocals[name] = tmpDict[classificationOne] * tmpDict[classificationTwo]
                
                    # Take the transpose of each non-self-self reciprocal to fill the heatmap
                    reciprocals[nameT] =  reciprocals[name].T   
            
                # Apparently, I am not allowed to turn the base.index into astype('categorical'); therefore, I turn the columns into base.index
                reciprocals[name].columns = reciprocals[name].columns.astype(str)
                reciprocals[name].index = reciprocals[name].index.astype(str)
                
                reciprocals[nameT].columns = reciprocals[nameT].columns.astype(str)
                reciprocals[nameT].index = reciprocals[nameT].index.astype(str)
    
                print('Assigned', name) 
                print('Assigned', nameT)
    
    return reciprocals 

###### fix mirroring of SPRING plot 

def mirror_spring(path_to_coordinates):
    with open(path_to_coordinates,'r') as infile:
      
        for line in infile:
            line = line.split(',')
            out.append(','.join([line[0], line[1], str(-float(line[2]))]))
        
    with open(path_to_coordinates, 'w') as outfile: 
        outfile.write('\n'.join(out))
        

####### color dataframe cells for saving heatmaps to excel 

# From Rapolas Zilionis
def color_dataframe_cells(
    frame,
    cmap = mpl.cm.get_cmap('RdBu_r'),
    vmin = None,
    vmax = None,
    ):
    
        
    """colors cells of dataframe by their values
    input:
        - frame: pandas dataframe
        - cmap: colormap to use, e.g mpl.cm.get_cmap('RdBu_r')
        - vmin: min value to saturate colormap at
        - vmax: max value to saturate colormap at
        
    output: pandas "styler" object with cells colored. This "styler" object  does not have the
    full functionality of a pandas dataframe.
    
    Example of use (including saving to excel):
    color_dataframe_cells(my_dataframe).to_excel('table.xlsx')
    
    Depends on custom function value_to_color
    
    """
    
    if vmin is None:
        vmin = frame.min().min()
    if vmax is None:
        vmax = frame.max().max()
        
    return frame.style.applymap(lambda x: 'background-color: %s'%value_to_color(x,vmin,vmax,cmap=cmap))
    
###### start matplotlib figure #####
# From Rapolas Zilionis 

def startfig(w=4,h=2,rows=1,columns=1,wrs=None,hrs=None,frameon=True,return_first_ax=True):

    '''
    for initiating figures, w and h in centimeters
    example of use:
    a,fig,gs = startfig(w=10,h=2.2,rows=1,columns=3,wr=[4,50,1],hrs=None,frameon=True)
    hrs - height ratios
    wrs - width ratios
    frameon - whether first axes with frame
    
    returns:
    if return_first_ax=True
    a,fig,gs
    else
    fig,gs
    '''
    
    ratio = 0.393701 #1 cm in inch
    myfigsize = (w*ratio,h*ratio)
    fig = plt.figure(figsize = (myfigsize))
    gs = mpl.gridspec.GridSpec(rows, columns ,width_ratios=wrs,height_ratios=hrs)
    if return_first_ax==True:
        a = fig.add_subplot(gs[0,0],frameon=frameon)
        return a,fig,gs
    else:
        return fig,gs

######## showspines #########
# From Rapolas Zilionis 

def showspines(an_axes,top=False,right=False,bottom=False,left=False):
    """
    for specifying which spines to make visible in a plot.
    input: 
        an_axes - matplotlib axes object
    returns: nothing

    """
    #after reinstalling conda, top and left switches places...
    [i for i in an_axes.spines.items()][3][1].set_visible(top) #top line
    [i for i in an_axes.spines.items()][1][1].set_visible(right) #right line
    [i for i in an_axes.spines.items()][2][1].set_visible(bottom) #bottom line
    [i for i in an_axes.spines.items()][0][1].set_visible(left) #left line
    an_axes.tick_params(bottom=bottom,right=right,left=left,top=top)
    
    
####### value_to_color #########

def value_to_color(value,vmin,vmax,cmap=mpl.cm.get_cmap('RdBu_r'),string_color='#FFFFFF'):
    
    """takes a value (float or int) and turns a hex code (string).
    input:
        - value: float or integer
        - vmin: min value in full range
        - vmax: max value in full range
        - cmap: colormap
        - string_color: what color to return is string given as "value",
          default is white (#FFFFFF)
        
        
    output:
        - hex code (string)
    """
    if type(value)==str:
        return string_color
    
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    rgb = cmap(norm(float(value)))
    return "#{:02x}{:02x}{:02x}".format(int(rgb[0]*255),int(rgb[1]*255),int(rgb[2]*255))


####### adata cell percent abudance #########
def calculate_percent(adata, replicate_col = '', state_col = '', dummy_col = '', states = []):
    count = adata.obs.loc[:, [replicate_col, state_col, dummy_col]].groupby([replicate_col, state_col]).count()
    count.columns = ['count']
    percent_dict = {}
    for replicate in count.index.get_level_values(replicate_col).unique():
        percent = count.loc[replicate]/count.loc[replicate].sum()*100.
        percent.columns = ['percentage']
        percent_dict[replicate] = percent
    
    # Create a dict with keys = replicates and values = df reindexed in the right order and with NA = 0 
    percent_dict_ordered = {}
    for replicate,value in percent_dict.items():
        percent_dict_ordered[replicate] = value.reindex(states).fillna(0)

    # Convert value df to 1D numpy array 
    percent_dict_ordered_flat = {}
    for replicate in percent_dict_ordered:
        percent_dict_ordered_flat[replicate] = percent_dict_ordered[replicate].values.flatten() 
        
    return(percent_dict_ordered_flat)


####### view color dict ########

def view_color_dict(color_dict, size =1, save = ''):
    import matplotlib.ticker as ticker
    pal = color_dict.values()
    pal_names = color_dict.keys()
    n = len(pal)
    size = 1
    f, ax = plt.subplots(1, 1, figsize=(n * size, size))
    ax.imshow(np.arange(n).reshape(1, n),
              cmap=mpl.colors.ListedColormap(list(pal)),
              interpolation="nearest", aspect="auto")
    ax.set_xticks(np.arange(n))
    ax.set_yticks([-.5, .5])
    # Ensure nice border between colors
    ax.set_xticklabels([name for name in pal_names], rotation=90)
    # The proper way to set no ticks
    ax.yaxis.set_major_locator(ticker.NullLocator())
    if save: 
        plt.savefig(save)