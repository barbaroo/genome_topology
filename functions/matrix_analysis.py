# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 11:29:02 2021

@author: scalvinib
"""

import numpy as np
import PIL.Image
import matplotlib.pyplot as plt
from plotting_tools import save_figures
import pandas as pd

def open_matrix(file_name):
    im = PIL.Image.open(file_name)
    imarray = np.array(im)
    mat=np.copy(imarray)
    return mat


def PlotMatrix_SelectFraction(select_relation, matrix, path_savefig = False):
    mat_parallel = Select_entangled_fraction(matrix, mode= select_relation)

    plot=plt.figure()
    plt.title('{} ONLY'.format(select_relation.upper()));
    plt.imshow(mat_parallel, cmap='Blues');
    
    if path_savefig:
        save_figures(plot,path_savefig, name_file = 'Matrix_only', method = select_relation)
    
def Select_entangled_fraction(mat, mode= None):
    if (mode == 'Parallel'):
        mat[mat == 1] = 0
        mat[mat == 7] = 0
        mat[mat == 4] = 0
        mat[mat != 0] = 1
        
    elif (mode == 'Cross'):
        mat[mat == 1] = 0
        mat[mat == 7] = 0
        mat[mat == 2] = 0
        mat[mat == 3] = 0
        mat[mat == 5] = 0
        mat[mat == 6] = 0
        mat[mat != 0] = 1 
                
    else:   
        mat[mat == 1] = 0
        mat[mat == 7] = 0
        mat[mat != 0] = 1 
        
    return mat  


def length_L_pattern(imarray):    
    y=imarray.shape[0]
    x=imarray.shape[1]
    length_pattern=np.zeros(y-1)

    for i in range(y-1):
        step=0
        for j in range(i+1,x):
            if (imarray[i,j]==1):
                step=step+1        
            length_pattern[i]=step
    return length_pattern


def write_topology_matrix(matrix, pathfile, namefile = 'top_matrix' ,format_output = 'feather'):
     
    size = matrix.shape[0]
    mat_zeros = matrix - np.ones((size, size))
    mat_upper = np.triu(mat_zeros)
    
    #get nonzero indexes for non-series matrix elements
    non_zero_indexes = np.nonzero(mat_upper)
    non_zero_indexes = np.transpose(non_zero_indexes)
    
    #find non zero (non-series) matrix elements
    mat_values = np.zeros(len(non_zero_indexes))
    mat_values = matrix[non_zero_indexes[:,0], non_zero_indexes[:,1]]
    
    #save sparse representation in dataframe
    sparse_matrix = {'Index1': non_zero_indexes[:,0], 'Index2': non_zero_indexes[:,1], 'Values': mat_values}
    sparse_matrix = pd.DataFrame(sparse_matrix)
    
    #Print dataframe to file
    
    if (format_output == 'csv'):
        sparse_matrix.to_csv('{}/{}.csv'.format(pathfile, namefile))
    elif (format_output == 'feather'):
        sparse_matrix.to_feather('{}/{}.feather'.format(pathfile, namefile))
        
    return 'matrix printed to file'



def read_topology_matrix(file):
    
    #open file
    if file.endswith('.feather'):
        df = pd.read_feather(file)
    elif file.endswith('.csv'):
        df = pd.read_csv(file)
        
    size = df['Index1'].max() + 1


    #compose matrix and fill in 2D array
    new_matrix = np.ones((size, size))
    
    for row in df.iterrows():
        x = int(row[1]['Index1'])
        y = int(row[1]['Index2'])
        value = row[1]['Values']
        new_matrix[x, y] = value
    
        if value == 3.:
            new_matrix[y, x] = 2.
        elif value == 5.:
            new_matrix[y, x] = 6.
    
        elif value == 6.:
            new_matrix[y, x] = 5.
        else:
            new_matrix[y, x] = value
    
    return new_matrix