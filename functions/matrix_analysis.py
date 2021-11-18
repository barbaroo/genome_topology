# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 11:29:02 2021

@author: scalvinib
"""

import numpy as np
import PIL.Image
import matplotlib.pyplot as plt
from plotting_tools import save_figures

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