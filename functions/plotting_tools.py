# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 14:32:40 2021

@author: scalvinib
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
import warnings

def set_layout(SMALL_SIZE = 13, MEDIUM_SIZE = 14, BIGGER_SIZE = 20):


    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)
    
    
def save_figures(plot,path, name_file, method=None):
    if method:
        name_figure = '{}/{}_{}'.format(path, method, name_file)
    else:
        name_figure = '{}/{}'.format(path, name_file)
        
    plot.savefig('{}.jpg'.format(name_figure), bbox_inches='tight')
    plot.savefig('{}.eps'.format(name_figure), bbox_inches='tight')
    
    
def change_width(ax, new_value) :
    for patch in ax.patches :
        current_width = patch.get_width()
        diff = current_width - new_value

        # we change the bar width
        patch.set_width(new_value)

        # we recenter the bar
        patch.set_x(patch.get_x() + diff * .5)
        
        
        
def plot_vertical_line(threshold,y):
    thresh_line_y=np.linspace(-0.5,2,len(y))
    thresh_line_x=np.ones(len(y))*threshold
    return thresh_line_x,thresh_line_y



def matrix_plot(mat,protid):

    #create custom colormap
    newcolors = np.array([[218/255, 219/255, 228/255,1], #grey - 
                      [131/255, 139/255, 197/255,1],    #purple S
                      [172/255,200/255,247/255,1],       #blue P
                      [174/255,213/255,129/255,1],         #mint green P-1
                      [186/255, 155/255, 201/255,1],        #red purple - X
                      [172/255, 200/255, 247/255,1],        #blue P
                      [174/255, 213/255, 129/255,1],        #mint green
                      [131/255,139/255, 197/255,1]])        #purple S
    newcmp = ListedColormap(newcolors)

    fig, ax = plt.subplots()
    color = plt.get_cmap(newcmp, 8)

    #plot data
    pngmat = plt.imshow(mat,cmap=color,vmin = np.min(mat)-.5, vmax = np.max(mat)+.5)
    ax.set_title(protid)
    ax.tick_params(labelleft = False,labelbottom = False,bottom = False,left= False)
    cbar = fig.colorbar(pngmat)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cbar.ax.set_yticklabels(['-','S','P','P-1','X','CP','CP-1','CS'])