# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 14:32:40 2021

@author: scalvinib
"""
import matplotlib.pyplot as plt

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