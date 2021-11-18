# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 16:38:42 2021

@author: scalvinib
"""

import numpy as np


def name_chromosomes(chrom_number):
    name_chromosomes=['chr{}'.format(num) for num in np.arange(1,chrom_number,1)]
    chrx = 'chrX'
    chrx = list(chrx.split(" "))
    name_chromosomes += chrx
    return name_chromosomes


def MiddlePointLoci(contacts_chr):
    diffA = (contacts_chr['end_A']-contacts_chr['start_A'])//2
    diffB = (contacts_chr['end_b']-contacts_chr['start_B'])//2
    
    contacts_chr['index A']=contacts_chr['start_A']+ diffA
    contacts_chr['index B']=contacts_chr['start_B']+ diffB
    
    return contacts_chr

def Map_Topology_to_Structure(contacts_chr, bins, mode = 'Sum'):
    
    trace_structural=np.zeros(len(bins))
    
    if (mode == 'Mean'):
        
        for t in range(len(bins)-1):
            select_chr= contacts_chr[(contacts_chr['index A']< bins[t+1]) & (contacts_chr['index A']> bins[t])]
        
            if (len(select_chr['Matrix trace']) == 0):
                trace_structural[t]= 0
            elif(len(select_chr['Matrix trace']) == 0):
                trace_structural[t]=select_chr['Matrix trace'][0]
            else:
                trace_structural[t]=np.mean(select_chr['Matrix trace'])
            
    if (mode == 'Sum'):
        
        for t in range(len(bins)-1):
            select_chr= contacts_chr[(contacts_chr['index A']< bins[t+1]) & (contacts_chr['index A']> bins[t])]
            trace_structural[t]=np.sum(select_chr['Matrix trace'])
      
    return trace_structural