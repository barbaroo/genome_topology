# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 16:16:24 2022

@author: 31649
"""

import numpy as np



import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import networkx as nx
import matplotlib.pyplot as plt

#from read_HiC import name_chromosomes
import seaborn as sns

import sys
path=r'C:\Users\31649\Documents\genome analysis\genome_topology\functions'
sys.path.append(path)

#from plotting_tools import set_layout
from genome_topology import normalize_psc
from genome_topology import get_matrix
#from genome_topology import fractal_dimension
from genome_topology import make_graph

#Load data

resolution = '80'
quantile_thresh = 0.99
path = 'data/zoomify processed'
samples = ['Control1', 'Control2', 'Treated1', 'Treated2']
#control1 = pd.DataFrame()
#control2 = pd.DataFrame()
#treated1 = pd.DataFrame()
#treated2 = pd.DataFrame()

frames = [None]*4

for ind, sample in enumerate(samples):
    df =  pd.read_csv(f'{path}/{sample}{resolution}', sep = '\t', header = None, 
                      names = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'count', 'balanced'],
                     dtype={'chrom1': 'str'})
    
    
    frames[ind] = df.dropna().reset_index(drop=True)
    sample_col = [sample] * len(frames[ind])
    frames[ind]['Sample'] = sample_col
    
control1 = frames[0]
control2 = frames[1]
treated1 = frames[2]
treated2 = frames[3]   
    
    
#Find threshold and filter data
type_counts = 'count'
contacts = pd.DataFrame()
total_data =( list(control1[type_counts]) + list(control2[type_counts]) 
             + list(treated1[type_counts]) + list(treated2[type_counts]))
data = {'counts': total_data}
data = pd.DataFrame(data)
threshold = data['counts'].quantile(quantile_thresh)

for frame in frames:
    frame = frame[frame[type_counts]>= threshold]
    frame = pd.DataFrame(frame)
    frames = [contacts, frame]
    contacts = pd.concat(frames)

    
new_index = np.linspace(1, len(contacts), len(contacts), dtype = int)
contacts['Index'] = new_index
contacts = contacts.set_index('Index')    
    
#Set parameters for binning
res = float(resolution)
n_start_bins = 150
start = n_start_bins* res * 1000
increment = 50* res * 1000
max_bin = np.max(control1['end1'])
cutoffs = np.arange(start, max_bin , increment, dtype = int)


#calculate cumulative CT parameters
chromosomes = control1['chrom1'].unique()
path_df = 'results/cumulative'
#samples = ['Control1', 'Control2', 'Treated1', 'Treated2']
#samples = ['Control1']





for sample in samples:
    
    N_contacts=np.zeros(len(chromosomes))
    clustering = np.zeros((len(cutoffs), len(chromosomes)))
    P=np.zeros((len(cutoffs), len(chromosomes)))
    S=np.zeros((len(cutoffs), len(chromosomes)))
    X=np.zeros((len(cutoffs), len(chromosomes)))
    Dim_fractal = np.zeros(len(chromosomes))
    r2_fractalfit = np.zeros(len(chromosomes))
   
    
    
    
    contacts_sample = contacts[contacts['Sample'] == sample]
    
    for t, chrom in enumerate(chromosomes):
        print(chrom)
        contacts_chr= contacts_sample[(contacts_sample['chrom1']== chrom) & (
            contacts_sample['chrom2']==chrom)]
        
        for j, cutoff in enumerate(cutoffs):
            print(cutoff)
            data_sel = contacts_chr[(contacts_chr['end1']<= cutoff)&(
                contacts_chr['end2']<= cutoff)]
        
            length1 = (data_sel['end1'] -  data_sel['start1'])//2
            data_sel['position1']=  pd.Series(data_sel['start1'] + length1)
    
            length2 = (data_sel['end2'] -  data_sel['start2'])//2
            data_sel['position2']=  pd.Series(data_sel['start2'] + length2)
    
            index = [data_sel['position1'], data_sel['position2']]
            counts = np.array(data_sel['count'])
            index=np.array(index)
            index= np.transpose(index)

            N_contacts[t]= len(index)
    
            mat, psc = get_matrix(index, chrom)
            #write_topology_matrix(mat, pathfile, namefile = f'{chrom}_{selection_sample}_{resolution}')
            P[j, t], S[j, t], X[j, t]=normalize_psc(psc, N_contacts[t])
    
            #Dim_fractal[t], r2_fractalfit[t] = fractal_dimension(mat, plot_fig=0)
    
            G=make_graph(index)
            clustering[j, t] = nx.average_clustering(G)
            
        df = {'Cut-off': cutoffs, 'Parallel (%)': P[:,t], 'Series (%)': S[:,t], 
              'Cross (%)': X[:,t], 'Clustering coeff': clustering[:,t]}
          #, 'Clustering coeff': clustering,
          #    'Fractal dimension': Dim_fractal, 'R2': r2_fractalfit,
          #    'N contacts': N_contacts}
        df = pd.DataFrame(df)

        df.to_feather(f'{path_df}/cumulative_{sample}_{chrom}_{resolution}_{quantile_thresh}.feather ')
    