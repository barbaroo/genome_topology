
"""
Spyder Editor

This is a temporary script file.
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
from genome_topology import fractal_dimension
from genome_topology import make_graph


#functions

def write_topology_matrix(matrix, pathfile, namefile = 'top_matrix' ,
                          format_output = 'feather'):
     
    size = mat.shape[0]
    mat_zeros = mat - np.ones((size, size))
    mat_upper = np.triu(mat_zeros)
    
    #get nonzero indexes for non-series matrix elements
    non_zero_indexes = np.nonzero(mat_upper)
    non_zero_indexes = np.transpose(non_zero_indexes)
    
    #find non zero (non-series) matrix elements
    mat_values = np.zeros(len(non_zero_indexes))
    mat_values = mat[non_zero_indexes[:,0], non_zero_indexes[:,1]]
    
    #save sparse representation in dataframe
    sparse_matrix = {'Index1': non_zero_indexes[:,0],
                     'Index2': non_zero_indexes[:,1], 'Values': mat_values}
    sparse_matrix = pd.DataFrame(sparse_matrix)
    
    #Print dataframe to file
    
    if (format_output == 'csv'):
        sparse_matrix.to_csv('{}/{}.csv'.format(pathfile, namefile))
    elif (format_output == 'feather'):
        sparse_matrix.to_feather('{}/{}.feather'.format(pathfile, namefile))
        
    return 'matrix printed to file'


#resolutions = ['20', '40','80', '160']
resolutions = ['20']

#['5120','2560','1280','640', '320','160','80','40', '20']
               #, '40', '80', '160', '320', '640', '1280', '2560', '5120',
               #'10240']

quantile_thresh = 0.988
path = 'data/zoomify processed'
path_df = 'results'
pathfile = 'results/matrices'

samples = ['Control1', 'Control2', 'Treated1', 'Treated2']
#samples = ['Control2']

plotting = 0

for resolution in resolutions:
    
    frames = [None]*4
    
    for ind, sample in enumerate(samples):
        print(resolution)
        df =  pd.read_csv(f'{path}/{sample}{resolution}', sep = '\t', header = None, 
              names = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2',
                       'count', 'balanced'], dtype={'chrom1': 'str'})
        print('here')
        frames[ind] = df.dropna().reset_index(drop=True)
        sample_col = [sample]* len(frames[ind])
        frames[ind]['Sample'] = sample_col


    control1 = frames[0]
    control2 = frames[1]
    treated1 = frames[2]
    treated2 = frames[3]
    
    #chromosomes = control1['chrom1'].unique()
    chromosomes = '1'
    
    
    type_counts = 'balanced'
    frames = [control1, control2, treated1, treated2]
    contacts = pd.DataFrame()
    total_data = (list(control1[type_counts]) + list(control2[type_counts]) 
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


    #Plot selected data
    if plotting:
        sns.histplot(x = 'balanced', hue = 'Sample', data = contacts, 
                     log_scale = True, bins= 100)
        plt.close()
        
        
    
    for selection_sample in samples:
        
        df = pd.DataFrame()
        
        selected_contacts = contacts[contacts['Sample'] == selection_sample]
        
        N_contacts=np.zeros(len(chromosomes))
        P=np.zeros(len(chromosomes))
        S=np.zeros(len(chromosomes))
        X=np.zeros(len(chromosomes))
        Dim_fractal = np.zeros(len(chromosomes))
        r2_fractalfit = np.zeros(len(chromosomes))
        clustering = np.zeros(len(chromosomes))

        for t, chrom in enumerate(chromosomes):
            print(chrom)
            contacts_chr= selected_contacts[(selected_contacts['chrom1'
                            ]== chrom) & (selected_contacts['chrom2']==chrom)]
    
            length1 = (contacts_chr['end1'] -  contacts_chr['start1'])//2
            position1 = pd.Series(contacts_chr['start1'] + length1)
            contacts_chr['position1']= position1
    
            length2 = (contacts_chr['end2'] -  contacts_chr['start2'])//2
            position2 = pd.Series(contacts_chr['start2'] + length2)
            contacts_chr['position2']=  position2
    
            index = [contacts_chr['position1'], contacts_chr['position2']]
            index=np.array(index)
            index= np.transpose(index)

            N_contacts[t]= len(index)
    
            mat, psc = get_matrix(index, chrom)
            #write_topology_matrix(mat, pathfile, 
             #           namefile = f'{chrom}_{selection_sample}_{resolution}')
            P[t], S[t], X[t]=normalize_psc(psc, N_contacts[t])
    
            Dim_fractal[t], r2_fractalfit[t] = fractal_dimension(mat, plot_fig=0)
    
            G=make_graph(index)
            clustering[t] = nx.average_clustering(G)
    
            if plotting:
                plt.figure(figsize=(5, 5))
                plt.imshow(mat[:300, :300])
                plt.title('{}'.format(chrom))
                
        plt.close()        
        #Save data        
        df = {'Chromosome': chromosomes, 'Parallel (%)': P, 'Series (%)': S, 
              'Cross (%)': X, 'Clustering coeff': clustering,
              'Fractal dimension': Dim_fractal, 'R2': r2_fractalfit,
              'N contacts': N_contacts}
        df = pd.DataFrame(df)

        df.to_feather(
            f'{path_df}/processed_{selection_sample}_{resolution}_{quantile_thresh}.feather ')
        