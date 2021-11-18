# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 15:16:05 2021

@author: scalvinib
"""
#functions for topological analysis of genome

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx





#read the pdb
def open_pdb(ID_chr, path):
    file_string= '{}/chrom{}.pdb'.format(path, ID_chr)
    #print(file_string)
    pdb = open(file_string, "r")
    ind=0
    for line in pdb:
        if ind < 0:
            print(line)
        else:
            break
        ind= ind+1 
        
    return pdb 


#retrieve coordinates of the chosen chromosome
def select_chrom(chrom, path):
    chr_vec=np.array(['chr a', 'chr b', 'chr c','chr d', 'chr e', 'chr f', 
    'chr g', 'chr h', 'chr i', 'chr j', 'chr k', 'chr l','chr m', 'chr n',
    'chr o', 'chr p', 'chr q', 'chr r', 'chr s', 'chr t'])      
    ID_chr=chr_vec[chrom][4]
    pdb=open_pdb(ID_chr, path)
    counter=0
    coord=[]
    for line in pdb:        
        if (line[17:22] != chr_vec[chrom]):
             break
        x=float(line[31:38])
        y=float(line[39:46])
        z=float(line[47:54])
        coord.append([x,y,z])
        counter= counter +1
    pdb.close()     
    coord=np.array(coord)
    
    
    return counter, coord 


#finds contacts in the chromosome based on geometrical distance
def geom_distance(x_frame, r_cutoff, nb):
    mat_dist=np.zeros((x_frame.shape[0], x_frame.shape[0]))
    indexes = []
    n_contact=0
    for t1 in range (x_frame.shape[0]):
        for t2 in range(t1+nb+1,x_frame.shape[0]):
            x=(x_frame[t2,0]-x_frame[t1,0])**2
            y=(x_frame[t2,1]-x_frame[t1,1])**2
            z=(x_frame[t2,2]-x_frame[t1,2])**2
            dist=x+y+z
            if (dist<= r_cutoff**2):
                n_contact=n_contact+1
                mat_dist[t1,t2]=1.0
                mat_dist[t2,t1]=1.0
                indexes.append([t1,t2])     
            
    return mat_dist, n_contact, np.array(indexes)



#retrieve topology matrix and topological fractions
def topology_matrix(index_loops, neighbours):
    #print(index_loops)
    S=0
    P=0
    Pp=0
    C=0
    CS=0
    CP=0
    top_mat=np.zeros((len(index_loops), len(index_loops)))
    
    for t in range(len(index_loops)):
        i=index_loops[t,0]       
        j=index_loops[t,1]
        
        
        for t2 in range(t+neighbours+1, len(index_loops)):
     
            k=index_loops[t2,0]
            l=index_loops[t2,1]
            
            #series
            if (j <k):
                S=S+1
                top_mat[t, t2]=1.0
                top_mat[t2, t]=1.0
                
            #parallel    
            elif (i>k and j<l):
                P=P+1
                top_mat[t, t2]=2.0
                top_mat[t2, t]=3.0
            
            #5: CP
            #6: CP-1    
            elif (i==k and j<l):
                top_mat[t, t2]=5.0
                top_mat[t2, t]=6.0
                CP=CP+1
           
            elif (i==k and l<j):
                top_mat[t, t2]=6.0
                top_mat[t2, t]=5.0
                CP=CP+1
               
            elif (k>i and j==l):
                top_mat[t,t2]=6.0
                top_mat[t2,t]=5.0
                CP=CP+1
                
            elif(i>k and l==j):
                top_mat[t,t2]=5.0
                top_mat[t2,t]=6.0
                CP=CP+1
                             
            #inverse parallel
            elif (k>i and l<j):
                Pp=Pp+1
                top_mat[t, t2]=3.0
                top_mat[t2, t]=2.0
                
            elif (j ==k):
                top_mat[t,t2]=7
                top_mat[t2,t]=7
                CS=CS+1

                
                
            if (k>i and k<j and j<l):
                C=C+1
                top_mat[t, t2]=4.0
                top_mat[t2, t]=4.0
            elif (i>k and i< l and j> l):
                C=C+1
                top_mat[t, t2]=4.0
                top_mat[t2, t]=4.0
                
              
    return S, P, Pp, C,CP, CS, top_mat  


def normalize_psc(psc, n_cont):
    size_mat= n_cont*n_cont-n_cont
    half_mat=size_mat/2
    p = psc[1]/half_mat
    s = psc[2]/half_mat
    x = psc[3]/half_mat
    return p,s,x

#retrieve topology matrix, Duane's version
def get_matrix(index,protid):
    
    if index.shape == (0,):
        print('Error - index empty')
        mat = np.zeros((len(index), len(index)),dtype = 'int')
        psc = [protid,0,0,0]
        return mat,psc

    if np.shape(index)[1] == 2:

        #create a numerical and character matrix based on the amount of nonzero values found in the previous function
        mat = np.zeros((len(index), len(index)),dtype = 'int')

        #Change the values based on the type of connection
        

        P = 0
        S = 0
        X = 0

        for x in range(0,len(index)):
            i = index[x,0]
            j = index[x,1]
            for y in range(x+1,len(index)):
                k = index[y,0]
                l = index[y,1]
                #series
                if (j < k):
                    S=S+1
                    mat[x, y]=1
                    mat[y, x]=1

                #parallel    
                elif (i>k and j<l):
                    P=P+1
                    mat[x, y]=2
                    mat[y, x]=3
                
                #5: CP
                #6: CP-1    
                elif (i==k and j<l):
                    mat[x, y]=5
                    mat[y, x]=6
                    P += 1
            
                elif (i==k and l<j):
                    mat[x, y]=6
                    mat[y, x]=5
                    P += 1
                
                elif (k>i and j==l):
                    mat[x,y]=6
                    mat[y,x]=5
                    P += 1
                    
                elif(i>k and l==j):
                    mat[x,y]=5
                    mat[y,x]=6
                    P += 1
                #inverse parallel
                elif (k>i and l<j):
                    P += 1
                    mat[x, y]=3
                    mat[y, x]=2
                #CS
                elif (j ==k):
                    mat[x,y]=7
                    mat[y,x]=7
                    S += 1
                #Cross
                if (k>i and k<j and j<l):
                    X += 1
                    mat[x, y]=4
                    mat[y, x]=4
                elif (i>k and i< l and j> l):
                    X += 1
                    mat[x, y]=4
                    mat[y, x]=4

        psc = [protid,P,S,X]

        return mat,psc

    elif np.shape(index)[1] == 4:

        mat = np.zeros((len(index), len(index)),dtype = 'int')
    
        P = 0
        S = 0
        X = 0
        I = 0
        T = 0
        L = 0

        for x in range(len(index)):
            chain1 = False
            
            i = index[x][0]
            j = index[x][1]
            chaini = index[x][2]
            chainj = index[x][3]
            
            if chaini == chainj:
                chain1 = True
                
            for y in range(x+1,len(index)):
                chain2 = False
                
                k = index[y][0]
                l = index[y][1]
                chaink = index[y][2]
                chainl = index[y][3]
                
                set1 = set([chaini,chainj])
                set2 = set([chaink,chainl])

                if chaink == chainl:
                    chain2 = True
                    
                if chain1 and chain2:
                    if chaini == chaink:
                        #series
                        if j < k:
                            S += 1
                            mat[x,y] = 2
                            mat[y,x] = 2
                        #parallel
                        elif k < i and j < l:
                            P += 1
                            mat[x,y] = 1
                            mat[y,x] = 1
                            
                        elif i < k and l < j:
                            P += 1
                            mat[x,y] = 1
                            mat[y,x] = 1
                            
                        elif (i==k and j<l):
                            mat[x, y]=1
                            mat[y, x]=1
                            P += 1
                
                        elif (i==k and l<j):
                            mat[x, y]=1
                            mat[y, x]=1
                            P += 1

                        elif (k>i and j==l):
                            mat[x,y]=1
                            mat[y,x]=1
                            P += 1

                        elif(i>k and l==j):
                            mat[x,y]=1
                            mat[y,x]=1
                            P += 1
                        #CS
                        elif j == k:
                            S += 1
                            mat[x,y] = 2
                            mat[y,x] = 2
                        #Cross
                        if (k>i and k<j and j<l):
                            X += 1
                            mat[x, y]=3
                            mat[y, x]=3
                        elif (i>k and i< l and j> l):
                            X += 1
                            mat[x, y]=3
                            mat[y, x]=3
                        #Independent
                    else:
                        I += 1
                        mat[x,y] = 4
                        mat[y,x] = 4

                #I - multiple chains
                elif not set1.intersection(set2):
                    I += 1
                    mat[x,y] = 4
                    mat[y,x] = 4
                
                #T
                elif chain1 and set1.intersection(set2) :
                    T += 1
                    mat[x,y] = 5
                    mat[y,x] = 5
                elif chain2 and set1.intersection(set2):
                    T += 1
                    mat[x,y] = 5
                    mat[y,x] = 5
                elif ~chain1 and ~chain2 and set1.intersection(set2) and set1 != set2:
                    T += 1
                    mat[x,y] = 5
                    mat[y,x] = 5
                #L
                elif ~chain1 and ~chain2 and set1 == set2:
                    L += 1
                    mat[x,y] = 6
                    mat[y,x] = 6
                else:
                    print('error - ',i,chaini,j,chainj,k,chaink,l,chainl)
                
        stats = [protid,P,S,X,I,T,L]
        return mat,stats    


#Build network
def make_graph(index):
    G = nx.Graph()
    array_of_tuples = map(tuple, index)
    tuple_of_tuples =tuple(array_of_tuples)
    G.add_edges_from(tuple_of_tuples)
    return G


#Box counting function for box counting method
# From https://github.com/rougier/numpy-100 (#87)
def boxcount(Z, k):
    S = np.add.reduceat(np.add.reduceat(Z, np.arange(0, Z.shape[0], k), axis=0),
    np.arange(0, Z.shape[1], k), axis=1)

    # We count non-empty (0) and non-full boxes (k*k)
    return len(np.where((S > 0) & (S < k*k))[0])



#Calculate fractal dimension with box counting method
def fractal_dimension(Z, plot_fig = 0):
    threshold=np.mean(Z)
    # Only for 2d image
    assert(len(Z.shape) == 2)


    # Transform Z into a binary array
    Z = (Z > threshold)

    # Minimal dimension of image
    p = min(Z.shape)

    # Greatest power of 2 less than or equal to p
    n = 2**np.floor(np.log(p)/np.log(2))

    # Extract the exponent
    n = int(np.log(n)/np.log(2))

    # Build successive box sizes (from 2**n down to 2**1)
    sizes = 2**np.arange(n, 1, -1)

    # Actual box counting with decreasing size
    counts = []
    for size in sizes:
        counts.append(boxcount(Z, size))

    # Fit the successive log(sizes) with log (counts)

    
    coeffs = np.polyfit(np.log(sizes), np.log(counts), 1)
    rsquared=r_squared(np.log(sizes), -np.log(counts), 1)
    D="%.3f" % round(coeffs[0], 3)
    
    
    #Plot fit
    if plot_fig:
       plt.figure()
       plt.scatter(np.log(sizes), -np.log(counts))
       plt.plot(np.log(sizes),-np.log(sizes)*coeffs[0]-coeffs[1], color='m', label='D = {}'.format(D))
       plt.xlabel('log(S)')
       plt.ylabel('-log(N)')
       plt.title('BOX COUNTING METHOD')
       plt.legend()

    return -coeffs[0], rsquared


#evaluates fractal dimension fit
def r_squared(x, y, degree):
    results = {}

    coeffs = np.polyfit(x, y, degree)

     # Polynomial Coefficients
    results['polynomial'] = coeffs.tolist()

    # r-squared
    p = np.poly1d(coeffs)
    # fit values, and mean
    yhat = p(x)                         # or [p(z) for z in x]
    ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    results['determination'] = ssreg / sstot
    #plt.plot(x,y)
    #plt.plot(x, (x*coeffs[0]+coeffs[1]))
    return results['determination']


