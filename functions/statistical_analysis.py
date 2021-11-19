# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 14:15:57 2021

@author: scalvinib
"""


import scipy.stats as ss
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy

def statistical_comparison(dist1, dist2):
    shapiro1=ss.shapiro(dist1)
    shapiro2=ss.shapiro(dist2)
    levene=ss.levene(dist1, dist2)
    
    normality = (np.array([shapiro1[1], shapiro2[1]]) > 0.05).all()
    equal_variance = levene[1]> 0.05
    
    if normality:
        print('The distributions are normal: {}, {}'.format(shapiro1[1], shapiro2[1]))
        
        if equal_variance:
            print('The distributions have equal variance: {}'.format(levene[1]))
            comparison_test = ss.ttest_ind(dist1, dist2)
            comparison = comparison_test[1] > 0.05
            string_test='Ttest'
            print('{}: {} = {}'.format(string_test, comparison, comparison_test[1]))
            
        else:
            print('The distributions do not have equal variance: {}'.format(levene[1]))
            comparison_test = ss.ttest_ind(dist1, dist2, equal_var = False)
            comparison = comparison_test[1] > 0.05
            string_test='Welch test'
            print('{}: {} = {}'.format(string_test, comparison, comparison_test[1]))
    else:
        print('the distributions are not normal: {}, {}'.format(shapiro1[1], shapiro2[1]))
        
        if equal_variance:
            print('The distributions have equal variance: {}'.format(levene[1]))
            comparison_test = ss.mannwhitneyu(dist1, dist2)
            comparison = comparison_test[1] > 0.05
            string_test='Mannwhitneyu test'
            print('{}: {} = {}'.format(string_test, comparison, comparison_test[1]))
            
        else:
            print('The distributions do not have equal variance: {}'.format(levene[1]))
            comparison_test = ss.kstest(dist1, dist2)
            comparison = comparison_test[1] > 0.05
            string_test='Kolmogorov test'
            print('{}: {} = {}'.format(string_test, comparison, comparison_test[1]))
    
    return shapiro1[1], shapiro2[1], levene[1], comparison_test[1], string_test



def comparison_distributions(dist1, dist2, path_results=None):
    test_results=statistical_comparison(dist1, dist2)
    round_test_results=[]
    for result in test_results[:-1]:
        result_round = "%.3f" % round(result, 3)
        
        if (result_round == '0.000'):
            result_round= "%.1E" % result
            
        round_test_results.append(result_round)
        
    if path_results:
        df = {'Shapiro1': [test_results[0]], 'Shapiro2': [test_results[1]], 'Levene': [test_results[2]], 
         '{}'.format(test_results[4]) : [test_results[3]]}
        df=pd.DataFrame(df)
        df.to_csv('{}_statistics.csv'.format(path_results))
        
    return round_test_results


def correlate (var1, var2):
    var1=np.asarray(var1.astype(np.float16))
    var2=np.asarray(var2.astype(np.float16))
    corr=ss.pearsonr(var1,var2)
    
    corr_coeff="%.2f" % round(corr[0], 2)
    p_coeff= "%.3f" % round(corr[1], 3)
    if (p_coeff=='0.000'):
        p_coeff= "%.1E" % corr[1]
    
    return corr_coeff, p_coeff

def find_peaks_distribution(length,n_peaks, threshold=0.00000):
    plt.figure()
    ax = sns.kdeplot(length)
    x = ax.lines[0].get_xdata() # Get the x data of the distribution
    y = ax.lines[0].get_ydata() # Get the y data of the distribution
    peak= scipy.signal.find_peaks(y, threshold)
    array=np.array(peak[0])
    for t in range(len(array)):
        plt.scatter(x[array[t]],y[array[t]])
    return x,array