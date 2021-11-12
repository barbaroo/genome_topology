# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 16:05:29 2021

@author: scalvinib
"""
import scipy as sy
import numpy as np

def lin_fit(x,y):    
    guess1= 0.0
    guess2=0.0
    guess=[guess1, guess2]

    errfunc2 = lambda p, Series_two, y: ( line(x, *p) - y)**2
    optim, success= sy.optimize.leastsq(errfunc2, guess[:], args=( x, y))
    fit=line(x, optim[0] , optim[1])
    return fit

def line(x, a , b):
    y= a*x+b
    return y

    
def exponential(x,A,B,C):
    y=B*(1-np.exp(-C*(x-A)))
    return y    
