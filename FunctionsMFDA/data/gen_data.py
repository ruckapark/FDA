# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 17:08:26 2023

Create four clusters of noisy data  with three dimensions for each cluster

These clusters should be in normalised and non-normalised form

@author: Admin
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import signal

if __name__ == '__main__':
    
    X = pd.DataFrame(index = range(200),columns = range(20))
    Y = pd.DataFrame(index = range(200),columns = range(20))
    Z = pd.DataFrame(index = range(200),columns = range(20))
    
    #%% Dimension 1 - impulse signal
    
    #generate signal with three categories
    #sys = [a,b,c#
    #impulse: a*x''(t) + b*x'(t) + c*x(t)
    
    for i in range(10):
        
        sys = ([1.0],[1.0 + np.random.random()/10,2.0 + np.random.random()/10,1.0])
        T = np.linspace(0,7,200)
        t,y = signal.impulse(sys,T = T)
        X[i] = y + np.random.normal(scale = 0.02,size = len(y))
        
    for i in range(5):
        
        sys = ([1.0],[1.0 + np.random.random()/10,2.0 + np.random.random()/10,3.0])
        T = np.linspace(0,7,200)
        t,y = signal.impulse(sys,T = T)
        X[i+10] = y + np.random.normal(scale = 0.02,size = len(y))
        
    for i in range(5):
        
        sys = ([1.5],[1.0 + np.random.random()/10,2.0 + np.random.random()/10,2.0])
        T = np.linspace(0,7,200)
        t,y = signal.impulse(sys,T = T)
        X[i+15] = y + np.random.normal(scale = 0.02,size = len(y))
        
    xplot,xaxe = plt.subplots(4,5,sharex = True,sharey = True, figsize = (20,12))
    for i in range(20):
        xaxe[i//5,i%5].plot(X[i])
    plt.suptitle('X dimension time series')
    
    #%% Dimension 2 - 3rd order polynomial
    
    def poly3(a,b,c,d,t):
        return np.array(a*t**3 + b*t**2 + c*t + d*np.ones(len(t)),dtype = float)
    
    for i in range(5):
        t = np.linspace(0,2,200)
        a,b,c,d = 0.3,0.5+(np.random.random()-1)/10,0.7,-0.3+(np.random.random()-1)/10
        Y[i] = poly3(a,b,c,d,t) + np.random.normal(scale = np.random.random()/30,size = len(t))
    
    for i in range(5):
        t = np.linspace(0,2,200)
        a,b,c,d = 0.3,-0.5+(np.random.random()-1)/10,0.7,0.2+(np.random.random()-1)/10
        Y[i+5] = poly3(a,b,c,d,t) + np.random.normal(scale = np.random.random()/30,size = len(t))
        
    for i in range(10):
        t = np.linspace(0,2,200)
        a,b,c,d = -0.3,0.5+(np.random.random()-1)/10,-0.7,0.7+(np.random.random()-1)/10
        Y[i+10] = poly3(a,b,c,d,t) + np.random.normal(scale = np.random.random()/30,size = len(t))
        
    yplot,yaxe = plt.subplots(4,5,sharex = True,sharey = True, figsize = (20,12))
    for i in range(20):
        yaxe[i//5,i%5].plot(Y[i])
    plt.suptitle('Y dimension time series')
    
    #%% Dimension 3 - gaussian pulse
    
    #generate gauss pulse - two variants
    for i in range(10):
        t = np.linspace(-1,1,200)
        Z[i] = signal.gausspulse(t, fc=1.5+(np.random.random()-0.5)) + np.random.normal(scale = 0.02,size = len(t))
    
    for i in range(10):
        t = np.linspace(-1,1,200)
        Z[i+10] = signal.gausspulse(t, fc=10+(2*np.random.random()-1)) + np.random.normal(scale = 0.02,size = len(t))
    
    zplot,zaxe = plt.subplots(4,5,sharex = True,sharey = True, figsize = (20,12))
    for i in range(20):
        zaxe[i//5,i%5].plot(Z[i])
    plt.suptitle('Z dimension time series')
    
    #%% Save data
    X.to_csv('xdata_raw.csv',header = False,index = False)
    Y.to_csv('ydata_raw.csv',header = False,index = False)
    Z.to_csv('zdata_raw.csv',header = False,index = False)
    
    #%% Standardise
    def norm_df(df):
        return (df - df.min().min())/(df.max().max() - df.min().min())
    Xn = norm_df(X)
    Yn = norm_df(Y)
    Zn = norm_df(Z)
    
    #%% Save standardised data
    
    Xn.to_csv('xdata_norm.csv',header = False,index = False)
    Yn.to_csv('ydata_norm.csv',header = False,index = False)
    Zn.to_csv('zdata_norm.csv',header = False,index = False)