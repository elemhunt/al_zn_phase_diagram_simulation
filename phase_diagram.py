#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 21:14:53 2019

@author: elemhunt
"""

import matplotlib.pyplot as plt
import numpy as np
import Gibbs_Module as GM
from scipy.spatial import ConvexHull

#%% Creating x and G Matrices


T = np.linspace(300,1000.,1000)


x = np.linspace(0.001,0.999,1000)

XVals1 = []
XVals2 = []
XVals3 = []
XVals4 = []
XVals5 = []
XVals6 = []


for i in range(len(T)):
    
  
    GLiq =[]
    GFcc = []
    GHcp =[]    
    for j in range(len(x)):
        GLiq.append(GM.GlobGLiq(x[j],T[i])) #blue line
        GFcc.append(GM.GlobGFcc(x[j],T[i])) #orange line
        GHcp.append(GM.GlobGHcp(x[j],T[i])) #Green line
        
    xcoor = np.hstack((x,x,x))
    allG = np.hstack((GLiq,GFcc,GHcp))
    Allpoints = np.column_stack((xcoor,allG))
        
    hull = ConvexHull(Allpoints)
   
    
    vert = hull.vertices
    jumps = np.zeros(len(vert))
    
    
    
#%% Max 1    
    
    for k in range(1,len(vert)):
        jump =abs(vert[k]-vert[k-1]) 
        jumps[k] = jump
    
    for l in range(len(vert)):
        for m in ([1,2,3]):
            if jumps[l] == len(x)*m:
                jumps[l] = 0
    
    Max = np.argmax(jumps)
    
#%%    Max #2
    
    jumps2 = np.zeros(len(vert))
    
    for n in range(Max):
        max2 = []
        jumps2[n] = jumps[n]
        Max2 = np.argmax(jumps2)
        if (jumps[Max]> jumps[Max2]):
            max2.append(Max2)

#%%    Compiling X values
            
    X1 = xcoor[vert[Max]] 
    X2 = xcoor[vert[Max-1]]
    XVals1.append(X1)  #BLUE LINE
    XVals2.append(X2)  #ORANGE LINE
    
    
    for q in range(len(max2)):   
        X3 = xcoor[vert[max2[q]]]
        X4 =xcoor[vert[max2[q]-1]]
        XVals3.append(X3) #GREEN LINE
        XVals4.append(X4) #RED LINE
        
#%% Ploting

        
plt.plot(XVals1, T,'.')
plt.plot(XVals2, T, '.')
plt.plot(XVals3, T,'.')
plt.plot(XVals4, T,'.')

#plt.grid(True)
plt.ylabel('Temperature (K)')
plt.xlabel('Atomic Zinc Concentration')
plt.title('Equilibrium Phase Diagram of Aluminum-Zinc')