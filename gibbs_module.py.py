#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 16:09:17 2019
@author: elemhunt
"""
'''
Gibbs Energy Equations for Liquid, FCC and HCP states
'''


#%%
def GlobGLiq(xB,T):
    import math
    #xB is the conc Zn
    #xA is the conc of Al
    L0 = 10465.55-3.39259*T
    R = 8.314
    
    #Al Liq
    if (298.15 < T <= 700): 
        GAlLiq = 3028.895 + (125.251188*T) - (24.3671976*T*math.log(T)) - ((1.884662E-3)*T**(2)) - ((0.877664E-6)*T**(3)) + (74092.*T**(-1)) + ((79.34E-21)*T**(7))
        
    elif (700 < T <= 933.473):
        GAlLiq = -271.194 + (211.206596*T) - (38.5844296*T*math.log(T)) + ((18.531982E-3)*T**(2)) - ((5.764227E-6)*T**(3)) + (74092.*T**(-1)) + ((79.34E-21)*T**(7))
        
    elif (933.473 < T <= 2900):
        GAlLiq = -795.991 + (177.430209*T) - (31.748192*T*math.log(T))
        
    # Zn Liq
    if (298.15 < T <= 692.677):
        GZnLiq = -128.565 + (108.177019*T) - (23.701314*T*math.log(T)) - ((1.712034E-3)*T**(2)) - ((1.264963E-6)*T**(3)) -((358.949E-21)*T**(7))
        
    elif (692.677 < T <= 1700):
        GZnLiq = -3620.385 + (161.60854*T) - (31.38*T*math.log(T))    
        
    Func = ((1.-xB)*GAlLiq)+(xB*GZnLiq) + (R*T)*((1.-xB)*math.log(1.-xB)+xB*math.log(xB)) + ((1.-xB)*xB)*(L0)
    return Func
#%%    
        




def GlobGFcc(xB,T):
    import math
    #xB is the conc Zn
   #xA is the conc of Al
    L0 = 7297.48+0.47512*T 
    L1 = 6612.88-4.59110*T
    L2 = -3097.19+3.30635*T
    R = 8.314
    # Al fcc
    if (298.15 < T <= 700):
        GAlFfc = -7976.15 + (137.093038*T) - (24.3671976*T*math.log(T)) - ((1.884662E-3)*T**(2)) - ((0.877664E-6)*T**(3)) + (74092.*T**(-1))
        
    elif(700 < T <= 933.473):
        GAlFfc = -11276.24 + (223.048446*T) - (38.5844296*T*math.log(T)) + ((18.531982E-3)*T**(2)) - ((5.764227E-6)*T**(3)) + (74092.*T**(-1))
        
    elif(933.473 < T <= 2900):
        GAlFfc = -11278.361 + (188.684136*T) - (31.748192*T*math.log(T)) - ((1230.622E25)*T**(-9))
        
    # Zn fcc
    if (298.15 < T <= 692.677):
        GZnFfc = -4315.967 + (116.900389*T) - (23.701314*T*math.log(T)) - ((1.712034E-3)*T**(2)) - ((1.264963E-6)*T**(3))
        
    elif (692.677 < T <= 1700):
        GZnFfc = -8100.726 + (170.775964*T) - (31.38*T*math.log(T)) + ((470.47E24)*T**(-9))

    Func = ((1.-xB)*GAlFfc)+(xB*GZnFfc) + (R*T)*((1.-xB)*math.log(1.-xB)+xB*math.log(xB)) + ((1.-xB)*xB)*(L0 + L1*((1.-xB)-xB) + L2*((1.-xB)-xB)**2)
    return Func
#%%
















def GlobGHcp(xB,T):
    import math
    #xB is the conc Zn
    #xA is the conc of Al
    L0 = 18820.95-8.95255*T
    L1 = 1.0E-6+0.00*T
    L2 = 1.0E-6+0.00*T
    L3 = -702.79+0.00*T
    R = 8.314
    # Al hcp
    if (298.15 < T <= 700):
        GAlHcp = -2495.15 + (135.293038*T) - (24.3671976*T*math.log(T)) - ((1.884662E-3)*T**(2)) - ((0.877664E-6)*T**(3)) + (74092.*T**(-1))
        
    elif (700 < T <= 933.473):
        GAlHcp = -5795.24 + (221.248446*T) - (38.5844296*T*math.log(T)) + ((18.531982E-3)*T**(2)) - ((5.764227E-6)*T**(3))+ (74092.*T**(-1))
        
    elif (933.473 < T <= 2900):
        GAlHcp = -5797.361 + (186.884136*T) - (31.748192*T*math.log(T)) - ((1230.622E25)*T**(-9))
        
    # Zn hcp
    if (298.15 < T <= 692.677):
        GZnHcp = -7285.787 + (118.470069*T) - (23.701314*T*math.log(T))- ((1.712034E-3)*T**(2)) - ((1.264963E-6)*T**(3))
        
    elif (692.677 < T <= 1700):
        GZnHcp = -11070.546 + (172.345644*T) - (31.38*T*math.log(T)) + ((470.47E24)*T**(-9)) 
        
    Func = ((1.-xB)*GAlHcp)+(xB*GZnHcp) + (R*T)*((1.-xB)*math.log(1.-xB)+xB*math.log(xB)) + ((1.-xB)*xB)*(L0 + L1*((1.-xB)-xB) +L2*((1.-xB)-xB)**2 + L3*((1.-xB)-xB)**3)
    return Func
#%%