# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 10:52:42 2022

@author: Office
"""
#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import seaborn as sns
sns.reset_defaults()
# %matplotlib qt
#%%
##############################################
##############################################



#### dimensions and no. of rebars
b = 1000 #mm
D = 150 #mm



d1 = 135 #mm
d2 = 127 #mm



#### concrete and steel properties
fy = 500 #MPa
Es = 200000 #MPa
fck = 25 #MPa


# ey = 0.87*fy/Es+0.002 #yield strain

# #### stress strain curve of HYSD bars
# ss_curve = np.array([[0, 0], [0.00174, 347.8], [0.00195, 369.6], [0.00226, 391.3], [0.00277, 413], [0.00312, 423.9], [0.00417, 434.8]]) #stress strain curve





##############################################
##############################################
#%%


def area_moment(M11,M22,M12):
    M1 = M11+np.abs(M12)
    M2 = M22+np.abs(M12)
    
########procedure for top reinforcement    
    idx = M1>=0
    M1b = idx*M1
    M2b = M2*idx+(M22+np.abs(M12**2/M11))*(1-idx)
    idx = M2b>=0
    M2b = idx*M2b
    M1b = M1b*idx+(M11+np.abs(M12**2/M22))*(1-idx)
    idx = M1b>=0
    M1b = idx*M1b


########procedure for bottom reinforcement  
    M1 = M11-np.abs(M12)
    M2 = M22-np.abs(M12)  
    idx = M1<=0
    M1t = idx*M1
    M2t = M2*idx+(M22-np.abs(M12**2/M11))*(1-idx)
    idx = M2t<=0
    M2t = idx*M2t
    M1t = M1t*idx+(M11-np.abs(M12**2/M22))*(1-idx)
    idx = M1t<=0
    M1t = idx*M1t
    
  
    return [M1t,M2t,M1b,M2b]

def area_design(M,b,d,fck,fy):
    R = M/(b*d**2)*10**6
    xu = 1.202 *(1-(1-4.597*R/fck)**0.5) *d
    Ast= M*10**6/(0.87*fy*d*(1-0.416*xu/d))
    return Ast


#%% save with units kN, m

df = pd.ExcelFile(r'area_forces.xlsx', engine = "openpyxl") # row 2 contains the units
df_elmforces = pd.read_excel(df,'Element Forces - Area Shells', skiprows = [0,2])
df_elmstresses = pd.read_excel(df,'Element Stresses - Area Shells', skiprows = [0,2])

#%%
M11 = df_elmforces["M11"]
M22 = df_elmforces["M22"]
M12 = df_elmforces["M12"]

[M1t,M2t,M1b,M2b] = area_moment(df_elmforces["M11"],df_elmforces["M22"],df_elmforces["M12"])
WA_moment = pd.DataFrame({"M1t": np.abs(M1t), "M2t": np.abs(M2t),"M1b":M1b,"M2b":M2b})

WA_moment["Ast_1t"] = area_design(WA_moment["M1t"],b,d1,fck,fy) 
WA_moment["Ast_2t"] = area_design(WA_moment["M2t"],b,d2,fck,fy) 
WA_moment["Ast_1b"] = area_design(WA_moment["M1b"],b,d1,fck,fy) 
WA_moment["Ast_2b"] = area_design(WA_moment["M2b"],b,d2,fck,fy) 

df_elmstresses["Smax"] = (df_elmstresses["S23 Average"]**2+df_elmstresses["S13 Average"]**2)**0.5
#%%
df_elmforces =pd.concat([df_elmforces,WA_moment],axis =1, join = "inner")

#%%

ans = df_elmforces[df_elmforces['Joint']=='244']
