# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 18:28:52 2022

@author: prani
"""
# from IPython import get_ipython
# get_ipython().magic('reset -sf')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import multiprocessing
from multiprocessing import Pool
from functools import partial
import time


## add a load combination

#%%
# global As
#%%
def beam_design(Mu, Vu, Tu, D, nc, d, b, fck,  fy, Es):
    [Ast, Asc, Asv_tor2, Asv_tor4, Asv_shear, Asvsv2,Asvsv4,Me1,Me2,Veq] = [0,0,0,0,0,0,0,0,0,0]
    dc=nc
    xumax = d* (0.0035/(0.0055+0.87*fy/Es)) #mm
    Mulim = 0.362*fck*b*xumax*(d-0.416*xumax) #Nmm
    Vu = np.abs(Vu)
    Tu = np.abs(Tu)
    if Mu<0:
        loc = -1 #-1 for top

        Mu = np.abs(Mu)
    else:
        loc = 1
    
    
    # equivalent bending moment
    Mt = Tu*((1+D/b)/1.7)
    
    Me1 = Mt+Mu ## for tension steel
    Me2 = Mt-Mu ## for compression steel
    
    if Me1 < Mulim:
        R = Me1/(b*d**2)
        xu = 1.202 *(1-(1-4.597*R/fck)**0.5) *d
        Ast= Me1/(0.87*fy*d*(1-0.416*xu/d))
        Ast = max(Ast,(0.85/fy)*b*d)
    
    else: 
# doubly reinforced    
        DAst = (Me1- Mulim)/(0.87*fy*(d-dc)) ##mm
        Astlim = Mulim/(0.87*fy*d*(1-0.416*xumax/d))
        Ast = DAst + Astlim
        esc = (xumax-dc)*0.0035/xumax
        if esc >= ss_curve[6,0]:
            fsc = fy/1.15;
        else:
            ans = ss_curve[:,0]-esc
            lb = ss_curve[ans == ans[ans>0].min(),:]
            ub = ss_curve[ans == ans[ans<0].max(),:]
            fsc = lb[0,1]+(esc-lb[0,0])*(ub[0,1]-lb[0,1])/(ub[0,0]-lb[0,0])
    
        Asc = 0.87*fy*DAst/(fsc - 0.447 *fck)
        rho_t = Ast/b/d
        rho_c = Asc/b/d
    

        
    if Me2 > 0:
        R2 = Me2/(b*d**2)
        xu = 1.202 *(1-(1-4.597*R2/fck)**0.5) *d
        Ast2= Me2/(0.87*fy*d*(1-0.416*xu/d))
        Ast2 = max(Ast2,(0.85/fy)*b*d)
        Asc= Asc+Ast2

# design for shear
    tau_v = Vu/b/d  ## MPa
    beta = max((0.8*fck)/(6.89*(Ast/b/d)),1)
    tau_c = 0.85*(0.8*fck)**0.5*((1+5*beta)**0.5-1)/(6*beta)
    tau_cmax = 0.62*fck**0.5
    
    VyRlim = tau_cmax*b*d/1000
    
    Ve = Vu+1.6*Tu/b
    tau_ve = Ve/b/d
    
    # 2 legged hoop 
    b1 = b-2*50
    d1 = d-2*50 
    b2 = b1-2*50
    
    Veq = tau_ve*b*d
    if tau_ve>tau_c:
        Asv_tor2 = (Tu/(b1*d1*(fy/1.15)) + Vu/(2.5*d1*(fy/1.15)))*1000
        Asv_tor4 = (Tu/((fy/1.15)*(b1*d1+b2*d1)) + Vu/(2.5*d1*(fy/1.15)))*1000
        Asv_shear = (tau_ve-tau_c)*b/(fy/1.15)*1000
        Asvsv2 = max(Asv_tor2,Asv_shear)  #mm2/m
        Asvsv4 = max(Asv_tor4,Asv_shear)  #mm2/m
    else:
        Asvsv2 = 0.4/(fy/1.15)*1000 #mm2/m
        Asvsv4 = 0.4/(fy/1.15)*1000 #mm2/m

    if Asvsv2 < 0.4/(fy/1.15):
        Asvsv2 = 0.4/(fy/1.15) *1000 #mm2/m
        Asvsv4 = 0.4/(fy/1.15) *1000 #mm2/m
    
    if Ast > 0.04*b*D:
        Ast = "max lim"
    
    if Asc > 0.04*b*D:
        Asc = "max lim"
   
       
    return [loc*Ast, Asc, Asv_tor2, Asv_tor4, Asv_shear, Asvsv2, Asvsv4, Me1/10**6, Me2/10**6, Veq/1000]

def sleepy_man(df_elmforces, df_sections, df_frameassigns,i):
    [Ast, Asc, Asv_tor2, Asv_tor4, Asv_shear, Asvsv2,Asvsv4,Me1,Me2,Veq] = [0,0,0,0,0,0,0,0,0,0]
    # As=pd.DataFrame(columns = ['beam id','elm station','Ast','Asc','Asv_tor2','Asv_tor4','Asv_shear','Asvsv2','Asvsv4','Me1', 'Me2', 'Veq'])

    
    ll1 = df_elmforces.iloc[i]['Unique Name'] ## unique number of the beam
    ll2 = df_frameassigns[df_frameassigns['UniqueName'] == ll1]['Section Property'].iloc[0]
    # ll2 = df_frameassigns[df_frameassigns['UniqueName'] == ll1]['Section Property'] ##frame section assignment
    ll3 = df_sections[df_sections['Name'] == ll2] ##section data
    D = ll3["Depth"].iloc[0]
    b = ll3["Width"].iloc[0]
    d = D- nc ## effective depth of steel
    dc = nc ## mm depth of compression steel
    P = df_elmforces.iloc[i]['P']*1000 #N
    V2 = df_elmforces.iloc[i]['V2']*1000 #N
    V3 = df_elmforces.iloc[i]['V3']*1000 #N
    T = df_elmforces.iloc[i]['T']*10**6 #Nmm
    M2 = df_elmforces.iloc[i]['M2']*10**6 #Nmm
    M3 = df_elmforces.iloc[i]['M3']*10**6 #Nmm
    
    [Ast, Asc, Asv_tor2, Asv_tor4, Asv_shear, Asvsv2,Asvsv4,Me1,Me2,Veq] = beam_design(M3, V2, T, D, nc, d, b, fck,  fy, Es)
    # As = As.append({'beam id':ll1,'elm station':df_elmforces.iloc[i]['Station'] , 'Ast':Ast,'Asc':Asc,'Asv_tor2':Asv_tor2,'Asv_tor4':Asv_tor4,'Asv_shear':Asv_shear,'Asvsv2':Asvsv2,'Asvsv4':Asvsv4,'Me1':Me1, 'Me2': Me2, 'Veq': Veq}, index = [i])
    
    
    As = pd.DataFrame([[ll1,df_elmforces.iloc[i]['Station'],Ast, Asc, Asv_tor2, Asv_tor4, Asv_shear, Asvsv2,Asvsv4,Me1,Me2,Veq]],columns = ['beam id','elm station','Ast','Asc','Asv_tor2','Asv_tor4','Asv_shear','Asvsv2','Asvsv4','Me1', 'Me2', 'Veq'], index=[i])
    return As
    # return [Ast, Asc, Asv_tor2, Asv_tor4, Asv_shear, Asvsv2,Asvsv4,Me1,Me2,Veq]
#%%
ss_curve = np.array([[0, 0], [0.00174, 347.8], [0.00195, 369.6], [0.00226, 391.3], [0.00277, 413], [0.00312, 423.9], [0.00417, 434.8]])
Es = 2*10**5 ## MPa young's modulus of steel
fck = 25 ## MPa compressive strength of concrete
fy = 500 ## MPa tensile strength of steel
nc = 60 ## mm nominal cover to beam

#%%


l = 7750 ## mm length of beam
# D = 550 ## mm overal depth of beam
# b = 450 ## mm width of beam


#%% save with units kN, m
df = pd.ExcelFile(r'beam_data.xlsx', engine = "openpyxl") # row 2 contains the units
#%% 
df_elmforces = pd.read_excel(df,'Element Forces - Beams', skiprows = [0,2])
df_sections = pd.read_excel(df,'Frame Sec Def - Conc Rect', skiprows = [0,2])
df_frameassigns = pd.read_excel(df,'Frame Assigns - Sect Prop', skiprows= [0,2])
# df_rebar = pd.read_excel(df,'Mat Prop - Rebar Data', skiprows= [0,2])
# df_concrete = pd.read_excel(df,'Mat Prop - Concrete Data', skiprows= [0,2])

beam_num  = df_elmforces['Unique Name']
beam_num = beam_num.drop_duplicates()
#%%
nsample = np.shape(df_elmforces)[0]
As=pd.DataFrame(columns = ['beam id','elm station','Ast','Asc','Asv_tor2','Asv_tor4','Asv_shear','Asvsv2','Asvsv4','Me1', 'Me2', 'Veq'])




tic= time.time()
process_list =[]

if __name__ == '__main__':
    pro = Pool()
    func = partial(sleepy_man,df_elmforces, df_sections, df_frameassigns )
    As = As.append(pro.map(func,[i for i in range(nsample)]))
    pro.close()
    pro.join()

    
toc = time.time()   


print('Done in {:.4f} seconds'.format(toc-tic))    
    
  
# df_elmforces = df.parse('Element Forces - Beams')
# df_sections = df.parse('Frame Sec Def - Conc Rect')
#%% inclusion of NBC for shear, moment and torsion



# df_frameassigns['Label'] == df_elmforces.iloc[0]['Unique Name']
# df_frameassigns[df_frameassigns['UniqueName'] == ll1]
# df_elmforces[df_elmforces['Unique Name'] == ll1]

#%% slenderness ratios
# lD_ratio = l/d ## should be below 26 for continuous, below 7 for cantilever
# cdist = min(60*b,250*b**2/d)  ## should be higher than l

df_elmforces =pd.concat([df_elmforces,As],axis =1, join = "inner")


#%%
rho_st=pd.DataFrame(columns = ['beam id','elm station','top','bot','Asvsv'])

for j in beam_num:
    ans = df_elmforces[df_elmforces['Unique Name']==j]
    elem_st = ans['Station']
    elem_st = elem_st.drop_duplicates()


    for i in elem_st:
        rho_st = rho_st.append({'beam id':j,'elm station':i , 'top':max(np.abs(min(ans[ans['Station']==i]['Ast'])),np.abs(max(ans[ans['Station']==i]['Asc'])))
                        ,'bot':max(np.abs(max(ans[ans['Station']==i]['Ast'])),np.abs(max(ans[ans['Station']==i]['Asc']))),'Asvsv':max(ans[ans['Station']==i]['Asvsv2'])}, ignore_index=True)
       

#%%
# %matplotlib qt
import seaborn as sns
sns.set_theme()
no=223
ans = df_elmforces[df_elmforces['Unique Name']==no]
plt.rcParams.update({'font.size': 10,'font.family': 'serif'})

fig, axes = plt.subplots(2, 3, figsize=(12.8,4.8))
                         
# plt.figure(figsize =[6.4,4.8],tight_layout = True)
sns.lineplot(ax= axes[0,0], x= ans['Station'],y=ans['Ast'], style=ans['Output Case'],ci =None, legend = False)

sns.lineplot(ax= axes[0,1],x= ans['Station'],y=ans['Asc'],  style=ans['Output Case'],ci =None , legend = False)

# plt.figure(figsize =[6.4,4.8],tight_layout = True)
sns.lineplot(ax= axes[0,2],x= ans['Station'],y=ans['Asvsv2'],  style=ans['Output Case'],ci =None, legend = False)

# plt.legend(bbox_to_anchor=(-0.75, 1), loc='lower center', borderaxespad=0, ncol = 6, frameon=False)
# plt.text(4.5,-0.01,['beam',no])

# fig, axes = plt.subplots(1, 3, figsize=(12.8,4.8))
                         
# plt.figure(figsize =[6.4,4.8],tight_layout = True)
sns.lineplot(ax= axes[1,0], x= ans['Station'],y=ans['M3'], style=ans['Output Case'],ci =None, legend = False)

sns.lineplot(ax= axes[1,1],x= ans['Station'],y=ans['V2'],  style=ans['Output Case'],ci =None , legend = False)

# plt.figure(figsize =[6.4,4.8],tight_layout = True)
sns.lineplot(ax= axes[1,2],x= ans['Station'],y=ans['T'],  style=ans['Output Case'],ci =None)

plt.legend(bbox_to_anchor=(-0.75, 2.2), loc='lower center', borderaxespad=0, ncol = 6, frameon=False)

plt.text(4.5,-0.01,['beam',no])
# #%% limiting moment of resistance
# import seaborn as sns
# sns.set_theme()
# for i in range(0,np.shape(beam_num)[0]):
#     ans = df_elmforces[df_elmforces['Unique Name']==beam_num.iloc[i]]
    
#     plt.rcParams.update({'font.size': 10,'font.family': 'serif'})
    
#     fig, axes = plt.subplots(1, 3, figsize=(12.8,4.8))
                             
#     # plt.figure(figsize =[6.4,4.8],tight_layout = True)
#     sns.lineplot(ax= axes[0], x= ans['Station'],y=ans['Ast'], hue=ans['Output Case'],ci =None, legend = False)
    
#     sns.lineplot(ax= axes[1],x= ans['Station'],y=ans['Asc'], hue=ans['Output Case'],ci =None , legend = False)
    
#     # plt.figure(figsize =[6.4,4.8],tight_layout = True)
#     sns.lineplot(ax= axes[2],x= ans['Station'],y=ans['Asvsv'], hue=ans['Output Case'],ci =None)
    
#     plt.legend(bbox_to_anchor=(-0.75, 1), loc='lower center', borderaxespad=0, ncol = 6, frameon=False)
    
#     plt.text(4.5,-0.01,['beam', beam_num.iloc[i]])

