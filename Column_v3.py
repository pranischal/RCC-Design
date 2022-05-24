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
#%%
##############################################
##############################################

#### dimensions and no. of rebars
b = 600 #mm
D = 700 #mm
As = np.array([1963.5,981.75,981.75,1963.5]) #include for no. of rows of reinforcing steel
nc = 60 #mm, cover to the centroid of the outermost longitudinal bar

#### set moment
moment = "M2"  # either M2 or M3

#### concrete and steel properties
fy = 500 #MPa
Es = 200000 #MPa
fck = 25 #MPa


ey = 0.87*fy/Es+0.002 #yield strain

#### stress strain curve of HYSD bars
ss_curve = np.array([[0, 0], [0.00174, 347.8], [0.00195, 369.6], [0.00226, 391.3], [0.00277, 413], [0.00312, 423.9], [0.00417, 434.8]]) #stress strain curve





##############################################
##############################################
#%%


def stress_calc(es, fy, fck, ss_curve):
    i = np.size(es)
    fs = np.zeros(i) # no. of rebar layers
    fc = np.zeros(i) 
    for i in range(0,i):
        if np.abs(es[i]) >= ss_curve[6,0]:
            fs[i] = fy/1.15;
        
        else:
            ans = ss_curve[:,0]-np.abs(es[i])
            ub = ss_curve[ans == ans[ans>0].min(),:]
            lb = ss_curve[ans == ans[ans<=0].max(),:]
            fs[i] = lb[0,1]+(np.abs(es[i])-lb[0,0])*(ub[0,1]-lb[0,1])/(ub[0,0]-lb[0,0])
               
                    
        if es[i]<=0:
            fc[i] = 0
            fs[i] = -1*fs[i]
        elif es[i]>=0.002:
            fc[i] = 0.446 *fck
        else:
            fc[i] = 0.446*fck*(2*(es[i]/0.002)-(es[i]/0.002)**2)
    return [fs,fc]

    

#%%
x = np.zeros([100])
Pu_norm = np.zeros([100])
Mu_norm = np.zeros([100])


d = D-nc
dc = D-d
rho = As/b/D*100

i = np.size(As) # no. of rebar layers
y = np.size(As) 
sc = (d-dc)/(i-1) #spacing between two inner rebars, (equal spacing assumed)


x_loc = np.zeros(i) 
#y = -x+D/2
for count in range(0,i):
    x_loc[count] = nc + count*sc
    
y = -x_loc+D/2
#%% pure axial compression
esc= 0.002
rho_T  = np.sum(rho)
ans = ss_curve[:,0]-esc
lb = ss_curve[ans == ans[ans>0].min(),:]
ub = ss_curve[ans == ans[ans<0].max(),:]
fsc = lb[0,1]+(esc-lb[0,0])*(ub[0,1]-lb[0,1])/(ub[0,0]-lb[0,0])
P_norm  = 0.446 + rho_T/100/fck * (fsc -0.446*fck) ### Pu/fck/b/D normalized force, fsc is stress corresponding to esc = 0.002

#%% generation of interaction diagoram
i = 0
for k in np.linspace(0.01,4,100):
    # k = i/100*100
    xu = k*D

############### stress block parameters with neutral axis outside the section
# k=1.2
# xu = k *D #for construction of interaction diagram, it is enough to consider  upto k = 1.2

    if k>1:
        
        g = 0.446*fck*(4/(7*k-3))**2 # g is the difference in stress between higly compressed edge and least compressed edge
        Asb = 0.446*fck*D*(1-(4/21)*(4/(7*k-3))**2) # area of the compression stress block
        Msb = 0.446*fck*D**2/2-(8/49)*g*D**2 # moment about highly compressed edge
        xsb = Msb/Asb # position of the centroid of the stress block from the highly compressed edge
        
        C1 = 0.446*(1-(4/21)*(4/(7*k-3))**2)
        C2 = xsb/D
        
        es = 0.0035*(xu-x_loc)/xu
        # subscript 1 denotes the side closer to highly compressed edge
        

        [fs,fc] = stress_calc(es, fy, fck, ss_curve)
                    
        Pu_norm[i] = C1 + np.sum(rho/100/fck*(fs-fc))  # Pu/fck/b/D
        Mu_norm[i] = C1*(0.5-C2) + np.sum(rho/100/fck*(fs-fc)*y/D)  # Mu/fck/b/d/D**2
    
#################  when neutral axis lies within the section
    else:
    # k = xu/D
        # est = -0.0035*(1-k)/k
        # es = (0.0035+np.abs(est))*(D-x_loc)/(D) - np.abs(est)
        es = -0.0035*x_loc/k/D+0.0035

        [fs,fc] =  stress_calc(es, fy, fck, ss_curve)
     
         
        
        Pu_norm[i] = 0.36*k + np.sum(rho/100/fck*(fs-fc)) 
        Mu_norm[i] = 0.36*k*(0.5-0.416*k)+ np.sum(rho/100/fck*(fs-fc)*y/D)
        

    i = i +1

Pu = Pu_norm*b*D*fck/1000
Mu = Mu_norm*b*D**2*fck/1000/1000
#%%
#%% save with units kN, m
df = pd.ExcelFile(r'column_forces.xlsx', engine = "openpyxl") # row 2 contains the units
df_elmforces = pd.read_excel(df,'Element Forces - Columns', skiprows = [0,2])

# %matplotlib qt

sns.reset_defaults()
plt.rcParams.update({'font.size': 10,'font.family': 'serif'})
fig = plt.figure(figsize =[4.8,6.4],tight_layout = True)

plt.plot(Mu,Pu)
plt.plot(np.abs(df_elmforces[moment]),(-1*df_elmforces['P']),"x")
# plt.plot(0,P_norm,'x')

plt.grid(True, color = '#ababab', which = 'major', axis = 'y', linestyle = '--', linewidth = 0.75)
# plt.grid(True, color = '#d6d6d6', which = 'minor', axis = 'y', linestyle = ':', linewidth = 0.75)
plt.minorticks_on()
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['bottom'].set_position(('data',0))
plt.gca().spines['bottom'].set(zorder = 2)
plt.xlabel([moment+' ,kNm'])
plt.ylabel('$P$, kN')
plt.xlim(0,)

plt.savefig("P-M"+moment+".svg")
plt.show()

# df = pd.ExcelFile(r'column_forces.xlsx', engine = "openpyxl") # row 2 contains the units
# df_elmforces = pd.read_excel(df,'Element Forces - Columns', skiprows = [0,2])

# # %matplotlib qt
# import seaborn as sns
# sns.reset_defaults()
# plt.rcParams.update({'font.size': 10,'font.family': 'serif'})
# fig = plt.figure(figsize =[6.4,4.8],tight_layout = True)

# plt.plot(Mu_norm,Pu_norm)
# plt.plot(np.abs(df_elmforces[moment])/fck/b/D**2*10**6,(-1*df_elmforces['P'])/fck/b/D*1000,"x")
# # plt.plot(0,P_norm,'x')

# plt.grid(True, color = '#ababab', which = 'major', axis = 'y', linestyle = '--', linewidth = 0.75)
# # plt.grid(True, color = '#d6d6d6', which = 'minor', axis = 'y', linestyle = ':', linewidth = 0.75)
# plt.minorticks_on()
# plt.gca().spines['top'].set_visible(False)
# plt.gca().spines['right'].set_visible(False)
# plt.gca().spines['bottom'].set_position(('data',0))
# plt.gca().spines['bottom'].set(zorder = 2)
# plt.xlabel('$M/f_{ck}bD^2$')
# plt.ylabel('$P/f_{ck}bD$')
# plt.xlim(0,)
# plt.show()

# plt.hlines(0,0,0.15)
# plt.vlines(0,-0.2,1)
