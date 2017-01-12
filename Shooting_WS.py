# -*- coding: utf-8 -*-
"""
Created on Sun Dec 18 22:34:50 2016
HW 6, extra credit
Calculates the energy and wavefunctions with arbitrary potential profile using the shooting method.
@author: wstreyer
"""
#%%Import Libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from timeit import default_timer as timer

#%% Physical Constants
q = 1.6E-19 # C
m0 = 9.11E-31 # kg
hbar = 1.05E-34 # Js

#%% Structure Definition
#Layer Thicknesses (nm)
zunit = 1.0e-9
thickness = np.array([20, 8, 0.5, 4, 20]) * zunit

#Location of Interfaces
zj = np.zeros((len(thickness) + 1))
for kk in range(0, len(thickness)):
    zj[kk+1] = zj[kk] + thickness[kk]

#%%Layer Potentials (J)
V0 = 0.5 #barrier height, eV
potential = np.array([1, 0, 1.5, 0.1, 1]) * V0 * q  
    
#%%Layer Masses (kg)
#AlxGa1-xAs
x = 0.3
me = lambda x: 0.067 + 0.083 * x
mass = np.array([me(0.3), me(0.0), me(0.3), me(0.0), me(0.3)]) * m0  
 
#%% Build a mesh for z, V, and m
dz = 0.05e-9; #nm
zj = zj / dz   
    
z = np.arange(dz, sum(thickness), dz)
V = np.ones(np.shape(z))
m = np.ones(np.shape(z))    
    
#Assign correct potential and mass to each region
for ii in range(0, len(zj) - 1):
    region = np.arange(int(zj[ii]) - 1, int(zj[ii + 1]) - 1, 1)
    V[region] = potential[ii] 
    m[region] = mass[ii]
   
#%% Add E-field
#Linear example (i.e. uniform electri  field)
#EF = -0.002 * z / zunit * q;     #V/nm = kV/um = MV/mm = 10 MV/cm

#Quadratic example (approx. QHO)
EF = 0.001/8*(z/zunit-40)**2*q     #V/nm = kV/um = MV/mm = 10 MV/cm
V = V + EF   

#%% Plot the potetial
plt.figure()    
plt.plot(z / zunit, V / q * 1000, 'b-', lw = 2.0)
plt.xlabel('Position (nm)')
plt.ylabel('Potential (meV)')
plt.show()

#%% Determine approximate eigen-values
#Form a 1D mesh in Energy Space
N1 = 500
E = np.linspace(min(V), max(V), N1)

#Transform potential into k**2 space
W = 2*m*V/hbar**2

#Initial Solutions for Psi
#Psi @ z = 0 must be zero
#Psi[1] can be any value, it only effects the magnitude of Psi, which will be normalized anyway
Y = np.ones(len(z))
Y[0] = 0
Y[1] = 1

#The key to this method is using the sign of the open end of Psi to refine the precision of eigen-values
Y_sign = np.zeros(np.shape(E));

for ii in np.arange(0, len(E)):
    #Transform E into k^2 space
    e = 2*m*E[ii]/hbar**2
    
    #Iteratively solve finite difference form of Schroedinger equation (Numerov Approximation)
    for jj in  np.arange(2, len(z)):
        Y[jj] = dz**2*(W[jj-1] - e[jj-1])*Y[jj-1] + 2*Y[jj-1] - Y[jj-2]
    
    #Find sign of the open end of Psi
    Y_sign[ii] = np.sign(Y[-1])

#Determine energies where open end of Psi changes sign
ndx = np.where(abs(np.diff(Y_sign)) > 0)
En = E[ndx]
print(En/q*1000)

#%% Create a new E-array that zooms in on each En value
N2 = 10     #this value can be quite small for better speed and still see good convergence
dE = np.diff(E)
dE = dE[0]
EE = np.zeros((len(En), N2))
for kk in np.arange(0, len(En)):
    EE[kk, :] = np.linspace(En[kk]-dE, En[kk]+dE, N2)

#Set a color index
cc = [cm.hsv(x) for x in np.linspace(0, 1, len(En))]
    
#Initialize final eigen-functions
Psi = np.zeros((len(z), len(En)))

#Determine correct Psi for each En
#The recursive equation here could see improved performance with Cython or PyPy
for kk in  np.arange(0, len(En)):
    Y_sign = np.zeros(np.shape(EE[kk, :]))
    
    #Refine the value of En until Psi converges
    start = timer()
    converged = 0
    while not converged:
        for ii in np.arange(0, len(EE[kk, :])):
            #Transform EE into k^2 space
            e = 2*m*EE[kk,ii]/hbar**2
            
            #Iteratively solve finite difference form of Schroedinger equation (Numerov Approximation)
            for jj in np.arange(2, len(z)):
                Y[jj] = dz**2*(W[jj-1] - e[jj-1])*Y[jj-1] + 2*Y[jj-1] - Y[jj-2]
            
            #Find sign of the open end of Psi
            Y_sign[ii] = np.sign(Y[-1])
        
        #Determine energies where open end of Psi changes sign               
        ndx = np.where(abs(np.diff(Y_sign)) > 0)
        En[kk] = EE[kk, ndx[0][-1]]
        
        #There should be only one eigen-value per EE range, note if there are more
        if len(ndx[0]) > 1:
           print('multiple eigen-values: {}'.format(len(ndx[0])))
           text = 'eigen values: '
           for i in range(0,len(ndx[0]),1):
               text +=  '{:0.2f}, '
           text = text[:-2]   
           print(text.format(*EE[kk,ndx[0]]/q*1000) + ' meV' + ' dE = {:0.2E}.'.format(dE/q*1000)) 
        
        #Check for convergence, max(Psi) > Psi[-1], i.e. max(Psi) ~= Psi[-1]
        #R is a measure of how well the eigen-function converges
        #The best solution would have Psi[-1] -> 0, just like Psi[0]
        R = np.max(abs(Y))/abs(Y[-1])
        if R > 1:
           converged = 1
           end = timer()
           print('{:0.2f} meV - converged - R = {:0.2f} - {:0.2f} s'.format(En[kk]/q*1000, R, (end - start)))
        else:
            #If not converged, refine eigen-value search
            dE = np.diff(EE[kk, :])
            dE = dE[0]
            EE[kk, :] = np.linspace(En[kk] - dE, En[kk] + dE, N2)
            
            #If the resolution of the range is less than the data-type resolution,
            #consider Psi to be converged, and note the exception
            if np.where(np.diff(EE[kk, :]) == 0)[0].any():
               converged = 1
               print('{} meV - max resolution - R = {}'.format(En[kk]/q*1000, R))
    
    #Check if the open end of |Psi|**2 still diverges
    #This is done by looking at the sign of the derivative
    sign_diff = np.sign(np.diff(Y**2))
    
    #If Psi diverges, set Psi to zero after the slope changes from - to +
    if sign_diff[-1] > 0:
        ndx2 = np.where(np.diff(sign_diff) > 0)
        Y[ndx2[0][-1]:] = 0
    
    #Normalize Psi
    C = np.sqrt(np.sum(Y*Y))
    Psi[:, kk] = Y/C
    
    #Plot En and Psi_n over the band diagram using an appropriate scaling factor
    plt.plot(z/zunit, np.ones(np.shape(z))*En[kk]/q*1000, lw = 3.0, ls = '-', c = cc[kk])
    plt.plot(z/zunit, 5*1000*Psi[:, kk]**2 + En[kk]/q*1000, lw = 3.0, ls = '--', c = cc[kk])
   