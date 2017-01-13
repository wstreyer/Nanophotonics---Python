# -*- coding: utf-8 -*-
"""
Created on Sat Dec 17 23:07:09 2016
Will Streyer, UIN: 535574423
Homework 1, problem: Refraction
Function to calculate refraction at various interfaces with constant
indices of refraction.
Output of function should be arrays L, n, R
@author: wstreyer
"""

#%%Import libraries
import matplotlib.pyplot as plt
from numpy import sin, arcsin, radians, degrees, ones, linspace, real

#%%Refraction function
def refraction_WS(thetai, matA, matB, indA, indB):
    #Calculate angle of refraction (Snell's Law)
    indA = complex(indA)
    indB = complex(indB)
    thetar = arcsin((indA / indB) * sin(radians(thetai))) #degrees

    #Find Critical angle
    angc = degrees(real(arcsin(indB / indA))) #degrees
    
    #Plot Angle of Refraction
    plt.plot(thetai, degrees(real(thetar)), 'b-', lw = 2.0)
    #This uses the string concatenation operator +
    mytitle = ('Angle of Refraction' + '\n' + '\n' +
               matA + ' : n = ' + str(indA) + '   ' + matB + ' : n = ' + str(indB))
    plt.title(mytitle)
    plt.xlabel('Angle of Incidence (Degrees)')
    plt.ylabel('Angle of Refraction (Degrees)')
    plt.axis([0, 90, 0, 91])
    
    #Plot Critical Angle
    #ones(n) creates an array of ones of length n.
    #creates a text box on your plot relative to the origin of your current axes.
    plt.plot(ones(len(thetai)) * angc, degrees(real(thetar)), 'r-', lw = 2.0)
    mytext = 'Critical Angle = {:02.2f} deg'.format(angc) 
    plt.text(40, 20, mytext)
    
    plt.show()

#%% Set up the inputs (incident angle array)
thetai = linspace(0, 90, 901)

#%% Part a) %Material names for labels
matA = 'GaAs'
matB = 'Air'

#%% Material indices
#http://refractiveindex.info/?shelf=main&book=GaAs&page=Kachare
indA = 3.3
indB = 1

#%%Calculate and Plot angle of refraction and critical angle
plt.figure()
refraction_WS(thetai, matA, matB, indA, indB)

#%% Part b)
#Material names for labels
matA = 'AlAsSb'
matB = 'InP'

#%%Material indices
#http://scitation.aip.org/content/aip/journal/apl/66/4/10.1063/1.114050
#"3.22<n<3.38 for the AlAsSb layer, at a wavelength of 1600 nm"
indA = 3.30

#http://refractiveindex.info/?shelf=main&book=InP&page=Pettit
indB = 3.16

#%%Calculate and Plot angle of refraction and critical angle
plt.figure()
refraction_WS(thetai, matA, matB, indA, indB);

#% Part c)
#Material names for labels
matA = 'Fiber Core'
matB = 'Fiber Clad'

#Material Indices
indA = 1.46
indB = 1.44

#Calculate and Plot angle of refraction and critical angle
plt.figure()
refraction_WS(thetai, matA, matB, indA, indB);
