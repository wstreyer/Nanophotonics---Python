"""
Created on Thu Jan 12 21:08:50 2017
#Will Streyer, UIN: 535574423
#Homework 1, problem: Diffraction
#Calculates the intensity of light diffracted by a single slit onto a
#screen. The wavelenth is fixed at 1 micron and the distance to the screen
#is 30 cm. Displays the diffraction pattern. You can change the slit width
#and the number of points sources used to approximate the pattern.
@author: wstreyer
"""
#%% Import modules
import matplotlib.pyplot as plt
import numpy  as np

#%% Numerically calculate diffraction
# d = Distance to screen (m)
# a = Slit Width (m)
# l = Wavelength of light (m)
# N = Number of Point Sources
#y = %Partition of the screen width (m)
    # i.e y = np.linspace(-d, d, 100001)*2
    
def numDiffraction(y, d, a, l, N):  
    # Find angle between slit and 
    # each point on the screen
    theta = np.arctan(y/d)   #radians
    
    #Partition slit into N
    #point sources
    sources = np.linspace(-a/2, a/2, N)    #m 
    
    #Pre-allocate E-field array
    #This will hold the value of the
    #total field amplitude at each
    #point on the screen
    E = np.zeros(np.shape(y))
    
    #Determine the E-field at each point on the screen
    for ii in range(1, len(y), 1):
        
        #Sum all the phase differences for each point source
        for kk in range(1, len(sources), 1):
            
            #Find the path length difference
            #between the current point source and the first one
            dx = (sources[0]-sources[kk])*np.sin(theta[ii])
            
            #Find the phase difference
            #d_phase = k*d_pathlength
            phase = 2*np.pi*dx/l
            
            #To add two waves with a phase difference:
            #Let E1 = E0cos(wt - kx) and E2 = E0cos(wt - kx +d_p),
            #then E1 + E2 = 2E0cos(d_p/2)cos(wt-kx)
            #Since we are going to normalize the intensity, we will just keep
            #track of the cos(d_p/2) contribution to field amplitude
            E[ii] += np.cos(phase/2)

    #Find the intensity
    I = E**2
    
    #Find the max intensity
    I_max = max(I)
    
    #Normalize the Intensity
    I /= I_max
    return(I)

#%% Compute Fraunhoffer diffraction pattern
# The Fraunhoffer formula is the limit as
# the number of point sources goes to infinity
def Fraunhoffer(y, d, a, l):
    theta = np.arctan(y/d)
    beta = np.pi/l*a*np.sin(theta)
    I = np.sinc(beta/np.pi)**2
    return(I)

#%% Find the expected minima
# o is the number of +/- orders to calculate
def diffraction_minima(d, l, a, o):
    # Enumerate orders of minima.
    order = np.concatenate([np.arange(-o, 0, 1), np.arange(1, o+1, 1)])
    
    # Calculate the angle at which minima are expected
    theta = np.arcsin(order*l/a)
    
    # Convert that angle to a screen position
    y_n = d*np.tan(theta)
    
    #Remove imaginary minima
    y_n[np.where(np.abs(np.real(y_n))<1e-15)] = np.nan
    return(y_n)

#%% Constants
#Distance to screen, const
d = .3                   #m

# Slit Width, const
a = 4e-6                 #m  

# Wavelength
l = 1e-6               #m

# Number of Point Sources
N = 25     

#Partition screen width
y = np.linspace(-d, d, 1001)*2         #m 

#%% Calculate and plot the diffraction pattern
plt.figure()
I = numDiffraction(y, d, a, l, N)

#Plot the intensity vs. screen position in cm
plt.plot(y*100, I, 'b-', lw = 2.0)  
plt.title('Intensity Pattern of ' + str(l*1e6) + ' um Light Diffracted by a ' +
          str(a*1e6) + ' um slit' + '\n at a distance of ' + str(d*100) + 
          ' cm, Approximated by ' + str(N) + ' point sources')
plt.xlabel('Screen Position (cm)')
plt.ylabel('Normalized Insensity (au)')
plt.axis([y[0]*100, y[-1]*100, 0, 1])
plt.show()

# Calculate and plot the Fraunhoffer pattern
F = Fraunhoffer(y, d, a, l)
plt.plot(y*100, F, 'g-', lw = 2.0)

# Calculate and plot the minima
y_n = diffraction_minima(d, l, a, 3)
plt.plot(y_n*100, np.zeros(len(y_n)), 'ro', mfc = 'none', mec ='r', mew = 10.0, ms = 2.0)
plt.legend(['N-Point Sources', 'Fraunhoffer', 'Expected Minina'])
