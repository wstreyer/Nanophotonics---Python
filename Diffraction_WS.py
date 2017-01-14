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
from matplotlib.widgets import Slider, Button, RadioButtons

#%% Numerically calculate diffraction
# d = Distance to screen (m)
# a = Slit Width (m)
# l = Wavelength of light (m)
# N = Number of Point Sources
#y = %Partition of the screen width (m)
    # i.e y = np.linspace(-d, d, 100001)*2
    
def numDiffraction(d, a, l, N):
    #Partition screen width
    y = np.linspace(-d0, d0, 10001)*2       #m
    
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
    dx = np.zeros(np.shape(sources))
    phase = np.zeros(np.shape(sources))
    
    #Determine the E-field at each point on the screen
    for ii in range(1, len(y), 1):
        
        #Find the path length difference
        #between the current point source and the first one
        dx = (sources[0]-sources)*np.sin(theta[ii])
        
        #Find the phase difference
        #d_phase = k*d_pathlength
        phase = 2*np.pi*dx/l
        
        #To add two waves with a phase difference:
        #Let E1 = E0cos(wt - kx) and E2 = E0cos(wt - kx +d_p),
        #then E1 + E2 = 2E0cos(d_p/2)cos(wt-kx)
        #Since we are going to normalize the intensity, we will just keep
        #track of the cos(d_p/2) contribution to field amplitude
        E[ii] = sum(np.cos(phase/2))
        
    #Find the intensity
    I = E**2
    
    #Find the max intensity
    I_max = max(I)
    
    #Normalize the Intensity
    I /= I_max
    return(y, I)

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
    
#%% Initial Values
#Distance to screen, const
d0 = .3                   #m

# Slit Width, const
a0 = 4e-6                 #m  

# Wavelength
l0 = 1e-6               #m

# Number of Point Sources
N0 = 100     

#%% Calculate and plot the diffraction pattern
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)

y, I = numDiffraction(d0, a0, l0, N0)

#Plot the intensity vs. screen position in cm
h1, = plt.plot(y*100, I, 'b-', lw = 2.0)  
#plt.title('Intensity Pattern of ' + str(l*1e6) + ' um Light Diffracted by a ' +
#          str(a*1e6) + ' um slit' + '\n at a distance of ' + str(d*100) + 
#          ' cm, Approximated by ' + str(N) + ' point sources')
plt.title('Single-Slit Diffraction')
plt.xlabel('Screen Position (cm)')
plt.ylabel('Normalized Insensity (au)')
plt.axis([y[0]*100, y[-1]*100, 0, 1])

# Calculate and plot the Fraunhoffer pattern
F = Fraunhoffer(y, d0, a0, l0)
h2, = plt.plot(y*100, F, 'g-', lw = 2.0)

# Calculate and plot the minima
y_n = diffraction_minima(d0, l0, a0, 3)
h3, = plt.plot(y_n*100, np.zeros(len(y_n)), 'ro', mfc = 'none', mec ='r', mew = 10.0, ms = 2.0)
plt.legend(['N-Point Sources', 'Fraunhoffer', 'Expected Minina'])

axcolor = 'lightgoldenrodyellow'
axslitw = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
axwavl = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)

sslitw = Slider(axslitw, 'Slit Width (um)', 0.1, 30.0, valinit=a0*1e6)
swavl = Slider(axwavl, 'Wavelength (um)', 0.1, 10.0, valinit=l0*1e6)

def update(val):
    a = sslitw.val
    l = swavl.val
    y, I = numDiffraction(d0, a, l, N0)
    F = Fraunhoffer(y, d0, a, l)
    y_n = diffraction_minima(d0, l, a, 6)
    h1.set_ydata(I)
    h2.set_ydata(F)
    h3.set_data(y_n*100, np.zeros(len(y_n)))
#    fig.sca(h)
#    plt.plot(y*100, I, 'b-', lw = 2.0)
#    plt.plot(y*100, F, 'g-', lw = 2.0)
#    plt.plot(y_n*100, np.zeros(len(y_n)), 'ro', mfc = 'none', mec ='r', mew = 10.0, ms = 2.0)
    fig.canvas.draw_idle()
sslitw.on_changed(update)
swavl.on_changed(update)
plt.show()
