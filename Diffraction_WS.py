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
    y = np.linspace(-d0, d0, 1001)*0.1       #m
    
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
    I[0] = I[-1]
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

#%% Convert wavelength (nm) to approximate RGB values
#http://stackoverflow.com/questions/3407942/rgb-values-of-visible-spectrum
def wl2RGB(wl):
    R = 0.0
    G = 0.0
    B = 0.0
    if (wl >= 400.0) & (wl < 410.0):
        t = (wl - 400.0) / (410.0 - 400.0)
        R = +(0.33*t) - (0.20*t*t)
    elif (wl>=410.0) & (wl<475.0):
        t = (wl - 410.0)/(475.0 - 410.0)
        R = 0.14 - (0.13*t*t)
    elif ((wl>=545.0) & (wl<595.0)):
        t = (wl-545.0)/(595.0-545.0)
        R = +(1.98*t) - (t*t)
    elif (wl>=595.0) & (wl<650.0):
        t = (wl - 595.0)/(650.0 - 595.0)
        R = 0.98 + (0.06*t) - (0.40*t*t)
    elif (wl>=650.0) & (wl<700.0):
        t = (wl - 650.0)/(700.0 - 650.0)
        R = 0.65 - (0.84*t) + (0.20*t*t)
    if (wl>=415.0) & (wl<475.0):
        t = (wl - 415.0)/(475.0 - 415.0)
        G = +(0.80*t*t)
    elif (wl>=475.0) & (wl<590.0):
        t = (wl - 475.0)/(590.0 - 475.0) 
        G = 0.8 + (0.76*t) - (0.80*t*t)
    elif (wl>=585.0) & (wl<639.0):
        t = (wl - 585.0)/(639.0 - 585.0)
        G = 0.84 - (0.84*t)
    if (wl>=400.0) & (wl<475.0):
        t = (wl - 400.0)/(475.0 - 400.0)
        B = +(2.20*t) - (1.50*t*t)
    elif (wl>=475.0) & (wl<560.0):
        t = (wl-475.0)/(560.0-475.0)
        B = 0.7 - t + (0.30*t*t)
    return(R, G, B)
    
#%% Initial Values
#Distance to screen, const
d0 = .30                   #m

# Wavelength/Color
l0 = 0.638e-6               #m
r, g, b = wl2RGB(l0*1e9)

# Slit Width, const
a0 = 100e-6                 #m  

# Number of Point Sources
N0 = 100     

#%% Calculate and plot the diffraction pattern
fig = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.35)

y, I = numDiffraction(d0, a0, l0, N0)

#Plot the intensity vs. screen position in cm
ax1 = plt.subplot(211)
h1, = plt.plot(y*100, I, lw = 3.0, ls = '-', c = (r, g, b))  
#plt.title('Intensity Pattern of ' + str(l*1e6) + ' um Light Diffracted by a ' +
#          str(a*1e6) + ' um slit' + '\n at a distance of ' + str(d*100) + 
#          ' cm, Approximated by ' + str(N) + ' point sources')
plt.title('Single-Slit Diffraction')
plt.ylabel('Insensity (au)')
plt.axis([y[0]*100, y[-1]*100, 0, 1])

# Calculate and plot the Fraunhoffer pattern
F = Fraunhoffer(y, d0, a0, l0)
h2, = plt.plot(y*100, F, lw = 3.0, ls = '--', c = (r, g, b))
#plt.legend(['N-Point Sources', 'Fraunhoffer'])

#Construct RGBA Array for diffraction pattern
#Alpha is mapped to intensity value, while RGB is mapped from wl
R = np.ones((1, len(y)))*r
G = np.ones((1, len(y)))*g
B = np.ones((1, len(y)))*b
A = np.ones((1, len(y)))
for ii in range(len(y)):
    A[:, ii] = I[ii]**0.5
ax2 = plt.subplot(212)
RGBA = np.stack((R,G,B,A), axis = -1)
h3 = plt.imshow(RGBA, extent=[y[0]*100, y[-1]*100, 0, 1], aspect = 'auto')
plt.xlabel('Screen Position (cm)')
#ax2.set_axis_bgcolor((0,0,0))
#axcbar = plt.axes([0.92,0.35,0.01,0.25])
#plt.colorbar(h3, cax = axcbar)

#Add Slider Widgets
axcolor = 'lightgoldenrodyellow'
axslitw = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
axwavl = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
axsrcs = plt.axes([0.25, 0.20, 0.65, 0.03], axisbg=axcolor)

sslitw = Slider(axslitw, 'Slit Width (um)', 1, 200, valinit=a0*1e6)
swavl = Slider(axwavl, 'Wavelength (nm)', 350, 750, valinit=l0*1e9)
ssrcs = Slider(axsrcs, 'Point Sources', 2, 100, valinit=N0)

#Update plots when sliders change
def update(val):
    #Get new values
    a = sslitw.val*1e-6
    l = swavl.val*1e-9
    r, g, b = wl2RGB(l*1e9)
    N = int(ssrcs.val)
    
    #Calcuate/Update new numerical diffraction
    y, I = numDiffraction(d0, a, l, N)
    h1.set_ydata(I)
    h1.set_color((r, g, b))
    
    #Calcuate/Update new analytical diffraction
    F = Fraunhoffer(y, d0, a, l)
    h2.set_ydata(F)
    h2.set_color((r, g, b))

    #Construct RGBA Array for diffraction pattern
    #Alpha is mapped to intensity value, while RGB is mapped from wl
    R = np.ones((1, len(y)))*r
    G = np.ones((1, len(y)))*g
    B = np.ones((1, len(y)))*b
    A = np.ones((1, len(y)))
    for ii in range(len(y)):
        A[:, ii] = I[ii]**0.5
    RGBA = np.stack((R,G,B,A), axis = -1)
    h3.set_data(RGBA)
    plt.draw()
    #fig.canvas.draw_idle()

#Apply updates when sliders move
sslitw.on_changed(update)
swavl.on_changed(update)
ssrcs.on_changed(update)

#Add linear/log Radio buttons
rax = plt.axes([0.025, 0.75, 0.10, 0.10], axisbg=axcolor)
radio = RadioButtons(rax, ('linear', 'log'), active=0)
def ax_type(label):
    ax1.set_yscale(label, nonposy='clip')
    ax1.autoscale(enable = True, axis = 'y', tight = 'none')
    plt.draw()
    #fig.canvas.draw_idle()

radio.on_clicked(ax_type)

#Add a reset button
resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
def reset(event):
    sslitw.reset()
    swavl.reset()
    ssrcs.reset()
    plt.draw()

button.on_clicked(reset)

#Show updated plots
plt.show()
