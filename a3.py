# Unitless 1D FDTD solution to wave equation. 
# Author- Kabir Bansod
# Date: 10/09/14


import numpy as np							# Loading required libraries.
import matplotlib.pyplot as plt
import time
import matplotlib.animation as animation

def get_source():							# Gets the type of source to be used for simulation
    """Gets the source to be used in the calculations from 
       the user"""
    print ("Pick a source to be used in 1D FDTD simulation of wave prpogation in vacuum.")
    print ("1: Sinusoidal source of frequency 1Hz.")
    print ("2: Impulse.")
    print ("3: Gaussian source.")
    return input("Enter your choice: ") 

def get_courantfactor():						# Gets the user input on the courant factor to be used.
    """ Gets the value of Courant factor to be used"""
    return input("Enter the value of courant factor to be used :")

def main():
    #defining dimensions							
    
    xdim=720								# Dimensions(Length; since we are simlating for 1d) of medium.
    time_tot = 1000							# Time for which the simulation runs.
    xsource = xdim/2							# Location of the source.
    
    # Stability factor
    S = get_courantfactor()						
    print ("Note: The value of Courant factor should be between 0 and 1 for stabilty")

    #Speed of light
    c=1.0								# Since we are simulating wave for unitless/normalized values.
    epsilon0=1.0						
    mu0=1.0

    delta =1.0								# Space step.
    deltat = S*delta/c							# Time step.

    Ez = np.zeros(xdim)							# Arrays that hold the values of Electric field and 
    Hy = np.zeros(xdim)							# magnetic field respectively.

    epsilon = epsilon0*np.ones(xdim)					# Arrays that hold the permittivity and permeability values of the 
    mu = mu0*np.ones(xdim)						# medium (medium charecteristics).
    
    #epsilon[3*xdim/4] =100						# Impurity (point impurity)

    source_choice = get_source()					# Get the type of source to be used from the user.

    fig , axes = plt.subplots(2,1)					# Figure with 2 subplots- one for E field and one for H.
    axes[0].set_xlim(len(Ez))						
    axes[0].set_ylim(-3,3)
    axes[0].set_title("E Field")
    axes[1].set_xlim(len(Hy))
    axes[1].set_ylim(-3,3)
    axes[1].set_title("H Field")
    line, = axes[0].plot([],[])						# Get line objects of the axes in the plot in order to modify it
    line1, = axes[1].plot([],[])					# from animate function.
    def init():
        """ Initializes the plots to have zero values."""
        line.set_data([],[])
        line1.set_data([],[])
        return line,
    
    if (source_choice==1):						# When the source is sinusoidal of frequency 1Hz.
        def animate(n, *args, **kwargs):
            """ Animation function that sets the value of line objects as required.
                  Argument n refers to the frame number: time in our case. """
            # Update Equations.
            Hy[0:xdim-1] = Hy[0:xdim-1]+(deltat/(delta*mu[0:xdim-1]))*(Ez[1:xdim]-Ez[0:xdim-1])
            Ez[1:xdim]= Ez[1:xdim]+(deltat/(delta*epsilon[1:xdim]))*(Hy[1:xdim]-Hy[0:xdim-1])
            Ez[xsource] = Ez[xsource] + np.sin(2*n*np.pi/180)		# Sinusoidal source.
            ylims = axes[0].get_ylim()					
            if (abs(np.amax(Ez))>ylims[1]):				# In case the value of Fields become more than that of scale,   
	        axes[0].set_ylim(-(np.amax(Ez)+2),np.amax(Ez)+2)	# sets the scale accordingly.
	        axes[1].set_ylim(-(np.amax(Ez)+2),np.amax(Hy)+2)
            line.set_data(np.arange(len(Ez)),Ez)			# Sets the updated values of fields.
            line1.set_data(np.arange(len(Hy)),Hy)
            #fig.suptitle('Time = %s'%str(n))
	    return line,
        ani = animation.FuncAnimation(fig, animate, init_func=init,	# Note: We assign a reference to the animation, otherwise the   
             frames=(time_tot), interval=10, blit=False, repeat =False)	# animation object will be considered for garbage collection.
        fig.show()


    elif (source_choice==2):						# When the source is an impulse of height 1 unit.
        axes[0].set_ylim(-0.1,0.1)
        axes[1].set_ylim(-0.1,0.1)
        Ez[xsource]=1.0							# Impulse source.
        def animate(n, *args, **kwargs):
            """ Same as above."""
            Hy[0:xdim-1] = Hy[0:xdim-1]+(deltat/(delta*mu[0:xdim-1]))*(Ez[1:xdim]-Ez[0:xdim-1])
            Ez[1:xdim]= Ez[1:xdim]+(deltat/(delta*epsilon[1:xdim]))*(Hy[1:xdim]-Hy[0:xdim-1])
            ylims = axes[0].get_ylim()
            if (abs(np.amax(Ez))>ylims[1]):				# Scaling axes.
	        axes[0].set_ylim(-(np.amax(Ez)+2),np.amax(Ez)+2)
	        axes[1].set_ylim(-(np.amax(Ez)+2),np.amax(Hy)+2)
            line.set_data(np.arange(len(Ez)),Ez)
            line1.set_data(np.arange(len(Hy)),Hy)
	    if n==0: Ez[xsource]=0					# Impulse source:end of impulse.
            return line,
        ani = animation.FuncAnimation(fig, animate, init_func=init, frames=(time_tot), interval=10, blit=False, repeat =False)
        fig.show()

    elif (source_choice==3):						# When the source is gaussian.
        def animate(n, *args, **kwargs):
            """Same as earlier."""
            Hy[0:xdim-1] = Hy[0:xdim-1] + (deltat/(delta*mu[0:xdim-1]))*(Ez[1:xdim]-Ez[0:xdim-1])
            Ez[1:xdim]= Ez[1:xdim] + (deltat/(delta*epsilon[1:xdim]))*(Hy[1:xdim]-Hy[0:xdim-1])
            Ez[xsource] = Ez[xsource] + 30.0*(1/np.sqrt(2*np.pi))*np.exp(-(n-80.0)**2/(100))	# Gaussian source.
            ylims = axes[0].get_ylim()
            if (abs(np.amax(Ez))>ylims[1]):				# Scaling the axes.
	        axes[0].set_ylim(-(np.amax(Ez)+2),np.amax(Ez)+2)
	        axes[1].set_ylim(-(np.amax(Ez)+2),np.amax(Hy)+2)
	    line.set_data(np.arange(len(Ez)),Ez)
            line1.set_data(np.arange(len(Hy)),Hy)
	    return line,
        ani = animation.FuncAnimation(fig, animate, init_func=init, frames=(time_tot), interval=5, blit=False, repeat =False)
        fig.show()
    else : 
        print (" The choice entered is not valid. Please enter a valid choice and try again.\n")

if __name__ == "__main__": main()
