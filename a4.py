# A time stepping leap frog solution to 1D Maxwell's equations 
# Author- Kabir Bansod
# Date: 10/09/14


import numpy as np							# Loading required libraries.
import matplotlib.pyplot as plt
import time
import matplotlib.animation as animation

def get_source():							
    """Gets the source to be used in the calculations from 
       the user"""
    print ("Pick a source to be used in 1D FDTD simulation of wave prpogation in vacuum.")
    print ("1: Sinusoidal source of frequency 1Hz.")
    print ("2: Impulse.")
    print ("3: Gaussian source.")
    return input("Enter your choice: ") 

def get_epsilon(epsilon0,x_dim):
    """ Returns an array containing the permittivity values of the medium""" 
    epsilon=np.ones(x_dim)*epsilon0
    return epsilon


def main():
    #defining dimensions							
    
    time_tot = 1000							# Time for which the simulation runs.
    #xsource = xdim/2							# Location of the source.
    
    # Stability factor
    S = 1						
    
    #Speed of light
    c = 299792458.0	
    epsilon0 = 8.854187817*10**(-12)					# Values of permittivity and permeability.	
    mu0 = 4*np.pi*10**(-7)

    n_x = 300								# Dimension(Length) of the simulation.
    n_xsrc = n_x/2							# Location of the source.
    
    Ez = np.zeros(n_x)							# Arrays that hold the values of Electric field and 
    Hy = np.zeros(n_x)							# magnetic field respectively.

    epsilon = get_epsilon(epsilon0,n_x)					# Medium charecteristics.
    mu = mu0*np.ones(n_x)
    sigma = np.ones(n_x)*4*10**(-4)
    sigma_star = np.ones(n_x)*4*10**(-4)


    source_choice = get_source()					# Get the type of source to be used from the user.

    if (source_choice==1):
        fr = input("Enter the frequency of the source in THz:")
        freq = float(fr*10**(12))
        min_wavelength = c/(freq*np.sqrt(np.amax(epsilon)/epsilon0))
        #print min_wavelength
        #print np.sqrt(np.amax(epsilon))*freq
        delta = min_wavelength/30
    else:
        delta = 10**(-6)

    deltat = S*delta/c							# Time step.

    A = ((mu-0.5*deltat*sigma_star)/(mu+0.5*deltat*sigma_star))		# Arrays that hold update coefficients, so that we dont need to
    B=(deltat/delta)/(mu+0.5*deltat*sigma_star)				# calculate them in every step of the loop.
    
    C=((epsilon-0.5*deltat*sigma)/(epsilon+0.5*deltat*sigma)) 
    D=(deltat/delta)/(epsilon+0.5*deltat*sigma)

    #epsilon[3*n_x/4] =100						# Impurity (point impurity)
    fig , axes = plt.subplots(2,1)					# Figure with 2 subplots- one for E field and one for H.
    axes[0].set_xlim(len(Ez))						
    #axes[0].set_ylim(-3,3)
    axes[0].set_title("E Field")
    axes[1].set_xlim(len(Hy))
    #axes[1].set_ylim(-0.001,0.001)
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
            Hy[0:n_x-1] = A[0:n_x-1]*Hy[0:n_x-1]+B[0:n_x-1]*(Ez[1:n_x]-Ez[0:n_x-1])
            Ez[1:n_x]= C[1:n_x]*Ez[1:n_x]+D[1:n_x]*(Hy[1:n_x]-Hy[0:n_x-1])
            Ez[n_xsrc] = Ez[n_xsrc] + np.sin(2.0*np.pi*freq*deltat*float(n))				# Sinusoidal source.
            Ez[0]=0						# Perfect Electric conductor (Dirichlet condition).
            Hy[n_x-1]=0						# Perfect Magnetic conductor (Dirichlet condition).
            ylims1 = axes[0].get_ylim()				# Gets the limits of axes to scale them.
            ylims2 = axes[1].get_ylim()
            e_max = abs(np.amax(Ez))
            h_max = abs(np.amax(Hy))
            if (e_max > ylims1[1]):				# Scaling axes.
	        axes[0].set_ylim(-(1.2*e_max),1.2*e_max)
            if ((h_max)>ylims2[1]):				# Scaling axes.
	        axes[1].set_ylim(-(1.2*h_max),1.2*h_max)
            line.set_data(np.arange(len(Ez)),Ez)			# Sets the updated values of fields.
            line1.set_data(np.arange(len(Hy)),Hy)
            #time.sleep(0.01)
	    return line,
        ani = animation.FuncAnimation(fig, animate, init_func=init, frames= time_tot, interval=10, blit=False, repeat =False)	
        # Note: We assign a reference to the animation, otherwise the                
        # animation object will be considered for garbage collection.
        fig.show()


    elif (source_choice==2):						# When the source is an impulse of height 1 unit.
        #axes[0].set_ylim(-0.1,0.1)
        axes[1].set_ylim(-0.001,0.001)
        Ez[n_xsrc]=1.0							# Impulse source.
        def animate(n, *args, **kwargs):
            """ Same as above."""
            Hy[0:n_x-1] = A[0:n_x-1]*Hy[0:n_x-1]+B[0:n_x-1]*(Ez[1:n_x]-Ez[0:n_x-1])
            Ez[1:n_x]= C[1:n_x]*Ez[1:n_x]+D[1:n_x]*(Hy[1:n_x]-Hy[0:n_x-1])
            Ez[0]=0						# Perfect Electric conductor (Dirichlet condition).
            Hy[n_x-1]=0						# Perfect Magnetic conductor (Dirichlet condition).
            ylims1 = axes[0].get_ylim()
            ylims2 = axes[1].get_ylim()
            e_max = abs(np.amax(Ez))
            h_max = abs(np.amax(Hy))
            if (e_max > ylims1[1]):				# Scaling axes.
	        axes[0].set_ylim(-(1.2*e_max),1.2*e_max)
            if ((h_max)>ylims2[1]):				# Scaling axes.
	        axes[1].set_ylim(-(1.2*h_max),1.2*h_max)
            line.set_data(np.arange(len(Ez)),Ez)
            line1.set_data(np.arange(len(Hy)),Hy)
	    if n==0: Ez[n_xsrc]=0					# Impulse source:end.
            return line,
        ani = animation.FuncAnimation(fig, animate, init_func=init, frames=(time_tot), interval=10, blit=False, repeat =False)
        fig.show()

    elif (source_choice==3):						# When the source is gaussian.
        def animate(n, *args, **kwargs):
            """Same as earlier."""
            Hy[0:n_x-1] = A[0:n_x-1]*Hy[0:n_x-1]+B[0:n_x-1]*(Ez[1:n_x]-Ez[0:n_x-1])
            Ez[1:n_x]= C[1:n_x]*Ez[1:n_x]+D[1:n_x]*(Hy[1:n_x]-Hy[0:n_x-1])
            Ez[n_xsrc] = Ez[n_xsrc] + 40.0*(1/np.sqrt(2*np.pi))*np.exp(-((n-80.0)*deltat)**2/(5*deltat)**2)
            # ^Gaussian source (soft source).
            Ez[0]=0						# Perfect Electric conductor (Dirichlet condition).
            Hy[n_x-1]=0						# Perfect Magnetic conductor (Dirichlet condition).
            ylims1 = axes[0].get_ylim()
            ylims2 = axes[1].get_ylim()
            e_max = abs(np.amax(Ez))
            h_max = abs(np.amax(Hy))
            if (e_max > ylims1[1]):				# Scaling axes.
	        axes[0].set_ylim(-(1.2*e_max),1.2*e_max)
            if ((h_max)>ylims2[1]):				# Scaling axes.
	        axes[1].set_ylim(-(1.2*h_max),1.2*h_max)
	    line.set_data(np.arange(len(Ez)),Ez)
            line1.set_data(np.arange(len(Hy)),Hy)
	    return line,
        ani = animation.FuncAnimation(fig, animate, init_func=init, frames=(time_tot), interval=5, blit=False, repeat =False)
        fig.show()
    else : 
        print (" The choice entered is not valid. Please enter a valid choice and try again.\n")

# Note that we are not updating values of E fields at the 0th position and H field at n_x position, so they will have 

if __name__ == "__main__": main()
