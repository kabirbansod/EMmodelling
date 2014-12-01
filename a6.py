# Unitless 1D FDTD solution to wave equation. 
# Author- Kabir Bansod
# Date: 10/09/14


import numpy as np							# Loading required libraries.
import matplotlib.pyplot as plt
import time
import matplotlib.animation as animation

class FdTd:
    def __init__(self):
        self.time_tot = 500						# Time for which the simulation runs.
    
    # Stability factor
        self.S = 1						
    
    #Speed of light
        self.c = 299792458						

    #defining dimensions							
        self.n_x = 300
        self.n_xsrc = self.n_x/6							# Location of the source.

        self.epsilon0 = 8.854187817*10**(-12)						
        self.mu0 = 4*np.pi*10**(-7)
    
        self.sigma = np.ones(self.n_x)*4*10**(-4)
        self.sigma_star = np.ones(self.n_x)*4*10**(-4)
        
        self.get_epsilon()

        self.source_choice = self.get_source()					# Get the type of source to be used from the user.

        if (self.source_choice==1):
            fr = input("Enter the frequency of the source in THz:")
            self.freq = fr*10**(12)
            min_wavelength = self.c/(self.freq*np.sqrt(np.amax(self.epsilon)/self.epsilon0))
            #print min_wavelength
            #print np.sqrt(np.amax(epsilon))*freq
            self.delta = min_wavelength/30
        else:
            self.delta = 10**(-6)

        self.deltat = self.S*self.delta/self.c							# Time step.
        print self.deltat

        self.h_frequency = 0.5/(self.deltat)
        self.number_steps_freq = 1001

        self.ar = np.zeros(self.n_x)
        self.dum1 = 0
        self.dum2 = 0
        return

    def get_epsilon(self):
        """ Returns an array containing the permittivity values of the medium""" 
        a =input("Enter Refractive index:")
        self.epsilon=np.ones(self.n_x)*self.epsilon0
        self.epsilon[self.n_x/3:2*self.n_x/3] = self.epsilon[self.n_x/3:2*self.n_x/3]*a
        self.mu = self.mu0*np.ones(self.n_x)
        self.mu[self.n_x/3:2*self.n_x/3] = self.mu[self.n_x/3:2*self.n_x/3]*2        
        return
    
    def get_source(self):
        """Gets the source to be used in the calculations from 
           the user"""
        print ("Pick a source to be used in 1D FDTD simulation of wave prpogation in vacuum.")
        print ("1: Sinusoidal source of frequency 1Hz.")
        print ("2: Impulse.")
        print ("3: Gaussian source.")
        return input("Enter your choice: ") 


    def plotFields(self):
        xdim = self.n_x
        xsrc = self.n_xsrc

        Ez = np.zeros(xdim)							# Arrays that hold the values of Electric field 
        Hy = np.zeros(xdim)								# and magnetic field respectively.

        A = ((self.mu-0.5*self.deltat*self.sigma_star)/(self.mu+0.5*self.deltat*self.sigma_star))
        B=(self.deltat/self.delta)/(self.mu+0.5*self.deltat*self.sigma_star)
    
        C=((self.epsilon-0.5*self.deltat*self.sigma)/(self.epsilon+0.5*self.deltat*self.sigma)) 
        D=(self.deltat/self.delta)/(self.epsilon+0.5*self.deltat*self.sigma)

        H_src_correction = -np.sqrt((self.epsilon[xsrc]/self.epsilon0)/(self.mu[xsrc]/self.mu0))
        n_src = np.sqrt(self.epsilon[xsrc]/self.epsilon0)
        delay_half_grid = n_src*self.delta/(2*self.c)
        half_time_step = self.deltat/2

        h_f = self.h_frequency
        n_steps = self.number_steps_freq
        frequencies = np.linspace(0,h_f,n_steps)
        
        exp_terms = np.exp(-1j*2*np.pi*self.deltat*frequencies)
        ref = np.ones(n_steps)*(0+0j)
        trn = np.ones(n_steps)*(0+0j)
        src = np.ones(n_steps)*(0+0j)

        fig , axes = plt.subplots(2,1)					# Figure with 2 subplots- one for E field and one for H.
        axes[0].set_xlim(len(Ez))						
        axes[0].set_title("E Field")
        axes[1].set_xlim(len(Hy))
        axes[1].set_title("H Field")
        line, = axes[0].plot([],[])					# Get line objects of the axes in the plot in order to modify it
        line1, = axes[1].plot([],[])					# from animate function.
        def init():
            """ Initializes the plots to have zero values."""
            line.set_data([],[])
            line1.set_data([],[])
            return line,

        if (self.source_choice==1):						# When the source is sinusoidal of frequency 1Hz.
            def sine_src(n):
                source = np.sin(2*np.pi*self.freq*self.deltat*n)
                return source
            def animate(n, *args, **kwargs):
                """ Animation function that sets the value of line objects as required.
                      Argument n refers to the frame number: time step in our case. """
                Esrc = sine_src(n)
                Hsrc = H_src_correction*sine_src(n + delay_half_grid + half_time_step)
                Ez[xsrc] = Ez[xsrc] + Esrc					# Sinusoidal source.
                Hy[0:xdim-1] = A[0:xdim-1]*Hy[0:xdim-1]+B[0:xdim-1]*(Ez[1:xdim]-Ez[0:xdim-1])
                Hy[xsrc-1] = Hy[xsrc-1] - B[xsrc-1]*Esrc
                Ez[1:xdim]= C[1:xdim]*Ez[1:xdim]+D[1:xdim]*(Hy[1:xdim]-Hy[0:xdim-1])
                #Ez[xsrc] = Ez[xsrc] - D[xsrc]*Hsrc
                ylims1 = axes[0].get_ylim()
                ylims2 = axes[1].get_ylim()
                e_max = abs(np.amax(Ez))
                h_max = abs(np.amax(Hy))
                Ez[xdim-1]= self.dum2+((self.S-1)/(self.S+1))*(Ez[xdim-2]-Ez[xdim-1])
                Ez[1]= self.dum1+((self.S-1)/(self.S+1))*(Ez[1]-Ez[0])
                Ez[0]= Ez[2]
                self.dum1 = Ez[2]
                self.dum2 = Ez[xdim-2]
              
                ref[0:n_steps] = ref[0:n_steps] + (exp_terms[0:n_steps]**(n))*Ez[0]
                trn[0:n_steps] = trn[0:n_steps] + (exp_terms[0:n_steps]**(n))*Ez[xdim-1]
                src[0:n_steps] = src[0:n_steps] + (exp_terms[0:n_steps]**(n))*Esrc

                if (e_max > ylims1[1]):				# Scaling axes.
	            axes[0].set_ylim(-(1.2*e_max),1.2*e_max)
                if ((h_max)>ylims2[1]):				# Scaling axes.
	            axes[1].set_ylim(-(1.2*h_max),1.2*h_max)
                line.set_data(np.arange(len(Ez)),Ez)			# Sets the updated values of fields.
                line1.set_data(np.arange(len(Hy)),Hy)
                #time.sleep(0.01)
	        return line,
            ani = animation.FuncAnimation(fig, animate, init_func=init, frames= self.time_tot, interval=10, blit=False, repeat =False)	
            # Note: We assign a reference to the animation, otherwise the                
            # animation object will be considered for garbage collection.
            fig.show()
            return ref, trn, src
    
        elif (self.source_choice==2):						# When the source is an impulse of height 1 unit.
            #axes[0].set_ylim(-0.1,0.1)
            axes[1].set_ylim(-0.002,0.002)
            def animate(n, *args, **kwargs):
                """ Same as above."""
                if (n==0):
                    Esrc = 1
                else: 
                    Esrc = 0
                Hsrc = H_src_correction*0
                Ez[xsrc] = Ez[xsrc] + Esrc
                Hy[0:xdim-1] = A[0:xdim-1]*Hy[0:xdim-1]+B[0:xdim-1]*(Ez[1:xdim]-Ez[0:xdim-1])
                Hy[xsrc-1] = Hy[xsrc-1] - B[xsrc-1]*Esrc
                Ez[1:xdim]= C[1:xdim]*Ez[1:xdim]+D[1:xdim]*(Hy[1:xdim]-Hy[0:xdim-1])
                Ez[xsrc] = Ez[xsrc] - D[xsrc]*Hsrc
                ylims1 = axes[0].get_ylim()
                ylims2 = axes[1].get_ylim()
                e_max = abs(np.amax(Ez))
                h_max = abs(np.amax(Hy))
                Ez[xdim-1]= self.dum2+((self.S-1)/(self.S+1))*(Ez[xdim-2]-Ez[xdim-1])
                Ez[1]= self.dum1+((self.S-1)/(self.S+1))*(Ez[1]-Ez[0])
                Ez[0]= Ez[2]
                self.dum1 = Ez[2]
                self.dum2 = Ez[xdim-2]

                ref[0:n_steps] = ref[0:n_steps] + (exp_terms[0:n_steps]**(n))*Ez[0]
                trn[0:n_steps] = trn[0:n_steps] + (exp_terms[0:n_steps]**(n))*Ez[xdim-1]
                src[0:n_steps] = src[0:n_steps] + (exp_terms[0:n_steps]**(n))*Esrc

                if (e_max > ylims1[1]):					# Scaling axes.
	            axes[0].set_ylim(-(1.2*e_max),1.2*e_max)
                if ((h_max)>ylims2[1]):					# Scaling axes.
	            axes[1].set_ylim(-(1.2*h_max),1.2*h_max)
                line.set_data(np.arange(len(Ez)),Ez)
                line1.set_data(np.arange(len(Hy)),Hy)
	        if n==0: 
                    Ez[xsrc] = Ez[xsrc] + 0				# Impulse source.
                return line,
            ani = animation.FuncAnimation(fig, animate, init_func=init, frames=(self.time_tot), interval=10, blit=False, repeat =False)
            fig.show()
            return ref, trn, src

        elif (self.source_choice==3):						# When the source is gaussian.
            def gauss_src(n):
                source = 1.0*(1/np.sqrt(2*np.pi))*np.exp(-((n-80.0)*self.deltat)**2/(5*self.deltat)**2)
                return source
            def animate(n, *args, **kwargs):
                """Same as earlier."""
                Esrc = gauss_src(n)
                Hsrc = H_src_correction*Esrc
                Ez[xsrc] = Ez[xsrc] + Esrc
                Hy[0:xdim-1] = A[0:xdim-1]*Hy[0:xdim-1]+B[0:xdim-1]*(Ez[1:xdim]-Ez[0:xdim-1])
                Hy[xsrc-1] = Hy[xsrc-1] - B[xsrc-1]*Esrc
                Ez[1:xdim]= C[1:xdim]*Ez[1:xdim]+D[1:xdim]*(Hy[1:xdim]-Hy[0:xdim-1])
                Ez[xsrc] = Ez[xsrc] - D[xsrc]*Hsrc
                ylims1 = axes[0].get_ylim()
                ylims2 = axes[1].get_ylim()
                e_max = abs(np.amax(Ez))
                h_max = abs(np.amax(Hy))
                Ez[xdim-1]= self.dum2+((self.S-1)/(self.S+1))*(Ez[xdim-2]-Ez[xdim-1])
                Ez[1]= self.dum1+((self.S-1)/(self.S+1))*(Ez[1]-Ez[0])
                Ez[0]= Ez[2]
                self.dum1 = Ez[2]
                self.dum2 = Ez[xdim-2]

                ref[0:n_steps] = ref[0:n_steps] + (exp_terms[0:n_steps]**(n))*Ez[0]
                trn[0:n_steps] = trn[0:n_steps] + (exp_terms[0:n_steps]**(n))*Ez[xdim-1]
                src[0:n_steps] = src[0:n_steps] + (exp_terms[0:n_steps]**(n))*Esrc

                if (e_max > ylims1[1]):				# Scaling axes.
	            axes[0].set_ylim(-(1.2*e_max),1.2*e_max)
                if ((h_max)>ylims2[1]):				# Scaling axes.
	            axes[1].set_ylim(-(1.2*h_max),1.2*h_max)
                line.set_data(np.arange(len(Ez)),Ez)
                line1.set_data(np.arange(len(Hy)),Hy)
	        return line,
            ani = animation.FuncAnimation(fig, animate, init_func=init, frames=(self.time_tot), interval=5, blit=False, repeat =False)
            fig.show()
            return ref, trn, src

        else : 
            print (" The choice entered is not valid. Please enter a valid choice and try again.\n")
    

def main():
    simulation = FdTd()
    ref ,trn, src = simulation.plotFields()
    fig2 , ax = plt.subplots(3,1)
    frequencies = np.linspace(1,simulation.h_frequency,simulation.number_steps_freq)
    ref = (ref)*simulation.deltat
    trn = (trn)*simulation.deltat
    src = (src)*simulation.deltat
    ref = (abs(ref/src))**2
    trn = (abs(trn/src))**2
    
    #print ref[1:10]
    #print np.amax(trn)
    ax[0].plot(frequencies, ref.real )
    ax[0].set_title("FFT of relection coefficient")
    ax[1].plot(frequencies, trn.real )
    ax[1].set_title("FFT of Transmission coefficient")
    ax[2].plot(frequencies, src.real)
    ax[2].set_title("FFT of Source")
    plt.show()
    fig3, axis = plt.subplots(1,1)
    axis.plot(frequencies,ref+trn)
    #axis.set_ylim(0,2)
    plt.show()
    
    
if __name__ == "__main__": main()
