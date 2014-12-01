# EE6506: Computational Electromagnetics
# Author- Kabir Bansod
# EE11B023
# Date: 10/09/14

# This program simulates the propagation of an electromagnetic wave in a 
# strip waveguide. The simulation is performed using Finite Difference 
# Time Domain form of Maxwell's Equations. The boundary used in the program
# is Mur's boundary. The simulation can be run for 2 or 3 strips. one strip 
# having the source while the other two having detectors.

import numpy as np							                            # Loading required libraries.
import matplotlib.pyplot as plt
import sys
import matplotlib.animation as animation


class Strip_wg:                                                         # Waveguide object.
    def __init__(self, a=2):
        self.xdim = 150                                                 # Dimensions.
        self.ydim = 150

        self.time_tot = 500                                            # Total time for which the simulation runs.

        self.xsource = self.xdim/2                                      # Location of source (X & Y co-ordinates).
        self.ysource = self.ydim-10
    
        self.S = 1/2**0.5                                               # Courant Factor.

        self.epsilon0 = 8.854187817*10**(-12)                           # Permittivity of free space.
        self.mu0 = 4*np.pi*10**(-7)                                     # Permeability of free space. 
        self.c = 299792458.0                                            # Velocity of light in vacuum.

        self.epsilonr=2.25                                              # Relative permeability of waveguide medium.

        self.delta = 10**(-7)                                           # Space step.
        self.deltat = self.delta*self.S/self.c                          # Time step.

        self.sigma =  4*10**(-4)*np.ones([self.xdim,self.ydim])         # Electric conductivity.
        self.sigma_star = 4*10**(-4)*np.ones([self.xdim,self.ydim])     # Magnetic conductivity.
 
        self.t = input("Enter the width of the strips: ")               # Thickness of strips.
        self.w = 1+input("Enter the separation between strips: ")       # Separation between the strips.

        # Location of strips.
        self.x_main_strip = np.arange(self.xsource-self.t/2,self.xsource+self.t/2)
        self.x_2nd_strip = np.arange(self.xsource+self.t/2+self.w,self.xsource+self.t/2+self.w+self.t)
        self.x_3rd_strip = np.arange(self.xsource-self.t/2-self.w-self.t,self.xsource-self.t/2-self.w)

        # Frequencies for which the detector calculates the Fourier Transforms.
        self.n_freq = 100
        self.wavelengths = np.arange(1530,1570)*10**-9
        self.frequencies= np.linspace(-(self.c/np.amin(self.wavelengths)),(self.c/np.amin(self.wavelengths)), self.n_freq)                             #self.c/self.wavelengths

        self.three_strps = False                                        # If three waveguides are to be used then this value becomes true.
        if (a==3):self.three_strps = True
        elif(a==2): pass
        else:                                                           # The number of strips should either be two or three.
            try:  
                raise NumberofStripsExceeded(a)
            except NumberofStripsExceeded as nbr:
                print 'The number of waveguides can be either 2 or 3, you entered:', nbr.value
                sys.exit(1)

        self.get_epsilon()                                              # Method that creates the medium (permittivity matrix).
        return
    
    
    def get_epsilon(self):
        """ Function that creates permittivity and permeability 
          array based on the thickness of strips and separation between them"""
        self.epsilon = np.ones([self.xdim,self.ydim])*self.epsilon0     
        self.mu = np.ones([self.xdim,self.ydim])*self.mu0
        self.epsilon[self.x_main_strip,:] =self.epsilon0*self.epsilonr 
        self.epsilon[self.x_2nd_strip,0:2*self.ydim/3]=self.epsilon0*self.epsilonr
        if (self.three_strps==True):
            self.epsilon[self.x_3rd_strip,0:2*self.ydim/3]=self.epsilon0*self.epsilonr

        return

    def show_strips(self):
        """ Displays the permittivity matrix of the medium. """
        plt.imshow(self.epsilon.T)
        plt.show()

    def set_thickness_and_separation(self, thickness, separation, relative_permittivity=2.25):
        """ This function changes the thickness of the strips and the separation between them. """
        self.t = thickness
        self.w = separation
        self.epsilonr = relative_permittivity 
        self.get_epsilon()
        return



    def plot_fields(self, simulation):
        """ This function simulates the propagation of EM waves in the waveguide. """
        
        # Creating Arrays having the update coeeficients so that we don have to calculate them again and again in the loop.
        A = ((self.mu-0.5*self.deltat*self.sigma_star)/(self.mu+0.5*self.deltat*self.sigma_star)) 
        B = (self.deltat/self.delta)/(self.mu+0.5*self.deltat*self.sigma_star)

        C = ((self.epsilon-0.5*self.deltat*self.sigma)/(self.epsilon+0.5*self.deltat*self.sigma)); 
        D = (self.deltat/self.delta)/(self.epsilon+0.5*self.deltat*self.sigma)  
        
        # Creating local variable so that we dont hav to access the object variables in the loop.
        n_x = self.xdim
        n_y = self.ydim
        xsrc = self.x_main_strip
        ysrc = self.ysource
        x_strp2 = self.x_2nd_strip

        number_of_frames = self.time_tot
        dt = self.deltat
        
        n_freq = self.n_freq

        frequency=self.frequencies
        exp_terms = np.exp(-1j*2*np.pi*self.deltat*frequency)
        
        # 2D Arrays that hold values of fields.
        Ez = np.zeros([n_x,n_y])
        Hy = np.zeros([n_x,n_y])
        Hx = np.zeros([n_x,n_y])

        strip1 = np.zeros(number_of_frames)   #np.ones((self.t,n_freq))*(0+0j)
        strip2 = np.zeros(number_of_frames)   #np.ones((self.t,n_freq))*(0+0j)
        source = np.zeros(number_of_frames)   #np.ones((self.t,n_freq))*(0+0j)

        p0 = 1
        p2=-0.5

        c0=(self.c/(2*self.S))*(1-(p0/self.S))
        c1=-(self.c/(2*self.S))*(1+(p0/self.S))
        c2=(self.c/(self.S**2))*(p0+(p2*(self.S**2)))
        c3=-(p2*self.c)/2
        c0efffor=-(c0/c1)
        c2efffor=-(c2/c1)
        c3efffor=-(c3/c1)
        c0=(self.c/(2*self.S))*(1+(p0/self.S))
        c1=-(self.c/(2*self.S))*(1-(p0/self.S))
        c2=-(self.c/(self.S**2))*(p0+(p2*(self.S**2)))
        c3=(p2*self.c)/2
        c1effrev=-(c1/c0)
        c2effrev=-(c2/c0)
        c3effrev=-(c3/c0)

        prev_xfor=np.zeros((1,n_y))
        prev_prev_xfor=np.zeros(prev_xfor.shape)
        prev_x_minus_1for=np.zeros((1,n_y))
        prev_prev_x_minus_1for=np.zeros(prev_x_minus_1for.shape)
        prev_yfor=np.zeros((n_x,1))
        prev_prev_yfor=np.zeros(prev_yfor.shape)
        prev_y_minus_1for=np.zeros((n_x,1))
        prev_prev_y_minus_1for=np.zeros(prev_y_minus_1for.shape)
        prev_xrev=np.zeros((1,n_y))
        prev_prev_xrev=np.zeros(prev_xrev.shape)
        prev_x_minus_1rev=np.zeros((1,n_y))
        prev_prev_x_minus_1rev=np.zeros(prev_x_minus_1rev.shape)
        prev_yrev=np.zeros((n_x,1))
        prev_prev_yrev=np.zeros(prev_yrev.shape)
        prev_y_minus_1rev=np.zeros((n_x,1))
        prev_prev_y_minus_1rev=np.zeros(prev_y_minus_1rev.shape)

        def init(p):
            global fig, ax, im
            fig = plt.figure()
            ax = plt.axes()
            im = ax.imshow(Ez, vmin=-p, vmax=+p)
            return

        def animate(n, *args, **kwargs):
            """ Function that updates the contour plots/frames. 
                     n refers to the frame number--time in our case"""
            Hx[1:n_x-3,1:n_y-3]=A[1:n_x-3,1:n_y-3]*Hx[1:n_x-3,1:n_y-3]-B[1:n_x-3,1:n_y-3]*(Ez[1:n_x-3,2:n_y-2]-Ez[1:n_x-3,1:n_y-3])
            Hy[1:n_x-3,1:n_y-3]=A[1:n_x-3,1:n_y-3]*Hy[1:n_x-3,1:n_y-3]+B[1:n_x-3,1:n_y-3]*(Ez[2:n_x-2,1:n_y-3]-Ez[1:n_x-3,1:n_y-3])

            Ez[2:n_x-3,2:n_y-3]=C[2:n_x-3,2:n_y-3]*Ez[2:n_x-3,2:n_y-3]+(Hy[2:n_x-3,2:n_y-3]-Hy[1:n_x-4,2:n_y-3]-Hx[2:n_x-3,2:n_y-3]+Hx[2:n_x-3,1:n_y-4])*D[2:n_x-3,2:n_y-3]
            
            Ez[n_x-3,2:n_y-3]=c0efffor*(Ez[n_x-4,2:n_y-3]+prev_prev_xfor[0,2:n_y-3])-prev_prev_x_minus_1for[0,2:n_y-3]+c2efffor*(prev_xfor[0,2:n_y-3]+prev_x_minus_1for[0,2:n_y-3])+c3efffor*(prev_x_minus_1for[0,1:n_y-4]+prev_x_minus_1for[0,3:n_y-2]+prev_xfor[0,1:n_y-4]+prev_xfor[0,3:n_y-2])

            prev_prev_xfor[:,:]=prev_xfor[:,:]
            prev_prev_x_minus_1for[:,:]=prev_x_minus_1for[:,:]
            prev_xfor[0,0:n_y]=Ez[n_x-3,0:n_y]
            prev_x_minus_1for[0,0:n_y]=Ez[n_x-4,0:n_y]

            Ez[1,2:n_y-3]=-prev_prev_xrev[0,2:n_y-3]+c1effrev*(Ez[2,2:n_y-3]+prev_prev_x_minus_1rev[0,2:n_y-3])+c2effrev*(prev_xrev[0,2:n_y-3]+prev_x_minus_1rev[0,2:n_y-3])+c3effrev*(prev_x_minus_1rev[0,1:n_y-4]+prev_x_minus_1rev[0,3:n_y-2]+prev_xrev[0,1:n_y-4]+prev_xrev[0,3:n_y-2])

            prev_prev_xrev[:,:]=prev_xrev[:,:]
            prev_prev_x_minus_1rev[:,:]=prev_x_minus_1rev[:,:]
            prev_xrev[0,0:n_y]=Ez[2,0:n_y]
            prev_x_minus_1rev[0,0:n_y]=Ez[1,0:n_y]

            Ez[2:n_x-3,n_y-3]=c0efffor*(Ez[2:n_x-3,n_y-4]+prev_prev_yfor[2:n_x-3,0])-prev_prev_y_minus_1for[2:n_x-3,0]+c2efffor*(prev_yfor[2:n_x-3,0]+prev_y_minus_1for[2:n_x-3,0])+c3efffor*(prev_y_minus_1for[1:n_x-4,0]+prev_y_minus_1for[3:n_x-2,0]+prev_yfor[1:n_x-4,0]+prev_yfor[3:n_x-2,0])

            prev_prev_yfor[:,:]=prev_yfor[:,:]
            prev_prev_y_minus_1for[:,:]=prev_y_minus_1for[:,:]
            prev_yfor[0:n_x,0]=Ez[0:n_x,n_y-3]
            prev_y_minus_1for[0:n_x,0]=Ez[0:n_x,n_y-4]

            Ez[2:n_x-3,1]=-prev_prev_yrev[2:n_x-3,0]+c1effrev*(Ez[2:n_x-3,2]+prev_prev_y_minus_1rev[2:n_x-3,0])+c2effrev*(prev_yrev[2:n_x-3,0]+prev_y_minus_1rev[2:n_x-3,0])+c3effrev*(prev_y_minus_1rev[1:n_x-4,0]+prev_y_minus_1rev[3:n_x-2,0]+prev_yrev[1:n_x-4,0]+prev_yrev[3:n_x-2,0])

            prev_prev_yrev[:,:]=prev_yrev[:,:]
            prev_prev_y_minus_1rev[:,:]=prev_y_minus_1rev[:,:]
            prev_yrev[0:n_x,0]=Ez[0:n_x,2]
            prev_y_minus_1rev[0:n_x,0]=Ez[0:n_x,1]

            Ez[1,1]=prev_prev_xrev[0,2]
            Ez[1,n_y-3]=prev_prev_xrev[0,n_y-4]
            Ez[n_x-3,1]=prev_prev_x_minus_1for[0,2]
            Ez[n_x-3,n_y-3]=prev_prev_x_minus_1for[0,n_y-4]

            #for m in range(self.t):
            #    strip1[m,0:n_freq] = strip1[m,0:n_freq] + (exp_terms[0:n_freq]**(n))*(Ez[xsrc[m],n_y-3])
            #    strip2[m,0:n_freq] = strip2[m,0:n_freq] + (exp_terms[0:n_freq]**(n))*(Ez[x_strp2[m],n_y-3])
            #    source[m,0:n_freq] = source[m,0:n_freq] + (exp_terms[0:n_freq]**(n))*(Ez[xsrc[m],ysrc])
            strip1[n] = Ez[xsrc[self.t/2],n_y-3]
            strip2[n] = Ez[x_strp2[self.t/2],n_y-3]
            source[n] = Ez[xsrc[self.t/2],ysrc]


            gauss_src= 1.0*np.exp(-(n-80.0)**2/(5**2))
            Ez[xsrc,ysrc]=gauss_src           #np.sin(2*np.pi*frequency[i]*n*dt)
            if(simulation==True):im.set_data(Ez.T)
            return

        if (simulation==True):
            init(0.9)
            ani = animation.FuncAnimation(fig, animate, frames= number_of_frames, interval=15, blit=False, repeat =False)	
            fig.show()
                # Note: We assign a reference to the animation, otherwise the                
                # animation object will be considered for garbage collection.
        #else:

        #strip1=strip1*self.deltat
        #strip2=strip2*self.deltat
        #source=source*self.deltat
        #strip1=(abs(strip1/source))**2
        #strip2=(abs(strip2/source))**2
      
        return strip1, strip2, source

class NumberofStripsExceeded(Exception):
    def __init__(self, num):
        self.value = num
    
    def __str__(self):
        return repr(self.value)

def main():
    #a = input("Number of strips to be used: ")
    a = 2                    
    s_wg = Strip_wg(a)                                          # Create an object of simulation parameters.
    #s_wg.set_thickness_and_separation(10,10)
    s_wg.show_strips()                                          # Display the arrangement of strips.
    str1, str2, src = s_wg.plot_fields(True)                    # Store the arrays having detector values returned after simulation. 
    #main_strip = np.ones(s_wg.n_freq)*(0+0j)
    #second_strip = np.ones(s_wg.n_freq)*(0+0j)
    #source_strip = np.ones(s_wg.n_freq)*(0+0j)
    #for i in range(s_wg.n_freq):                                # Taking the average along the thickness of the strips.
    #    for j in range(s_wg.t):
    #        main_strip[i]=main_strip[i] + str1[j,i]
    #        second_strip[i] = second_strip[i]+ str2[j,i]
    #        source_strip[i]=source_strip[i] + src[j,i]

    #main_strip=main_strip/s_wg.t
    #second_strip=second_strip/s_wg.t
    #source_strip=source_strip/s_wg.t
    
    fig2 , ax = plt.subplots(3,1)
    ax[0].plot(np.arange(s_wg.time_tot), str1 )
    ax[0].set_title("FFT of 1st strip")
    ax[1].plot(np.arange(s_wg.time_tot), str2 )
    ax[1].set_title("FFT of 2nd strip")
    ax[2].plot(np.arange(s_wg.time_tot), src)
    ax[2].set_title("FFT of Source")
    plt.show()
    

if __name__=="__main__": main()
