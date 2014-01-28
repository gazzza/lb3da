#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.integrate as sp
from scipy.interpolate import spline
import lb3da.matrix.mdf

Lx=128
Ly=96
Lz=96
Rpar=5

def PeriodFilter(array):
    

#    array[:,0] = np.where(array[:,0] >= 0.5*Lx, array[:,0] - Lx, array[:,0])
 
    array[:,1] = np.where(array[:,1] >= 0.5*Ly, array[:,1] - Ly , array[:,1])
    array[:,2] = np.where(array[:,2] >= 0.5*Lz, array[:,2] - Lz , array[:,2])
    return array

def CoordListFromFile(filename):
    """Reads lines from MD file and returns numpy array of coordinate tuples"""
    coordinates  = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if len(line) > 0:
                line = line.split()
                (x,y,z) = (map(np.float,(line[0],line[1],line[2]))) 
                coordinates.append((x,y,z))
    return np.array(coordinates)


def RawDists(array,matrix):
    """Subtracts the particle's coordinates from the other particle's coords"""
    raw_dists = array - matrix
    return np.abs(raw_dists)

def Distance(matrix):
    """Returns a numpy array of distances along each array row"""
    return  np.sqrt(matrix[:,0]**2 + matrix[:,1]**2 + matrix[:,2]**2)


def PairCorrelationFunction(filename):
    
    particles = mdf.MDMatrix(filename)
    allcoords = particles.x()

    init_shell = 15
    
    totalN = particles.shape[0]
    print 'Total N:' + str(totalN)

    coords = allcoords[np.where(Distance(allcoords - 64) < init_shell)]
#    coords = allcoords

    N = coords.shape[0]
    print 'Number of considered particles: ' + str(N)


    rho = np.float(totalN)/(Lx*Ly*Lz) 
    #print 'Particle number density: ' + str(rho)

    dr = 1.0
    maxR = int(np.ceil(min(Lx/2,Ly/2,Lz/2) - dr - init_shell)) ##Ensure that the radius considered does not extend beyond the system boundaries
    print 'Max R: ' + str(maxR)

    rvals  = [] #radius values list
    pvals  = [] #g(r) values list

    for r in np.nditer(np.arange(0,maxR+1,dr)):                                             #Iterate over an array of floats rather than ints
        total_particles_in_shell = 0                                                        #Initialise the number of particles in a shell
        for i in range(N):                                                                  #Loop over chosen particles
            arrayslice                = coords[i]                                           #Get the coordinates of the i'th particle
            raw_dists                 = RawDists(arrayslice,allcoords)                      #Find x,y,z distance between the considered particle and all others 
            filtered_dists            = PeriodicFilter(raw_dists)                           #Filter the distances to take into account periodic boundaries   
            dists                     = Distance(filtered_dists)                            #Find absolute distance between the considered particle and all others
            dists                     = dists[np.where(dists != 0.0)]                       #Ignore the particle self-distance. 
            shell                     = dists[(np.where((dists >= r) & (dists <= r + dr)))] #Find all particles in shell of between r and r + dr
            number_in_shell           = shell.size                                          #Number of particles in the shell
            total_particles_in_shell += number_in_shell                                     #Accumulate the number of particles in the shell

        
        shell_vol = (4*np.pi/3.)*((r + dr)**3 - r**3)                           #Volume of the spherical shell
        normalised_by_number_of_p        = total_particles_in_shell/np.float(N) #Divide the total number of particles in a shell by the number of particles consider                                                                                  ed to get the average number of particles in a shell
        normalised_by_shell_volume       = normalised_by_number_of_p/shell_vol  #Divide the average number of particles in a shell by the volume of the shell. """
        normalised_by_number_density     = normalised_by_shell_volume/rho       #Divide the average number of particles per unit volume of a shell by the number densi                                                                                ty of particles in your simulation domain

        # print 'Average number in shell/shell volume/density: ' + str(normalised_by_number_density)
        # print 'Average number in shell/shell volume: '         + str(normalised_by_shell_volume)
        # print 'Average number in shell: '                      + str(normalised_by_number_of_p)
        
        rvals.append(r)                             #Add the current r value to the r array
        pvals.append(normalised_by_number_density)  #Add the corresponding current g(r) value to the array


    return (rvals,pvals)
    
    plt.show()
    
if __name__ == "__main__":

    pcfavg = []
    mdlist = MDListMaker(1,400000)
    for files in mdlist[-100:]:
        pcf = PairCorrelationFunction(files)
        pcfavg.append(pcf[1])
    pcfavg = np.array(pcfavg)
    avgpcf = np.mean(pcfavg,axis = 0)
    
    r = np.array(pcf[0])/(2*Rpar)
#    fit = FittingFunc(pcf[0],avgpcf,10)
#    plt.plot(pcf[0],avgpcf,'b-')
#    plt.plot(fit[0],fit[1],'b-')
    
    
    #Spline plot
    r_new = np.linspace(0.8,np.max(r),1000)
    gr_smooth = spline(r,avgpcf,r_new)
    plt.plot(r,avgpcf,'ro') #Data
    plt.plot(r_new, gr_smooth) #Fit
    plt.yticks([1,np.floor(np.max(avgpcf)),np.ceil(np.max(avgpcf))])
    plt.ylim([0,np.ceil(np.max(avgpcf))])
    plt.xlabel(r'$r/D_{p}$',fontsize=20)
    plt.ylabel(r'$g(r)$',fontsize=20)
    plt.title(r'Radial Distribution Function - 1.4% minority fluid, $u_{s} = 0.001$')
    plt.savefig('/home/gary/Dropbox/PRL_paper/pcf_fr0.01_u0.001.pdf',format='pdf',dpi=1000)
    plt.show()




    
