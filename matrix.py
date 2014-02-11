#!/usr/bin/env python

import os
import numpy as np
import re
import pdb

from collections import OrderedDict


        
class MDMatrix(object):
    """Create an MDMatrix object e.g. matrix = MDMatrix(MDfile.create_matrix)"""
    _default_outputs = ['position','velocity','orientation','angular_velocity','id']
    _all_outputs = OrderedDict([('position',3),('velocity',3),('force',3),('orientation_q',4),
                               ('orientation',3),('angular_velocity',3),('torque',3),('id',1)])    
    def __init__(self,matrix,outputs_list=[]):
        self.outputs = sorted(self._default_outputs + outputs_list,key=self._all_outputs.keys().index)
        self.m = np.array(matrix)
        self.N = self.m.shape[0]
        self._indices_dict = self._indices()

    def _indices(self):
        """Calculates where along each row vector the relevant quantities 
        lie and returns a dict of the indices of each quantity"""
        index = 0
        indices_dict = {}
        for quantity in self.outputs:
            newindex = index + self._all_outputs[quantity]
            indices_dict[quantity] = (index,newindex)
            index = newindex
        return indices_dict

    def _4vector(self,quantity,pid=None):
        """Extracts 3 array values depending on the quantity 
        defined and additionally returns their norm"""
        try:
            i,j = self._indices_dict[quantity]
        except KeyError:
            print "You didn't set " + quantity + "in the output parameters"
            return None
            
        if pid is not None:
            quant = self.m[pid,i:j]
            return np.hstack((quant,np.linalg.norm(quant)))
        else:
            quant = self.m[:,i:j]    
            return np.hstack((quant,np.reshape(np.linalg.norm(quant,axis=1),(self.N,1))))

    def _3vector(self,quantity,pid=None):
        """Extracts 3 array values according to their precedence 
        in the all_outputs list and the quantities defined""" 
        try:
            i,j = self._indices_dict[quantity]
        except KeyError:
            print "You didn't set " + quantity + " in the output parameters"
            return None
        

        if pid is not None:
            return self.m[pid,i:j]
        else:
            return self.m[:,i:j]

    def x(self,pid=None):
        """Returns the x,y,z positions of the particle(s)
        
        Parameters
        ----------
        pid: integer, optional
        particle of interest - corresponds to rows in the MD file matrix
        
        Returns
        -------
        out: ndarray"""

        return self._3vector('position',pid)
    

    def v(self,pid=None):
        """Returns the vx,vy,vz,sqrt(vx^2+vy^2+vz^2) of the particle(s)
        
        Parameters
        ----------
        pid: integer, optional
        particle of interest - corresponds to rows in the MD file matrix
        
        Returns
        -------
        out: ndarray"""

        return self._4vector('velocity',pid)
    
    def f(self,pid=None):
        """Returns the fx,fy,fz,sqrt(fx^2+fy^2+fz^2) of the particle(s)
        
        Parameters
        ----------
        pid: integer, optional
        particle of interest - corresponds to rows in the MD file matrix
        
        Returns
        -------
        out: ndarray"""

        return self._4vector('force',pid)
    
    def t(self,pid=None):
        """Returns the tx,ty,tz,sqrt(tx^2+ty^2+tz^2) of the particle(s)
        
        Parameters
        ----------
        pid: integer, optional
        particle of interest - corresponds to rows in the MD file matrix
        
        Returns
        -------
        out: ndarray"""

        return self._4vector('torque',pid)
    
    def o(self,pid=None):
        """Returns the ox,oy,oz of the particle(s)
        
        Parameters
        ----------
        pid: integer, optional
        particle of interest - corresponds to rows in the MD file matrix
        
        Returns
        -------
        out: ndarray"""

        return self._3vector('orientation',pid)
    def oq(self,pid=None):
        """Returns the orientation in quaternion notation, o1,o2,o3,o4 of the particle(s)
        
        Parameters
        ----------
        pid: integer, optional
        particle of interest - corresponds to rows in the MD file matrix
        
        Returns
        -------
        out: ndarray"""

        return self._3vector('orientation_q',pid)
    
    def av(self,pid=None):
        """Returns the avx,avy,avz,sqrt(avx^2+avy^2+avz^2) of the particle(s)
        
        Parameters
        ----------
        pid: integer, optional
        particle of interest - corresponds to rows in the MD file matrix
        
        Returns
        -------
        out: ndarray"""

        return self._4vector('angular_velocity',pid)
    
    def id(self,pid=None):
        """Returns the particle(s) id(s). Remember that python arrays start from 0 and fortran arrays start from 1. 
        Hence, particle 1 in the LB3D MD output file is particle 0 in this analysis programme.
        
        Parameters
        ----------
        pid: integer, optional
        particle of interest - corresponds to rows in the MD file matrix
        
        Returns
        -------
        out: ndarray"""

        return self._3vector('id',pid)

    def _get_angle_matrix(self,pid=None):
        """This function simply divides each element of the orientation matrix
        by the norm of its corresponding row in the orientation matrix"""
        orientation_matrix = self.o(pid)
        if pid is not None:
            return orientation_matrix/np.linalg.norm(orientation_matrix)
        else:
            return orientation_matrix/np.reshape(np.linalg.norm(orientation_matrix,axis=1),(self.N,1))
    
    def angle(self,pid=None,radians=False):
        """Returns the orientation of the particle with respect to a particular axis i.e. 
        the dot product of the particle's orientation vector and the corresponding axis"""
        angle_matrix = self._get_angle_matrix(pid)
        cos_angle_matrix = np.arccos(angle_matrix)

        if radians:
            return cos_angle_matrix
        else:
            return np.rad2deg(cos_angle_matrix)

    def S(self):
        """Returns the S orientational order parameter, S=0.5*<3cos^2x -1> where x 
        is the angle between the INTERFACE NORMAL and the particle symmetry axis"""
        angle_matrix = self.angle(radians=True)
        s_avg_part = 3*np.cos(angle_matrix)**2 - 1
        s_avg = np.mean(s_avg_part,axis=0)
        full_s = 0.5*s_avg
        return full_s

    def Q(self):
        """Returns the Q orientational order parameter, Q=1.5*<sin^2x cos2y> where x 
        is the angle between the INTERFACE NORMAL and the particle symmetry axis"""
        angle_matrix = self.angle(radians=True)
        sint = np.sin(angle_matrix)**2
        cost = np.cos(2*angle_matrix)

        full_Q = np.dot(sint.T,cost)/sint.shape[0]
        
        return 1.5*full_Q


    def displacements(self,x=None):
        """This function will return the displacements of all the particles
        and a coordinate"""
        if x is not None:
            return self.x() - np.array(x)
        else:
            return self.x() - self.x(0)

    def distances(self,x=None):
        """This function will return the absolute distances between 
        the particles and a specified coordinate"""

        displacement_matrix = self.displacements(x)

        return np.reshape(np.linalg.norm(displacement_matrix,axis=1),(self.N,1))

    @staticmethod
    def _periodic_filter(array,dims,periodicity):
        """This function takes into account the periodic boundary conditions of the system.
        If the boundary is periodic, then whenever the particle distance is greater than half 
        the system dimension in that axis direction, we must subtract half the system size 
        in that direction to account for the fact that the periodic images are closer"""
        
        Lx,Ly,Lz = dims
        
        if periodicity[0]:
            array[:,0] = np.where(array[:,0] >= 0.5*Lx, array[:,0] - Lx, array[:,0])
        
        if periodicity[1]:
            array[:,1] = np.where(array[:,1] >= 0.5*Ly, array[:,1] - Ly, array[:,1])
        
        if periodicity[2]:
            array[:,2] = np.where(array[:,2] >= 0.5*Lz, array[:,2] - Lz, array[:,2])

        return array

    def pair_correlation_function(self,dims=(128,96,96),periodicity=(False,True,True),shell_r=15,dr=1.0):
        """Calculates the pair correlation function (pcf). This pcf differs from most in the literature 
        because it can handle non-periodic boundaries and finite-sized particles. 

        The algorithm should be used as follows:
        
        1) Select some particles close to the centre of the system by choosing shell_r parameter. 
        E.g. shell_r = 15 will find all particles within a distance of 15 of the centre of the system.

        2) Choose how many concentric shells you want by choosing dr. 
        Make sure you take into account your particle sizes when selecting dr. 
        E.g. Does choosing dr < 1.0 give meaningful results if your particles are radius 5?
        
        3) Enter the system size using the dim parameter and which boundaries are periodic by using the 
        periodicity parameter. E.g. periodicity=(False,True,False) sets periodic boundaries in the y direction
        and non-periodic boundaries in the x and z directions. 
        
        The algorithm then outputs the r values and corresponding 
        g(r) values in a tuple of numpy arrays for your pleasure
        
        Notes:
        
        1) The Max R parameter sets the largest concentric shell radius 
        that can be considered, which depends on your initial choice of shell_r.
        Hence, a shell_r as small as possible is preferable, but make sure you have 
        enough particles included in shell_r to give reasonable statistics. 

        2) The pcf, or g(r) values, should tend to 1 as r goes to infinity. However, if you 
        have even a single non-periodic boundary this will not necessarily be the case. """


        Lx,Ly,Lz = dims
        maxR = int(np.ceil(min(Lx/2,Ly/2,Lz/2) - dr - shell_r)) 
        print 'Max Shell Radius: ' + str(maxR)
        
        all_particle_coords = self.x()
        
        """This is some numpy trickery to select the rows of the allcoords matrix which match the conditions of bool_array"""
        bool_array = self.distances((Lx/2,Ly/2,Lz/2)) < shell_r
        bool_array = bool_array.reshape(self.N)
        particles_within_shell_r = all_particle_coords[bool_array,:]

        N = particles_within_shell_r.shape[0]
        print 'Number of particles selected by shell_r: ' + str(N)

        rvals  = [] 
        pvals  = [] 

        for r in np.nditer(np.arange(0,maxR+1,dr)):
            total_particles_in_shell = 0           
            for p in particles_within_shell_r:        
                "Compute the displacements between particle p and all other particles"
                displacements             = self.displacements(p)
                "Take into account the periodic boundaries"
                filtered_displacements    = MDMatrix._periodic_filter(displacements,dims,periodicity)
                dists                     = np.linalg.norm(filtered_displacements,axis=1)
                """Remove particle self-distance. Not necessary because it will never 
                be included in the shell variable, but just making it explicit"""
                dists                     = dists[np.where(dists != 0.0)]
                shell                     = dists[(np.where((dists >= r) & (dists <= r + dr)))]
                number_in_current_shell   = shell.size                                         
                total_particles_in_shell += number_in_current_shell                                    

            "This version of shell_vol should be used because it is EXACT. Do not use the approximate (4*np.pi/3.)*r**2*dr."
            shell_vol = (4*np.pi/3.)*((r + dr)**3 - r**3)                           
            
            "Normalise by number of particles expected in an ideal gas"
            normalised_by_no_of_particles = total_particles_in_shell/np.float(N) 
            normalised_by_shell_volume    = normalised_by_no_of_particles/shell_vol  
            normalised_by_number_density  = normalised_by_shell_volume/(np.float(self.N)/(Lx*Ly*Lz))

        
            rvals.append(r)                             
            pvals.append(normalised_by_number_density) 


        return (np.array(rvals),np.array(pvals))
    

    def orientation_correlation_function(self,dims=(96,256,256),periodicity=(False,True,True),shell_r=15,dr=1.0):
        """Calculates the orientation correlation function (ocf). This ocf differs from most in the literature 
        because it can handle non-periodic boundaries and finite-sized particles. 

        The algorithm should be used as follows:
        
        1) Select some particles close to the centre of the system by choosing shell_r parameter. 
        E.g. shell_r = 15 will find all particles within a distance of 15 of the centre of the system.

        2) Choose how many concentric shells you want by choosing dr. 
        Make sure you take into account your particle sizes when selecting dr. 
        E.g. Does choosing dr < 1.0 give meaningful results if your particles are radius 5?
        
        3) Enter the system size using the dim parameter and which boundaries are periodic by using the 
        periodicity parameter. E.g. periodicity=(False,True,False) sets periodic boundaries in the y direction
        and non-periodic boundaries in the x and z directions. 
        
        The algorithm then outputs the r values and corresponding 
        o(r) values in a tuple of numpy arrays for your pleasure
        
        Notes:
        
        1) The Max R parameter sets the largest concentric shell radius 
        that can be considered, which depends on your initial choice of shell_r.
        Hence, a shell_r as small as possible is preferable, but make sure you have 
        enough particles included in shell_r to give reasonable statistics. 

        Algorithm overview:

        1) Finds ID's of particles within a distance shell_r of the system centre.
        These particles are called reference particles.
        2) For each of these reference particles, calculates the angle difference between
        the reference particle and every other particle within a distance r.
        3) Averages the cosine of these differences
        4) Averages (3) over all the reference particles."""

        Lx,Ly,Lz = dims
        maxR = int(np.ceil(min(Ly/2,Lz/2) - dr - shell_r)) 
        print 'Max Shell Radius: ' + str(maxR)

        ids_reference_particles = np.where(self.distances((Lx/2,Ly/2,Lz/2)) < shell_r)[0]
        print "ID's of reference particles: " + str(ids_reference_particles)

        rvals = []
        ovals = []

        for r in np.nditer(np.arange(0,maxR+1,dr)):
    
            ref_p_cosines = []
            for ref_p in ids_reference_particles:
                ids_particles_within_r = np.where(self.distances(self.x(ref_p)) < r)[0]
                no_particles_within_r = ids_particles_within_r.shape[0]
                if no_particles_within_r <= 1:
                    break
        
                cosines_array = []
                for p in ids_particles_within_r:
                    cosine = np.cos(self.angle(ref_p,radians=True) - self.angle(p,radians=True))
                    if not p == ref_p:
                        cosines_array.append(cosine)
                cosines_avg = np.mean(cosines_array,axis=0)
                ref_p_cosines.append(cosines_avg)
        
            if not len(ref_p_cosines) == 0:
                ref_p_cosines_avg = np.mean(ref_p_cosines,axis=0)
                rvals.append(r)
                ovals.append(ref_p_cosines_avg)
            else:
                pass


        return (np.array(rvals),np.array(ovals))



    
    



    
