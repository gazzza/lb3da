#!/usr/bin/env python

import os
import numpy as np
import re

from collections import OrderedDict

class MDFile(object):
    """Create an MD file object e.g. bla = MDFile('md-filename-here')"""

    def __init__(self,mdfile):
        self.f = mdfile
        
    def path(self):
        """Returns the full path of the output directory"""
        return os.path.abspath(self.f)
    
    def create_matrix(self,pid=None):
        """Reads lines from a file and stores as a matrix"""
        matrix  = []
        f = open(self.f, 'r')
        for line in f:
            line = line.strip()
            if len(line) > 0:
                try:
                    matrix.append(map(np.float, line.split()))
                except ValueError:
                    pass
        f.close()
        if pid is not None:
            return matrix[pid]
        else:
            return matrix        
    
    def tstep(self):
        """Strips the tstep out of the filename and returns it as an int"""
        fname = os.path.split(self.f)[1]
        str = re.search('t\d\d\d\d\d\d\d\d-\d\d\d\d\d\d\d\d\d\d',fname).group()
        splitstr = str.split('-')[0]
        strip_zeros = splitstr.lstrip('t0')
        return int(strip_zeros)
        
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
    def _periodic_filter(array):
        Lx,Ly,Lz = (128,96,96)
    #    array[:,0] = np.where(array[:,0] >= 0.5*Lx, array[:,0] - Lx, array[:,0])
    
        array[:,1] = np.where(array[:,1] >= 0.5*Ly, array[:,1] - Ly , array[:,1])
        array[:,2] = np.where(array[:,2] >= 0.5*Lz, array[:,2] - Lz , array[:,2])

        return array

    def pair_correlation_function(self):
        allcoords = self.x()
        Lx,Ly,Lz = (128,96,96)
        shell_r = 15
        rho = np.float(self.N)/(Lx*Ly*Lz) 
        dr = 1.0
        maxR = int(np.ceil(min(Lx/2,Ly/2,Lz/2) - dr - shell_r)) 
        print 'Max R: ' + str(maxR)
        print 'Total N:' + str(self.N)

        bool_array = self.distances((64,64,64)) < shell_r
        bool_array = bool_array.reshape(self.N)
        coords = allcoords[bool_array,:]

        N = coords.shape[0]
        print 'Number of considered particles: ' + str(N)

        rvals  = [] 
        pvals  = [] 

        for r in np.nditer(np.arange(0,maxR+1,dr)):
            total_particles_in_shell = 0           
            for i in coords:                       
                raw_dists                 = self.displacements(i)
                filtered_dists            = MDMatrix._periodic_filter(raw_dists)
                dists                     = np.linalg.norm(filtered_dists,axis=1)
                dists                     = dists[np.where(dists != 0.0)]
                shell                     = dists[(np.where((dists >= r) & (dists <= r + dr)))]
                number_in_shell           = shell.size                                         
                total_particles_in_shell += number_in_shell                                    

        
            shell_vol = (4*np.pi/3.)*((r + dr)**3 - r**3)                           
            normalised_by_number_of_p        = total_particles_in_shell/np.float(N) 
            normalised_by_shell_volume       = normalised_by_number_of_p/shell_vol  
            normalised_by_number_density     = normalised_by_shell_volume/rho       

        
            rvals.append(r)                             
            pvals.append(normalised_by_number_density) 


        return (np.array(rvals),np.array(pvals))
    
        



    
