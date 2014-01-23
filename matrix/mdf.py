#!/usr/bin/env python

import os
import numpy as np
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


class MDMixin(object):
    default_outputs = ['position','velocity','orientation','angular_velocity','id']
    all_outputs = OrderedDict([('position',3),('velocity',3),('force',3),('orientation_q',4),
                               ('orientation',3),('angular_velocity',3),('torque',3),('id',1)])

    def __init__(self,matrix):
        self.outputs = self.default_outputs
        
    def md_outputs(self,outputs_list=[]):
        """Returns the list of quantities, in order, used in the output file. If no additional quantities are defined
        the function returns the default ordered list as per MDMatrix.all_outputs"""
        if 'reset' in outputs_list:
            self.outputs = self.default_outputs
        else:
            self.outputs = sorted(self.default_outputs + outputs_list,key=self.all_outputs.keys().index)
        return self.outputs

class MDMatrix(object):
    """Create an MDMatrix object e.g. matrix = MDMatrix(MDfile.create_matrix)"""
    
    def __init__(self,matrix):
        super(MDMixin,self).__init__()
        self.shape = self.m.shape
        self.m = 

    def md_outputs(self,outputs_list=[]):
        """Returns the list of quantities, in order, used in the output file. If no additional quantities are defined
        the function returns the default ordered list as per MDMatrix.all_outputs"""
        if 'reset' in outputs_list:
            self.outputs = self.default_outputs
        else:
            self.outputs = sorted(self.default_outputs + outputs_list,key=self.all_outputs.keys().index)
        return self.outputs

    def _indices(self):
        """Calculates where along each row vector the relevant quantities 
        lie and returns a dict of the indices of each quantity"""
        index = 0
        indices_dict = {}
        for quantity in self.outputs:
            newindex = index + self.all_outputs[quantity]
            indices_dict[quantity] = (index,newindex)
            index = newindex
        return indices_dict

    def _4vector(self,quantity,pid=None):
        """Extracts 3 array values depending on the quantity 
        defined and additionally returns their norm"""
        try:
            i,j = self._indices()[quantity]
        except KeyError:
            print "You didn't set " + quantity + "in the output parameters"
            return None
            
        if pid is not None:
            quant = self.m[pid,i:j]
            return np.hstack((quant,np.linalg.norm(quant)))
        else:
            quant = self.m[:,i:j]    
            return np.hstack((quant,np.reshape(np.linalg.norm(quant,axis=1),(self.shape[0],1))))

    def _3vector(self,quantity,pid=None):
        """Extracts 3 array values according to their precedence 
        in the all_outputs list and the quantities defined""" 
        try:
            i,j = self._indices()[quantity]
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



        






