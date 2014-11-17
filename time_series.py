#!/usr/bin/env python

import os
import numpy as np
import pdb

import matrix


class TimeSeries(object):
    def __init__(self,dir_ob):
        self.d = dir_ob
        self.ts = self._construct_time_series_object()
        
    def _construct_time_series_object(self):
        md_list = self.d.list('md-cfg','.asc')
        
        tsteps = []
        objects = []
        for f in md_list:
            f_ob = MDFile(self.d.path() + '/' + f)
            
            tsteps.append(f_ob.tstep())
            objects.append(matrix.MDMatrix(f_ob.create_matrix()))
        
        return (tsteps,objects)
