#!/usr/bin/env python

import os

class OutputDir(object):
    """Create an output directory object e.g. output_dir = OutputDir('path-to-file')"""
    def __init__(self,output_dir):
        self.dir = output_dir

    def path(self,string=None):
        """Returns the full path of the output directory"""
        full_path = os.path.abspath(self.dir)
        if string is not None:
            return os.path.join(full_path,string)
        else:
            return full_path
    
    def list(self,prefix=None,file_ext=None):
        """Returns a list of the contents of the output directory"""
        if prefix is None and file_ext is None:
            return os.listdir(self.path())
        elif prefix is not None and file_ext is None:
            return [file for file in self.list() if prefix in file.split('_')[0]]
        elif prefix is None and file_ext is not None:
            return [file for file in self.list() if (os.path.splitext(file)[1] == file_ext)]
        else: 
            return [file for file in self.list() if (prefix in file.split('_')[0] and os.path.splitext(file)[1] == file_ext)]

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
