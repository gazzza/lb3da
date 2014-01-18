#!/usr/bin/env python

import os

class OutputDir(object):
    """Create an output directory object e.g. output_dir = OutputDir('path-to-file')"""
    def __init__(self,output_dir):
        self.dir = output_dir

    def path(self):
        """Returns the full path of the output directory"""
        return os.path.abspath(self.dir)
    
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
            
    def filter(self,prefix,file_ext):
        """Returns a list of files which match prefix and suffix"""
        return [file for file in self.list if (prefix in file.split('_')[0] and os.path.splitext(file)[1] == file_ext)]

    def filter_ext(self,file_ext):
        """Returns a list of files which match the file extension"""
        return [file for file in self.list if (os.path.splitext(file)[1] == file_ext)]
