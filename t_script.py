#!/usr/bin/env python

import lb3da.matrix.mdfiles as mdf
import lb3da.files.outputlist as opf

import os
import matplotlib.pyplot as plt

output_dir = opf.OutputDir('test_directory/')
mdlist = [os.path.join(output_dir.path(),file) for file in output_dir.list('md-cfg','.asc')]


t_array = []
for file in mdlist:
    print file
    f = mdf.MDFile(file)
    m = mdf.MDMatrix(f.create_matrix())
    m.md_outputs(['force','torque'])
    print m.f()
    t_array.append(m.t(0)[2])

print t_array

