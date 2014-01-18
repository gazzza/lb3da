#!/usr/bin/env

import mdfiles


f = mdfiles.MDFile('md-cfg_output_t00003500-1413193586.asc')


m = mdfiles.MDMatrix(f.create_matrix())

m.md_outputs(['force','torque'])
