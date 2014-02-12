#!/usr/bin/env python

import nuro.files

def test_MDFile_create_matrix():

    ref_particle_1_list = [36.7,64.0,108.5,
                           0.0,0.0,0.0,
                           -0.0059716544,-4.4408921e-16,0.000874300658,
                           0.0,1.0,2.22044605e-16,
                           0.0,0.0,0.0,
                           1.02418074e-14,-0.00012134029,-2.84217094e-14,
                           1.0]
    ref_particle_2_list = [36.7,64.0,148.5,
                           0.0,0.0,0.0,
                           -0.0059716544,-8.8817842e-16,-0.000874300658,
                           0.0,1.0,2.22044605e-16,
                           0.0,0.0,0.0,
                           -8.29891711e-15,0.00012134029,-7.46069873e-14,
                           2.0]

    fob = nuro.files.MDFile('two_particles_interface.asc')
    m = fob.create_matrix()
    
    assert len(m) == 2
    assert ref_particle_1_list == m[0]
    assert ref_particle_2_list == m[1]
    assert len(m[0]) == 19
    assert len(m[1]) == 19

