#!/usr/bin/env python

import numpy as np

import nuro.matrix
import nuro.files

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
m = nuro.matrix.MDMatrix(fob.create_matrix(),['force','torque'])

def test_MDFile_x():
    ref_x = np.array([[36.7,64.0,108.5],[36.7,64.0,148.5]])
    ref_x1 = np.array([36.7,64.0,108.5])
    ref_x2 = np.array([36.7,64.0,148.5])

    x = m.x()
    x1 = m.x(0)
    x2 = m.x(1)
    
    "x tests"
    assert x.shape == (2,3)
    assert x.shape == ref_x.shape
    assert ref_x[0].all() == x[0].all()
    assert ref_x[1].all() == x[1].all()
    
    "x1 tests"
    assert x1.shape == (3,)
    assert x1.shape == ref_x1.shape
    assert x1.all() == ref_x1.all()

    "x2 tests"
    assert x2.shape == (3,)
    assert x2.shape == ref_x2.shape
    assert x2.all() == ref_x2.all()
    
    

def test_MDFile_v():
    ref_v  = np.array([[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0]])
    ref_v1 = np.array([0.0,0.0,0.0,0.0])
    ref_v2 = np.array([0.0,0.0,0.0,0.0])

    v = m.v()
    v1 = m.v(0)
    v2 = m.v(1)
    
    "v tests"
    assert v.shape == (2,4)
    assert v.shape == ref_v.shape
    assert ref_v[0].all() == v[0].all()
    assert ref_v[1].all() == v[1].all()
    
    "v1 tests"
    assert v1.shape == (4,)
    assert v1.shape == ref_v1.shape
    assert v1.all() == ref_v1.all()

    "v2 tests"
    assert v2.shape == (4,)
    assert v2.shape == ref_v2.shape
    assert v2.all() == ref_v2.all()

def test_MDFile_f():
    ref_f  = np.array([[-0.0059716544,-4.4408921e-16,0.000874300658,((-0.0059716544)**2 + (-4.4408921e-16)**2 + (0.000874300658**2))**0.5]
                       ,[-0.0059716544,-8.8817842e-16,-0.000874300658,((-0.0059716544)**2 + (-8.8817842e-16)**2 + (-0.000874300658**2))**0.5]])
    ref_f1 = np.array([-0.0059716544,-4.4408921e-16,0.000874300658,((-0.0059716544)**2 + (-4.4408921e-16)**2 + (0.000874300658)**2)**0.5])
    ref_f2 = np.array([-0.0059716544,-8.8817842e-16,-0.000874300658,((-0.0059716544)**2 + (-8.8817842e-16)**2 + (-0.000874300658)**2)**0.5])

    f = m.f()
    f1 = m.f(0)
    f2 = m.f(1)
    
    "f tests"
    assert f.shape == (2,4)
    assert f.shape == ref_f.shape
    assert ref_f[0].all() == f[0].all()
    assert ref_f[1].all() == f[1].all()
    
    "f1 tests"
    assert f1.shape == (4,)
    assert f1.shape == ref_f1.shape
    assert f1.all() == ref_f1.all()

    "f2 tests"
    assert f2.shape == (4,)
    assert f2.shape == ref_f2.shape
    assert f2.all() == ref_f2.all()

def test_MDFile_o():
    ref_o  = np.array([[0.0,1.0,2.22044605e-16],
                       [0.0,1.0,2.22044605e-16]])
    ref_o1 = np.array([0.0,1.0,2.22044605e-16])
    ref_o2 = np.array([0.0,1.0,2.22044605e-16])

    o = m.o()
    o1 = m.o(0)
    o2 = m.o(1)
    
    "o tests"
    assert o.shape == (2,3)
    assert o.shape == ref_o.shape
    assert ref_o[0].all() == o[0].all()
    assert ref_o[1].all() == o[1].all()
    
    "o1 tests"
    assert o1.shape == (3,)
    assert o1.shape == ref_o1.shape
    assert o1.all() == ref_o1.all()

    "o2 tests"
    assert o2.shape == (3,)
    assert o2.shape == ref_o2.shape
    assert o2.all() == ref_o2.all()


def test_MDFile_av():
    ref_av  = np.array([[0.0,0.0,0.0,0.0],
                        [0.0,0.0,0.0,0.0]])
    ref_av1 = np.array([0.0,0.0,0.0,0.0])
    ref_av2 = np.array([0.0,0.0,0.0,0.0])

    av = m.av()
    av1 = m.av(0)
    av2 = m.av(1)
    
    "av tests"
    assert av.shape == (2,4)
    assert av.shape == ref_av.shape
    assert ref_av[0].all() == av[0].all()
    assert ref_av[1].all() == av[1].all()
    
    "av1 tests"
    assert av1.shape == (4,)
    assert av1.shape == ref_av1.shape
    assert av1.all() == ref_av1.all()

    "av2 tests"
    assert av2.shape == (4,)
    assert av2.shape == ref_av2.shape
    assert av2.all() == ref_av2.all()
    
def test_MDFile_t():
    ref_t  = np.array([[1.02418074e-14,-0.00012134029,-2.84217094e-14,((1.02418074e-14)**2 + (-0.00012134029)**2 + (-2.84217094e-14)**2)**0.5],
                       [-8.29891711e-15,0.00012134029,-7.46069873e-14,((-8.29891711e-15)**2 + (0.00012134029)**2 + (-7.46069873e-14)**2)**0.5]])
    ref_t1 = np.array([1.02418074e-14,-0.00012134029,-2.84217094e-14,((1.02418074e-14)**2 + (-0.00012134029)**2 + (-2.84217094e-14)**2)**0.5])
    ref_t2 = np.array([-8.29891711e-15,0.00012134029,-7.46069873e-14,((-8.29891711e-15)**2 + (0.00012134029)**2 + (-7.46069873e-14)**2)**0.5])

    t = m.t()
    t1 = m.t(0)
    t2 = m.t(1)
    
    "t tests"
    assert t.shape == (2,4)
    assert t.shape == ref_t.shape
    assert ref_t[0].all() == t[0].all()
    assert ref_t[1].all() == t[1].all()
    
    "t1 tests"
    assert t1.shape == (4,)
    assert t1.shape == ref_t1.shape
    assert t1.all() == ref_t1.all()

    "t2 tests"
    assert t2.shape == (4,)
    assert t2.shape == ref_t2.shape
    assert t2.all() == ref_t2.all()
    
def test_MDFile_id():
    ref_pid = np.array([[1.0],[2.0]])
    ref_pid_1 = np.array([1.0])
    ref_pid_2 = np.array([2.0])

    pid = m.id()
    pid_1 = m.id(0)
    pid_2 = m.id(1)

    assert pid.shape == (2,1)
    assert pid.shape == ref_pid.shape
    assert ref_pid[0].all() == pid[0].all()
    assert ref_pid[1].all() == pid[1].all()
