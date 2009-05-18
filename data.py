#! /usr/bin/env python

# Copyright 2009 Thomas Neumann
#
# Redistribution of this file is permitted under
# the terms of the GNU Public License (GPL) version 2.

from enthought.tvtk.api import tvtk
import numpy as N
from model import Joint, Element, Construction, ElementMaterial

# ----------- 3D models ------------
# completely free joint: represented by a sphere
joint_model_movable = tvtk.SphereSource(phi_resolution=24, theta_resolution=24, radius=0.15).output
# unmovable joint: represented by a solid cube
joint_model_unmovable = tvtk.CubeSource(x_length=0.7, y_length=0.7, z_length=0.5).output
# ground in x direction
collect = tvtk.AppendPolyData()
collect.add_input(tvtk.CubeSource(x_length=0.7, y_length=0.2, z_length=0.25, center=(0, -0.5, 0)).output)
collect.add_input(tvtk.SphereSource(radius=0.11, center=(-0.18, -0.3, 0), phi_resolution=24, theta_resolution=24).output)
collect.add_input(tvtk.SphereSource(radius=0.11, center=(0.18, -0.3, 0), phi_resolution=24, theta_resolution=24).output)
collect.add_input(tvtk.ConeSource(center=(0,-0.05,0), direction=(0,1,0), height=0.3, radius=0.35, resolution=16).output)
joint_model_movable_x = collect.output
# ground in y direction
collect = tvtk.AppendPolyData()
collect.add_input(tvtk.CubeSource(x_length=0.2, y_length=0.7, z_length=0.25, center=(-0.5, 0, 0)).output)
collect.add_input(tvtk.SphereSource(radius=0.11, center=(-0.3, -0.18, 0), phi_resolution=24, theta_resolution=24).output)
collect.add_input(tvtk.SphereSource(radius=0.11, center=(-0.3, 0.18, 0), phi_resolution=24, theta_resolution=24).output)
collect.add_input(tvtk.ConeSource(center=(-0.05,0,0), direction=(1,0,0), height=0.3, radius=0.35, resolution=16).output)
joint_model_movable_y = collect.output


# ----------- Materials ------------
air = ElementMaterial(E = 100.0, density = 1.0, radius = 0.000001, maximum_stress = 1e+99, name='deleted')
steels = []
#radii = N.concatenate((N.arange(0.005, 0.05, 0.005), N.arange(0.05, 0.25, 0.01)))
radii = N.concatenate((N.arange(0.005, 0.05, 0.005), N.arange(0.05, 0.25, 0.01)))
for radius in radii:
    steels.append(ElementMaterial(E = 2.1e+11, density = 7700.0, maximum_stress = 10000000., radius = radius, name='steel'))

# rajan:
rajan_materials = []
# rajan defines cross sectional areas (in in^2)
cross_sectional_areas = N.array([1.62, 1.8, 2.38, 2.62, 2.88, 3.09, 3.13, 3.38, 3.63, 3.84, 
    4.18, 4.49, 4.8, 4.97, 5.12, 5.74, 7.22, 7.97, 
    11.5, 13.5, 13.9, 14.2, 15.5, 16.0, 28.8, 19.9, 
    22.0, 22.9, 26.5, 30.0, 33.5])
radii = N.sqrt(cross_sectional_areas*0.00064516 / N.pi)
for radius in radii:
    rajan_materials.append(ElementMaterial(E = 6.8947573e+10, density = 2767.9905, maximum_stress = 1.7236893e+08, radius=radius, name='aluminium'))

# ----------- Sample Constructions ------------
def bridge2():
    #   j6--j7--j8
    #  / |\/| \/| \
    # /  |/\| /\|  \
    #j1--j2-j3--j4--j5
    j1 = Joint(x=-6, y=0, movable_x=False, movable_y=False,
            mutatable_x=False, mutatable_y=False)
    j2 = Joint(x=-3, y=0, movable_x=True, movable_y=True,
            mutatable_x=True, mutatable_y=False, force_y=-1)
    j3 = Joint(x=0, y=0, movable_x=True, movable_y=True,
            mutatable_x=True, mutatable_y=False, force_y=-1)
    j4 = Joint(x=3, y=0, movable_x=True, movable_y=True,
            mutatable_x=True, mutatable_y=False, force_y=-1)
    j5 = Joint(x=6, y=0, movable_x=False, movable_y=False,
            mutatable_x=False, mutatable_y=False)
    j6 = Joint(x=-4, y=3, movable_x=True, movable_y=True,
            mutatable_x=True, mutatable_y=True)
    j7 = Joint(x=0, y=3, movable_x=True, movable_y=True,
            mutatable_x=True, mutatable_y=True)
    j8 = Joint(x=4, y=3, movable_x=True, movable_y=True,
            mutatable_x=True, mutatable_y=True)
    mat = steels[2]
    elems = [
        Element(joint1=j1, joint2=j6, material=mat),
        Element(joint1=j1, joint2=j2, material=mat),
        Element(joint1=j2, joint2=j6, material=mat),
        Element(joint1=j2, joint2=j3, material=mat),
        Element(joint1=j2, joint2=j7, material=mat),
        Element(joint1=j3, joint2=j6, material=mat),
        Element(joint1=j3, joint2=j4, material=mat),
        Element(joint1=j3, joint2=j7, material=mat),
        Element(joint1=j3, joint2=j8, material=mat),
        Element(joint1=j4, joint2=j7, material=mat),
        Element(joint1=j4, joint2=j8, material=mat),
        Element(joint1=j4, joint2=j5, material=mat),
        Element(joint1=j5, joint2=j8, material=mat),
        Element(joint1=j6, joint2=j7, material=mat),
        Element(joint1=j7, joint2=j8, material=mat),
        ]
    return Construction(joints=[j1,j2,j3,j4,j5,j6,j7,j8], elements=elems, force_magnitude=10000.)


def hang():
    j1 = Joint(x=-6, y=4, mutatable_x=False, mutatable_y=False)
    j2 = Joint(x=-2, y=4, mutatable_x=False, mutatable_y=False)
    j3 = Joint(x= 2, y=4, mutatable_x=False, mutatable_y=False)
    j4 = Joint(x= 6, y=4, mutatable_x=False, mutatable_y=False)
    j5 = Joint(x= 0, y=-4, mutatable_x=False, mutatable_y=False, movable_x=True, movable_y=True, force_y=1)
    mat = steels[4]
    elems = [
            Element(joint1=j1, joint2=j5, material=mat),
            Element(joint1=j2, joint2=j5, material=mat),
            Element(joint1=j3, joint2=j5, material=mat),
            Element(joint1=j4, joint2=j5, material=mat),
            ]
    return Construction(joints=[j1,j2,j3,j4,j5], elements=elems, available_element_materials=steels, element_deleted_material=air, force_magnitude=10000., width=14, height=10)
