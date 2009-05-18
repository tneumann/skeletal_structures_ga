#! /usr/bin/env python

# Copyright 2009 Thomas Neumann
#
# Redistribution of this file is permitted under
# the terms of the GNU Public License (GPL) version 2.

import pdb
from enthought.traits.api import *
import numpy as N
from time import time
import physics


class Joint(HasTraits):
    x = Float()
    y = Float()
    # TODO: max_x, min_x, max_y, min_y?
    movable_x = Bool(False)
    movable_y = Bool(False)

    # represent force only as direction, can set the magnitude via the attribute "joint_force_magnitude" in Construction
    force_x = Range(low=-1., high=1., value=0.)
    force_y = Range(low=-1., high=1., value=0.)

    mutatable_x = Bool(True)
    mutatable_y = Bool(True)

    max_displacement_x = Float(0.05)
    max_displacement_y = Float(0.05)

    @property
    def position(self):
        return (self.x, self.y)

    @property
    def force(self):
        return (self.force_x, self.force_y)

    @property
    def movable_xy(self):
        return self.movable_x and self.movable_y

    @on_trait_change('movable_+, force_+', post_init=True)
    def _kill_force_when_not_movable(self):
        if not self.movable_xy:
            self.force_x, self.force_y = 0, 0


class ElementMaterial(HasTraits):
    E = Float() # elasticity module in N/m^2 (Pa)
    density = Float() # kg/m^3
    radius = Float() # meter
    maximum_stress = Float() # N/m^2
    name = String()

    @property
    def maximum_stress_mpa(self):
        return self.maximum_stress * 1e-06


class Element(HasTraits):
    joint1 = Instance(Joint)
    joint2 = Instance(Joint)
    material = Instance(ElementMaterial)
    mutatable_material = Bool(True)
    deletable = Bool(True)


class SimulationResult(HasTraits):
    joint_displacements = CArray(dtype=N.double, shape=(None,2))
    element_strains = CArray(dtype=N.double)
    element_stresses = CArray(dtype=N.double)
    stability = Float()
    construction_mass = Float()

    status = String('not simulated', transient=True)
    time_preparation = Float()
    time_solution = Float()
    time_postprocessing = Float()
    time_conversion = Float()

    @property
    def okay(self):
        return self.status == 'success'

    @property
    def time(self):
        return self.time_preparation + self.time_solution + self.time_postprocessing + self.time_conversion

    @property
    def time_preparation_percent(self):
        return (self.time_preparation / self.time) * 100.

    @property
    def time_solution_percent(self):
        return (self.time_solution / self.time) * 100.

    @property
    def time_postprocessing_percent(self):
        return (self.time_postprocessing / self.time) * 100.

    @property
    def time_conversion_percent(self):
        return (self.time_conversion / self.time) * 100.

    #@property
    #def stability(self):
    #    return max(100 - (self.element_breakage.max() * 100.), 0)


class Construction(HasTraits):
    joints = List(Joint)
    elements = List(Element)
    width = Float(20.)
    height = Float(15.)
    weight_factor = Float(1.0)
    joint_force_magnitude = Float(10000.) # scaling of force vectors
    available_element_materials = List(Instance(ElementMaterial, ()))
    element_deleted_material = Instance(ElementMaterial)

    def build_full_connections(self, material):
        if material is None:
            raise ValueError, 'given material cannot be None'
        connectivity = {}
        for j in self.joints:
            connectivity[j] = []
        for e in self.elements:
            connectivity[e.joint1].append(e.joint2)
            connectivity[e.joint2].append(e.joint1)
        for j1 in self.joints:
            for j2 in self.joints:
                if not j2 in connectivity[j1] and j1 != j2:
                    self.elements.append( Element(joint1=j1, joint2=j2, material=material) )
                    connectivity[j1].append(j2)
                    connectivity[j2].append(j1)


    joint_positions = Property(CArray(dtype=N.double, shape=(None,2)), transient=True, depends_on='joints.[x,y]')

    @cached_property
    def _get_joint_positions(self):
        return N.array( [(j.x, j.y) for j in self.joints], N.double)

    element_index_table = Property(CArray(dtype=N.uint, shape=(None,2)), transient=True, depends_on='elements.joint+')

    @cached_property
    def _get_element_index_table(self):
        if len(self.elements) == 0 or len(self.joints) == 0:
            return N.zeros((0,2), N.uint)
        # build index table for joint connections
        index_for_joint = {}
        for i,joint in enumerate(self.joints):
            index_for_joint[joint] = i
        return N.array( [(index_for_joint[e.joint1], index_for_joint[e.joint2]) for e in self.elements], N.uint)
        

    joint_freedom_array = Property(CArray(dtype=N.bool, shape=(None,2)), depends_on='joints.movable_+')
    @cached_property
    def _get_joint_freedom_array(self):
        return N.array( [(j.movable_x, j.movable_y) for j in self.joints], N.bool)

    loads_array = Property(CArray(dtype=N.double, shape=(None,2)), depends_on='joint_force_magnitude, joints.force_+')
    @cached_property
    def _get_loads_array(self):
        return N.array( [(j.force_x, j.force_y) for j in self.joints], N.double) * self.joint_force_magnitude

    max_element_stress_array = Property(CArray(dtype=N.bool, shape=(None,2)), depends_on='elements.material.maximum_stress')
    @cached_property
    def _get_max_element_stress_array(self):
        return N.array( [e.material.maximum_stress for e in self.elements], N.double )

    max_joint_displacement_array = Property(CArray(dtype=N.bool, shape=(None,2)), depends_on='joints.max_displacement_+')
    @cached_property
    def _get_max_joint_displacement_array(self):
        return N.array( [(j.max_displacement_x, j.max_displacement_y) for j in self.joints], N.double )

    element_deletable_array = Property(CArray(dtype=N.bool), depends_on='elements.deletable')
    @cached_property
    def _get_element_deletable_array(self):
        return N.array( [e.deletable for e in self.elements] )


    simulation_result = Property(Instance(SimulationResult), transient=True,
            depends_on='joints, elements, joints.[x,y,movable+,force_+], elements.joint+, elements.material.+, joint_force_magnitude, weight_factor')

    @property
    def simulated(self):
        return self.simulation_result != None

    @cached_property
    def _get_simulation_result(self):
        result = SimulationResult()
        if len(self.elements) > 0 and len(self.joints) > 0:
            #print "simulate"
            before = time()
            # convert construction to appropriate numpy arrays for physical simulation
            jpos = self.joint_positions
            # build element arrays
            elements = self.element_index_table
            element_E = N.array( [e.material.E for e in self.elements], N.double)
            element_r = N.array( [e.material.radius for e in self.elements], N.double)
            element_p = N.array( [e.material.density for e in self.elements], N.double)
            #max_stresses = N.array( [e.material.maximum_stress for e in self.elements], N.double )
            tconv = time() - before
            jfreedom = self.get_adjusted_joint_freedom_array([e.material for e in self.elements])
            # analyse construction
            result.joint_displacements, result.element_strains, result.element_stresses, result.construction_mass, result.status, \
                (result.time_preparation, result.time_solution, result.time_postprocessing) = \
                        physics.analyse_truss(jpos, jfreedom, self.loads_array, elements, element_E, element_r, element_p, 
                                weight_factor=self.weight_factor)
            #result.element_breakage = N.abs(result.element_stresses) / N.array([e.material.maximum_stress for e in self.elements], N.double)
            result.stability = 1 - (N.max(N.abs(result.element_stresses) / self.max_element_stress_array))
        else:
            result.time_preparation, result.time_solution, result.time_postprocessing = 0.000001, 0.000001, 0.000001
        return result

    def _joints_changed(self, old, new):
        self._handle_joints(old, new)

    def _joints_items_changed(self, listevent):
        self._handle_joints(listevent.removed, listevent.added)

    def _handle_joints(self, removed, added):
        for r in removed:
            # remove elements that were connected to this joint
            elements_to_remove = []
            for e in self.elements:
                if e.joint1 == r or e.joint2 == r:
                    elements_to_remove.append(e)
            for e in elements_to_remove:
                self.elements.remove(e)

    def _elements_changed(self, old, new):
        self._handle_elements(old, new)

    def _elements_items_changed(self, listevent):
        self._handle_elements(listevent.removed, listevent.added)

    def _is_duplicate_element(self, a):
        for e in self.elements:
            if e != a and (e.joint1 == a.joint1 and e.joint2 == a.joint2) or (e.joint2 == a.joint1 and e.joint1 == a.joint2):
                return True
        return False

    def element_addable(self, a):
        return not self._is_duplicate_element(a) and self._element_material_available(a.material)

    def _element_material_available(self, material):
        return material in self.available_element_materials or material == self.element_deleted_material

    def _handle_elements(self, removed, added):
        # check if addable
        for a in added:
            if self._is_duplicate_element(a):
                raise DuplicateElement
            #self._check_element_material_change(a, a.material)
            #a.on_trait_change(self._check_element_material_change, 'material')

    #def _check_element_material_change(self, element, new_material):
    #    if not (new_material in self.available_element_materials or new_material == self.element_deleted_material):
    #        print "auto-add material: %f, %s", (new_material.radius, new_material.name)
    #        self.available_element_materials.append(new_material)

    def get_adjusted_joint_freedom_array(self, new_element_materials):
        # pin joint if only deleted members are connected to the joint in question
        # to do this check for each joint if it has at least on non-deleted element
        has_non_deleted_elements = N.zeros(len(self.joints), N.bool)
        for (j1,j2), mat in zip(self.element_index_table, new_element_materials):
            if mat != self.element_deleted_material:
                has_non_deleted_elements[j1] = True
                has_non_deleted_elements[j2] = True
        # joint is discarded as beeing "deleted" if it has no non-deleted element ( -> only deleted elements connected)
        jfreedom = self.joint_freedom_array.copy()
        jfreedom[~has_non_deleted_elements] = (False, False)
        return jfreedom


class DuplicateElement(Exception):
    pass

class ElementMaterialNotAvailable(Exception):
    pass
