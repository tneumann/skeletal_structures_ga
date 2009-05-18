#! /usr/bin/env python

# Copyright 2009 Thomas Neumann
#
# Redistribution of this file is permitted under
# the terms of the GNU Public License (GPL) version 2.

from enthought.traits.api import *
import numpy as N

import ga as GA
import physics
import model
import data


class World(HasTraits):
    # the basis construction, it should not be changed anymore when set and during the genetic algorithm
    construction = Instance(model.Construction)

    deleted_material_gene_places = Int(1)
    gene_index_for_element_material = Property(Instance(Dict), depends_on='construction.available_element_materials')

    joint_position_gene_is_x = Property(CArray(dtype=N.bool), depends_on='construction')
    joint_position_gene_is_y = Property(CArray(dtype=N.bool), depends_on='construction')

    joint_position_gene_size = Property(Int, depends_on='construction')
    element_material_gene_size = Property(Int, depends_on='construction')

    @cached_property
    def _get_joint_position_gene_size(self):
        c = 0
        for joint in self.construction.joints:
            if joint.mutatable_x:
                c += 1
            if joint.mutatable_y:
                c += 1
        return c

    @cached_property
    def _get_element_material_gene_size(self):
        return len( N.flatnonzero(N.array([e.mutatable_material for e in self.construction.elements]).ravel() ) )

    @cached_property
    def _get_joint_position_gene_is_x(self):
        is_x = []
        for joint in self.construction.joints:
            if joint.mutatable_x:
                is_x.append(True)
            if joint.mutatable_y:
                is_x.append(False)
        return N.array(is_x)

    @cached_property
    def _get_joint_position_gene_is_y(self):
        return N.logical_not(self.joint_position_gene_is_x)

    @cached_property
    def _get_gene_index_for_element_material(self):
        #print "getit"
        d = {}
        for i,mat in enumerate(self.construction.available_element_materials):
            d[mat] = i
        d[self.construction.element_deleted_material] = -1
        return d

    def construction_from_individual(self, individual, construction=None):
        # TODO: check if construction compatible to base construction
        if construction is None:
            old_construction = self.construction
            construction = self.construction.clone_traits()
            construction.element_deleted_material = old_construction.element_deleted_material
        for joint,(x,y) in zip(construction.joints, self.joint_positions_for_individual(individual)):
            joint.x = x
            joint.y = y
        for element, mat in zip(construction.elements, self.materials_for_individual(individual)):
            element.material = mat
        construction.from_individual = individual
        return construction

    def individual_from_construction(self, construction=None, individual=None):
        # TODO: check if construction compatible to base construction
        if construction is None:
            construction = self.construction
        if individual is None:
            individual = ConstructionIndividual(joint_positions_gene = N.zeros(self.joint_position_gene_size, N.double), 
                element_material_gene = N.zeros(self.element_material_gene_size, N.uint))
        i = 0
        for joint in construction.joints:
            if joint.mutatable_x:
                individual.joint_positions_gene[i] = joint.x
                i += 1
            if joint.mutatable_y:
                individual.joint_positions_gene[i] = joint.y
                i += 1
        i = 0
        for element in construction.elements:
            if element.mutatable_material:
                individual.element_material_gene[i] = self.gene_index_for_element_material[element.material]
                i += 1
        return individual


    joint_position_gene_locations = Property(CArray(dtype=N.uint), depends_on='construction.joints.mutatable_+')

    @cached_property
    def _get_joint_position_gene_locations(self):
        i = 0
        gene_locations = []
        for joint in self.construction.joints:
            if joint.mutatable_x:
                gene_locations.append(i)
            i += 1
            if joint.mutatable_y:
                gene_locations.append(i)
            i += 1
        return N.array(gene_locations, N.uint)

    def joint_positions_for_individual(self, individual):
        jp = self.construction.joint_positions.copy().ravel()
        if len(individual.joint_positions_gene) > 0:
            jp[self.joint_position_gene_locations] = individual.joint_positions_gene
        return jp.reshape((-1,2))

    element_material_gene_locations = Property(CArray(dtype=N.uint), depends_on='construction.elements.[deletable,mutatable_material]')

    @cached_property
    def _get_element_material_gene_locations(self):
        return N.array([i for i,element in enumerate(self.construction.elements) if element.mutatable_material], N.uint)

    def materials_for_individual(self, individual):
        mats = [e.material for e in self.construction.elements]
        for i,loc in enumerate(self.element_material_gene_locations):
            g = individual.element_material_gene[i]
            if g < 0:
                mats[loc] = self.construction.element_deleted_material
            else:
                mats[loc] = self.construction.available_element_materials[g]
        return mats



class ConstructionIndividual(GA.BaseIndividual):
    joint_positions_gene = CArray(dtype=N.double)
    element_material_gene = CArray(dtype=N.int)

    def on_birth(self, world):
        # genotype -> phenotype ;-)
        #self.construction = world.construction_from_individual(self)
        self.joint_positions = world.joint_positions_for_individual(self)
        self.element_materials = world.materials_for_individual(self)
        # TODO: can be faster via lookup for next 4 calls
        element_E = N.array( [material.E for material in self.element_materials], N.double)
        element_r = N.array( [material.radius for material in self.element_materials], N.double)
        element_p = N.array( [material.density for material in self.element_materials], N.double)
        jfreedom = world.construction.get_adjusted_joint_freedom_array(self.element_materials)
        # analyze the individual
        self.joint_displacements, self.element_strains, self.element_stresses, self.mass, self.simulation_status, \
            (time_preparation, time_solution, time_postprocessing) = \
                    physics.analyse_truss(self.joint_positions, jfreedom, world.construction.loads_array, 
                            world.construction.element_index_table, element_E, element_r, element_p)
        # calculate constraints
        self.stability_constraint = N.max( N.abs(self.element_stresses) / world.construction.max_element_stress_array ) - 1
        # calculate displacement constraints
        self.displacement_constraint = N.max( N.abs(self.joint_displacements) / world.construction.max_joint_displacement_array ) - 1
        # combine constraints into feasibility
        self.feasability = max(self.stability_constraint, self.displacement_constraint)
        self.feasible = self.stability_constraint <= 0 and self.displacement_constraint <= 0

    def clone(self):
        c = ConstructionIndividual()
        c.joint_positions_gene = self.joint_positions_gene.copy()
        c.element_material_gene = self.element_material_gene.copy()
        return c


class CopyConstructionInitializer(HasTraits):
    def __call__(self, world):
        return world.individual_from_construction()

class RandomInitializer(HasTraits):
    joint_offset = Float(1.0)
    copy_probability = Float(0.1)
    mutate_probability = Float(0.5)
    validator = Callable()
    def __call__(self, world):
        if N.random.random() < self.copy_probability:
            individual = CopyConstructionInitializer()(world)
        else:
            while True:
                individual = world.individual_from_construction()
                individual.joint_positions_gene = JointPositionOffsetMutation(maximum_offset=self.joint_offset)(self.mutate_probability, world, individual.joint_positions_gene)
                individual.element_material_gene = uniform_random_material_selection_mutation(self.mutate_probability, world, individual.element_material_gene)
                individual.on_birth(world)
                if (individual.feasability <= 10 or world.construction.simulation_result.stability < 0) and self.validator(world, individual):
                    break
        return individual


class ValidateConstruction(HasTraits):
    check_joint_overlap = Bool(True)
    check_simulatable = Bool(True)
    kill_stability_constraint = Bool(False)
    kill_displacement_constraint = Bool(False)
    joint_overlap_tolerance = Float(0.5)

    def __call__(self, world, individual):
        valid = True
        # check if all elements are deleted
        valid &= N.any(individual.element_material_gene >= 0)
        # check if joint positions are out of construction
        jp = individual.joint_positions
        valid &= N.all(jp[:,0] <  (world.construction.width/2.)) and \
                 N.all(jp[:,0] > -(world.construction.width/2.)) and \
                 N.all(jp[:,1] <  (world.construction.height/2.)) and \
                 N.all(jp[:,1] > -(world.construction.height/2.))
        # check if any element has been deleted but is not allowed to
        valid &= N.all((individual.element_material_gene >= 0) | world.construction.element_deletable_array[world.element_material_gene_locations])
        # check if any joints are at the same location
        if self.check_simulatable:
            valid &= individual.simulation_status == 'success'
        if valid and self.check_joint_overlap:
            # build matrix of manhatten distances of each joint to each other joint
            d = N.sum(N.abs(jp[:,N.newaxis,:] - jp[N.newaxis,:,:]), axis=2)
            # there should be no distance entry larger then the tolerance, only the diagonal should be zero distance 
            # (because there the same joints are compared)
            valid &= ((d <= self.joint_overlap_tolerance) == N.eye(len(jp), dtype=N.bool)).all()
        # check if all the elements have the "deleted element" material
            #valid &= any((material != self.element_deleted_material for material in individual.element_materials))
        if valid and self.kill_stability_constraint:
            valid &= individual.stability_constraint > 0 
        if valid and self.kill_displacement_constraint:
            valid &= individual.displacement_constraint > 0
        return valid


class MassOfConstructionFitness(HasTraits):
    name = String('construction mass')

    def __call__(self, world, individual):
        return individual.mass


class ConstructionSurvival(HasTraits):
    num_similar_individuals = Int(5)
    force_variety_topology = Bool(True)
    force_variety_size = Bool(False)
    force_variety_shape = Bool(True)
    shape_variety_distance = Float(0.5)
    topology_variety_distance = Int(1)

    unfeasible_constraint_penalty = Bool(True)

    def __call__(self, world, individuals):
        # find out worst fitness from all individuals that are feasible (no constraint violated)
        base_fitnesses_feasible = [i.base_fitness for i in individuals if i.feasible]
        if len(base_fitnesses_feasible) > 0:
            worst_feasible_fitness = max(base_fitnesses_feasible)
            for individual in individuals:
                individual.raw_fitness = individual.base_fitness
                if self.unfeasible_constraint_penalty and not individual.feasible:
                    #print "penalty ", individual.feasability
                    individual.raw_fitness = individual.raw_fitness + worst_feasible_fitness * (individual.feasability+1)
        else:
            #print "all unfeasible, choosing stability fitness function!"
            worst_fitness = max([i.base_fitness for i in individuals])
            for individual in individuals:
                individual.raw_fitness = (individual.feasability+1) + worst_fitness

        # TODO: document/cleanup/paremetrize
        # variety reward
        if self.force_variety_topology and len(individuals[0].element_material_gene) > 0:
            # topology (1=member present, 0=member absent)
            all_element_genes = N.array([i.element_material_gene < 0 for i in individuals], N.int)
            all_fitness = N.array([i.raw_fitness for i in individuals])
            bad_individual = N.zeros(len(individuals), N.bool)
            for i,individual in enumerate(individuals):
                # find duplicates to that topology
                duplicates = N.flatnonzero(N.abs(all_element_genes[i] - all_element_genes).sum(axis=1) < self.topology_variety_distance) # indices to duplicates
                # select worst (highest fitness) individuals from those duplicaates
                if len(duplicates) > self.num_similar_individuals:
                    worst_duplicates = duplicates[all_fitness[duplicates].argsort()[self.num_similar_individuals:]]
                    bad_individual[worst_duplicates] = True
            for individual_index in N.flatnonzero(bad_individual):
                individual = individuals[individual_index]
                individual.raw_fitness = individual.raw_fitness + 1e+6

        if self.force_variety_size and len(individuals[0].element_material_gene) > 0:
            all_element_genes = N.maximum(N.array([i.element_material_gene for i in individuals]), -1) # max(..,-1) because deleted elements are all < 0
            all_fitness = N.array([i.raw_fitness for i in individuals])
            bad_individual = N.zeros(len(individuals), N.bool)
            #print all_element_genes
            for i,individual in enumerate(individuals):
                # find duplicates to that topology
                duplicates = N.flatnonzero((all_element_genes[i] == all_element_genes).all(axis=1)) # indices to duplicates
                # select worst (highest fitness) individuals from those duplicaates
                if len(duplicates) > self.num_similar_individuals:
                    worst_duplicates = duplicates[all_fitness[duplicates].argsort()[self.num_similar_individuals:]]
                    bad_individual[worst_duplicates] = True
            for individual_index in N.flatnonzero(bad_individual):
                individual = individuals[individual_index]
                individual.raw_fitness = individual.raw_fitness + 1e+6
                #print "bad: ", individual.element_material_gene

        if self.force_variety_shape and len(individuals[0].joint_positions_gene) > 0:
            all_jp_genes = N.array([i.joint_positions_gene for i in individuals])
            all_fitness = N.array([i.raw_fitness for i in individuals])
            bad_individual = N.zeros(len(individuals), N.bool)
            for i,individual in enumerate(individuals):
                duplicates = N.flatnonzero((N.abs(individual.joint_positions_gene - all_jp_genes).sum(axis=1) < self.shape_variety_distance)) # indices to duplicates
                # select worst (highest fitness) individuals from those duplicaates
                if len(duplicates) > self.num_similar_individuals:
                    worst_duplicates = duplicates[all_fitness[duplicates].argsort()[self.num_similar_individuals:]]
                    bad_individual[worst_duplicates] = True
            for individual_index in N.flatnonzero(bad_individual):
                individual = individuals[individual_index]
                individual.raw_fitness = individual.raw_fitness + 1e+6
                #print "bad: ", individual.element_material_gene


def uniform_joint_position_crossover(crossover_rate, world, mummy, daddy):
    return GA.numpy_array_uniform_real_crossover(mummy, daddy, crossover_rate)

def uniform_element_material_crossover(crossover_rate, world, mummy, daddy):
    return GA.numpy_array_uniform_integer_crossover(mummy, daddy, crossover_rate)

def twopoint_element_material_crossover(crossover_rate, world, mummy, daddy):
    if crossover_rate < N.random.random():
        return GA.numpy_array_2point_crossover(mummy, daddy)
    else:
        return mummy.copy(), daddy.copy()

# we have individuals with different genes, so we need a composite crossover operator
class CompositeCrossover(HasTraits):
    joint_position_crossover = Callable()
    element_material_crossover = Callable()

    def __call__(self, crossover_rate, world, mummy, daddy, sister, brother):
        # apply crossover to joint positions gene
        sister.joint_positions_gene, brother.joint_positions_gene = \
            self.joint_position_crossover(crossover_rate, world, 
                    mummy.joint_positions_gene, daddy.joint_positions_gene)
        # apply crossover to element material gene
        sister.element_material_gene, brother.element_material_gene = \
            self.element_material_crossover(crossover_rate, world, 
                    mummy.element_material_gene, daddy.element_material_gene)


def uniform_random_material_selection_mutation(mutation_rate, world, gene):
    return GA.numpy_array_uniform_integer_mutation(gene, mutation_rate, 
            N.r_[ -world.deleted_material_gene_places : len(world.construction.available_element_materials) ])

class JointPositionOffsetMutation(HasTraits):
    maximum_offset = Float(10.0)
    def __call__(self, mutation_rate, world, gene):
        return GA.numpy_array_uniform_real_offset_mutation(gene, mutation_rate, self.maximum_offset)

class CompositeMutation(HasTraits):
    joint_position_mutation = Callable()
    element_material_mutation = Callable()

    def __call__(self, mutation_rate, world, individual):
        # apply mutation to joint positions gene
        individual.joint_positions_gene = self.joint_position_mutation(mutation_rate, world, individual.joint_positions_gene)
        # apply mutation to element material gene
        individual.element_material_gene = self.element_material_mutation(mutation_rate, world, individual.element_material_gene)


_validator = ValidateConstruction()
default_genetic_algorithm = GA.GeneticAlgorithm(
        fitness_function = MassOfConstructionFitness(),
        validator = _validator,
        survival = ConstructionSurvival(),
        init_individual = RandomInitializer(validator=_validator),
        selection_operator = GA.TournamentSelection(arena_size=3),
        crossover_operator = CompositeCrossover(
            joint_position_crossover = uniform_joint_position_crossover,
            element_material_crossover = uniform_element_material_crossover ),
            #element_material_crossover = twopoint_element_material_crossover ),
        mutation_operator = CompositeMutation(
            joint_position_mutation = JointPositionOffsetMutation(maximum_offset=5.0),
            element_material_mutation = uniform_random_material_selection_mutation ),
        world = World(
            deleted_material_gene_places = len(data.steels),
        ) )

