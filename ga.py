#! /usr/bin/env python

# Copyright 2009 Thomas Neumann
#
# Redistribution of this file is permitted under
# the terms of the GNU Public License (GPL) version 2.

import numpy as N
from enthought.traits.api import *
import random


class BaseIndividual(HasTraits):
    # saved in the base_fitness property is the value directly returned from the fitness function
    base_fitness = Float(0.0)

    # raw_fitness is the fitness after survival phase
    raw_fitness = Float(0.0)

    # the age property will be used when elitism is on and tells how many generations the individual survived
    age = Int(0)

    def clone(self):
        ''' clone the individual, this should be implemented by subclasses '''


class Population(HasTraits):
    number = Int()
    # all individuals, sorted by raw_fitness
    individuals = List()

    best = Property()
    worst = Property()

    def _get_best(self):
        return self.individuals[0]

    def _get_worst(self):
        return self.individuals[-1]

    @property
    def mean_fitness(self):
        return N.mean([i.raw_fitness for i in self.individuals])

    @property
    def mean_age(self):
        return N.mean([i.age for i in self.individuals])


def no_competition_survival(world, individuals):
    ''' each individual lives on its own without competition with other individuals '''
    for individual in individuals:
        individual.raw_fitness = individual.base_fitness


class FitnessProportionateSelection(HasTraits):
    def __init__(self, **kwargs):
        super(FitnessProportionateSelection, self).__init__(**kwargs)
        self._pop_in_cache = None

    def __call__(self, world, population):
        if self._pop_in_cache != population:
            fs = 0.0
            self._cumulative_individual_sums = []
            for i in population.individuals:
                fs += 1.0 / i.raw_fitness
                self._cumulative_individual_sums.append(fs)
            self._fitness_sum = fs
            self._pop_in_cache = population
        # select a random number between 0 and the sum of the raw fitness of the populations individuals
        r = randrange(0.0, self._fitness_sum, int=float)
        # find the individual by linear search
        cumfs = self._cumulative_individual_sums
        fs = cumfs[0]
        i = 0
        while fs < r:
            i += 1
            fs = cumfs[i]
        return population.individuals[i]

class TournamentSelection(HasTraits):
    arena_size = Int(5)

    def __call__(self, world, population):
        contestants = [population.individuals[random.randint(0, len(population.individuals)-1)] for a in xrange(self.arena_size)]
        return min(contestants, key = lambda individual: individual.raw_fitness)


class GeneticAlgorithm(HasTraits):
    # the world object can be any object which holds data that helps the genetic operators and fitness functions
    # it can contain further settings for the genetic algorithm or utility functions to evaluate individuals
    world = Any()

    # initialization function, takes as input the given "world", returns a individual, e.g. with random genes
    init_individual = Callable()

    # mutation operator takes the mutation_rate, world, and an individual and mutates the individual in-place
    mutation_operator = Callable()

    # crossover operator, takes crossover_rate, world, and 4 individuals (mummy, daddy, sister, brother). 
    # it should modify sister and brother via crossover operations in-place
    crossover_operator = Callable()

    # the selection operator is a function that takes the world and a population and returns 
    # one individual that should be selected in the reproduction phase
    selection_operator = Callable()

    # fitness function takes the world and an individual and returns it's fitness
    fitness_function = Callable()

    # the survival operator takes a world and a whole population of one generation and should set the raw_fitness properties 
    # of all the individuals based upon the base_fitness of the individual
    # in this function one can implement competition strategies or steer the evolution, e.g. by giving
    # duplicate individuals a fitness penalty
    survival = Callable(no_competition_survival)

    # the validation function can be empty, it takes the world and an individual and returns True 
    # if the individual is valid and should go into the population
    validator = Callable()

    # crossover and mutation probabilities
    crossover_rate = Range(0.0, 1.0, 0.7)
    mutation_rate = Range(0.0, 1.0, 0.1)

    # step event executed whenever a fresh generation has been evaluated
    on_step = Event()

    populations = List(Population)
    population_size = Int(20)
    current_population = Property(Instance(Population), depends_on='populations')
    num_steps = Int(0)
    num_keep_generations = Int(-1) # how many populations to keep. -1 == keep all
    current_best_raw_fitness = Property(Int, depends_on='populations')
    current_worst_raw_fitness = Property(Int, depends_on='populations')
    current_mean_raw_fitness = Property(Int, depends_on='populations')

    # if elitism is turned on, then children will compete with their parents, and good individuals will
    # not die but stay for as long as they are competed by better individuals
    # when elitism is off parents will be killed after reproduction
    elitism = Bool(True)

    inited = Bool(False)

    def _get_current_population(self):
        if len(self.populations) > 0:
            return self.populations[-1]
        else:
            return None

    def _get_current_best_raw_fitness(self):
        return self.current_population.best.raw_fitness

    def _get_current_worst_raw_fitness(self):
        return self.current_population.worst.raw_fitness

    def _get_current_mean_raw_fitness(self):
        return self.current_population.mean_fitness

    def init_evolution(self):
        ''' initialize empty population from initializer function and evaluate first population '''
        cur_pop = Population(individuals = [self.init_individual(self.world) for i in xrange(self.population_size)], number=1)
        self.populations = [cur_pop]
        # calculate base fitness for each of the individuals
        for individual in cur_pop.individuals:
            self._call_birth_callback(individual)
            individual.base_fitness = self.fitness_function(self.world, individual)
        self._calc_base_fitness(cur_pop)
        self.survival(self.world, cur_pop.individuals)
        self.num_steps = 0
        self.inited = True

    def evolution_step(self):
        ''' complete genetic algorithm cycle including selection and reproduction '''
        # this is the core of the genetic algorithm
        # this algorithm tries to minimize the fitness function
        cur_pop = self.current_population
        # make a new empty population
        new_pop = Population(number=cur_pop.number+1)
        # make love until new population is full, always get 2 children
        while len(new_pop.individuals) < self.population_size:
            # select 2 individuals for reproduction
            mummy = self.selection_operator(self.world, cur_pop)
            daddy = self.selection_operator(self.world, cur_pop)
            # make 2 children, which are clones of mummy and daddy first, change genes through crossover/mutation later
            sister = mummy.clone()
            brother = daddy.clone()
            # apply crossover operator
            self.crossover_operator(self.crossover_rate, self.world, mummy, daddy, sister, brother)
            # apply mutation operator
            self.mutation_operator(self.mutation_rate, self.world, sister)
            self.mutation_operator(self.mutation_rate, self.world, brother)
            # call birth callback for individuals if available
            self._call_birth_callback(sister)
            self._call_birth_callback(brother)
            # put children into new population if they are valid
            if self.validator is None or self.validator(self.world, sister):
                new_pop.individuals.append(sister)
            if self.validator is None or self.validator(self.world, brother):
                new_pop.individuals.append(brother)
        # if population size is odd then we may have made one child too much, so just randomly kill one
        if len(new_pop.individuals) > len(cur_pop.individuals):
            del new_pop.individuals[rand_list_index(len(cur_pop.individuals))]
        # calculate base fitness for each new child
        self._calc_base_fitness(new_pop)
        # now put back the parents into the new population when elitism is enabled
        if self.elitism:
            new_pop.individuals += cur_pop.individuals
            for individual in cur_pop.individuals:
                individual.age += 1
        self.survival(self.world, new_pop.individuals)
        # sort individuals by fitness
        new_pop.individuals.sort(key=lambda i: i.raw_fitness)
        # kill individuals so that population size stays the same
        new_pop.individuals = new_pop.individuals[:self.population_size]
        # done with this generation :-)
        self.populations.append(new_pop)
        if self.num_keep_generations != -1 and len(self.populations) > self.num_keep_generations:
            self.populations = self.populations[-self.num_keep_generations:]
        self.num_steps += 1
        self.on_step = new_pop

    def _calc_base_fitness(self, population):
        for individual in population.individuals:
            individual.base_fitness = self.fitness_function(self.world, individual)

    def _call_birth_callback(self, individual):
        if hasattr(individual, 'on_birth'):
            individual.on_birth(self.world)

# --- 
# some utility functions that can be used by custom-made mutation/crossover operators on numpy-array based genes
def numpy_array_2point_crossover(mummy, daddy):
    ''' do a 2 point crossover on 2 1-dimensional numpy arrays, return the 2 children arrays '''
    # select crossover points
    index1_ = rand_list_index(len(mummy))
    index2_ = rand_list_index(len(mummy))
    index1 = min(index1_, index2_)
    index2 = max(index1_, index2_)
    # do crossover
    sister  = N.concatenate((mummy[:index1], daddy[index1:index2], mummy[index2:]))
    brother = N.concatenate((daddy[:index1], mummy[index1:index2], daddy[index2:]))
    return sister, brother

def numpy_array_uniform_integer_crossover(mummy, daddy, crossover_probability):
    do_crossover = N.random.random(len(mummy)) < crossover_probability
    r = N.random.random(len(mummy))
    selector = do_crossover & (r<=0.5)
    sister  = N.where( selector, mummy, daddy)
    brother = N.where( selector, daddy, mummy)
    return sister, brother

def numpy_array_uniform_real_crossover(mummy, daddy, crossover_probability):
    do_crossover = N.random.random(len(mummy)) < crossover_probability
    r = N.random.random(len(mummy))
    sister  = N.where( do_crossover, mummy, r * mummy + (1-r) * daddy)
    brother = N.where( do_crossover, daddy, r * daddy + (1-r) * mummy)
    return sister, brother

def numpy_array_uniform_integer_mutation(gene, mutation_probability, possible_gene_values):
    mutant = possible_gene_values[ N.random.randint(len(possible_gene_values), size=len(gene)) ]
    do_mutate = N.random.random(len(gene)) <= mutation_probability
    return N.where( do_mutate, mutant, gene )

def numpy_array_uniform_real_offset_mutation(gene, mutation_probability, maximum_offset):
    offset = (N.random.random(len(gene)) * 2.0 - 0.5) * maximum_offset
    do_mutate = N.random.random(len(gene)) <= mutation_probability
    return N.where( do_mutate, gene + offset, gene )

def rand_list_index(list_len):
    return random.randint(0, list_len-1)


if __name__ == '__main__':
    # a sample genetic algorithm:
    # solve the following problem via a genetic algorithm:
    # the gene of the individual is an array of N elements with numbers from 0 to N-1
    # the target is to develop an individual that has the numbers 0..N-1 ordered in its gene
    n = 10

    # define individual
    class TestIndividual(BaseIndividual):
        gene = CArray(dtype=N.uint)
        def clone(self):
            return self.clone_traits(trait_names='gene')

    # fitness: distance to optimal solution
    def fitness(world, individual):
        return N.sum(N.abs(individual.gene - N.r_[:n]))

    # initialize individuals with random gene
    def init(world):
        return TestIndividual( gene = N.random.randint(n, size=n).astype(N.uint) )

    # define crossover
    def crossover(crossover_rate, world, mummy, daddy, brother, sister):
        sister.gene, brother.gene = numpy_array_2point_crossover(mummy.gene, daddy.gene)

    # define mutation
    def mutation(mutation_rate, world, individual):
        individual.gene = numpy_array_uniform_integer_mutation(individual.gene, mutation_rate, N.r_[:n])

    ga = GeneticAlgorithm(
            init_individual = init,
            selection_operator = FitnessProportionateSelection(),
            crossover_operator = crossover,
            mutation_operator = mutation,
            fitness_function = fitness,
            population_size = n*2, elitism = True)

    ga.init_evolution()
    while ga.current_population.best.raw_fitness > 0:
        print "generation %3d: best: %3d , worst: %3d , mean: %3d , best gene: %s, mean age: %.3f" % (ga.num_steps,\
                ga.current_population.best.raw_fitness, ga.current_population.worst.raw_fitness, ga.current_population.mean_fitness,
                str(ga.current_population.best.gene), ga.current_population.mean_age)
        ga.evolution_step()
    print "---------------------------------------------------------------"
    print "solution found after %d generations: %s" % (ga.num_steps, str(ga.current_population.best.gene))

