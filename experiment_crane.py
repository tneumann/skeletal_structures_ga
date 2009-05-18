import pickle
import numpy as N
import ga_truss as GAT
import screenshot

# load base construction
con = pickle.load(open('examples/crane_onlytopo.con'))

# setup GA
ga = GAT.default_genetic_algorithm
ga.set(
        population_size = 100, 
        crossover_rate = 0.6, 
        mutation_rate = 0.1
    )
ga.init_individual.set(
        copy_probability = 0, 
        joint_offset = 8, 
        mutate_probability = 0.2
    )
ga.selection_operator.arena_size = 3
ga.survival.set(
        force_variety_shape = True, 
        force_variety_topology = False, 
        force_variety_size = False, 
        num_similar_individuals = 5, 
        shape_variety_distance = 2.0
    )
ga.world.deleted_material_gene_places = 2

for load in [10000., 50000., 75000., 100000., 200000.]:
    # run GA
    con.joint_force_magnitude = load
    ga.world.construction = con
    ga.init_evolution()
    for gen in xrange(500):
        ga.evolution_step()
    # save results
    best = ga.current_population.best
    best_con = ga.world.construction_from_individual(best)
    pickle.dump(best_con, open('/tmp/load%08d.con' % load, 'w'))
    screenshot.figure(best_con, path = ('/tmp/load%08d.eps' % load), size_factor=0.2, radius_amplify=25)
    print( "load %10.2f: %10.2fkg, %3.2f%% stability" % (load, best_con.simulation_result.construction_mass, best_con.simulation_result.stability*100.) )
