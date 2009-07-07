# this script demonstrates the capabilities of the truss ga modules to be used in non-gui automated experiment scripts
# it executes a genetic algorithm for a construction with multiple loads
# the construction was built with the editor, and the GA settings were set in the GUI
# the construction was then saved to examples/crane_onlytopo.con
# this script loads the construction and the GA, and then changes
# the load applied to the construction, exectues the GA, 
# and saves the resulting constructions after 500 generations as an EPS graphic
import pickle
import screenshot

# load base construction and ga
con, ga = pickle.load(open('examples/crane_onlytopo.con'))

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

