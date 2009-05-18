from IPython.Shell import IPShellEmbed
import IPython
import IPython.Shell
import pdb
import numpy as N
import time

def elem_lengths(jpos, elements):
    edge = jpos[elements[:,1]] - jpos[elements[:,0]]
    return N.sqrt(edge[:,0]**2 + edge[:,1]**2)

def elem_masses(ls, r, p):
    ''' return element mass in kg,
    given the element lengths (ls),
    element radius (r) in meter,
    and element material density in kg/m^3 '''
    return p * N.pi * r**2 * ls

gravity=9.798

def elem_stiffness(ls, E, r):
    ''' return element stiffness k,
    given element lengths (ls) in meter, 
    element eleasticity module in N/m^3,
    element radius in meter
    '''
    return (E*N.pi*r**2) / ls

def analyse_truss(jpos, jfreedom, loads, elements, E, r, p, weight_factor=1.0):
    ''' analyze given truss physically and return joint displacements, given:
    jpos: array of joint coordinates (2D)
    jfreedom: array of joint freedoms, example: [[True,True], [False,True]] (joint 0 can move in all directions, joint 1 only in y direction)
    loads: array of force vectors that act at the joints
    elements: index array showing which joint is connected to which other joint, 
        like this: [[0,1], [2,0], [3,0]] (means: joint 0 and 1 is connected, joint 2 and 0 connected, and so on)
    E: array with element's elasticity modules in N/m^2
    r: array with element radius
    p: array with element material density
    elements are assumed to have circular cross section 
    uses the direct stiffness method and returns:
        - joint displacements
        - element strains
        - elements stresses,
        - mass of the truss,
        - status of analysis,
        - timing values '''

    before = time.time()

    loads = loads.copy() # we gonna modify those

    jdof_indices = N.r_[:len(jpos)*2].reshape((-1,2))
    if N.any(loads[(~jfreedom[:,0]) & (~jfreedom[:,1])] > (0,0)):
        #print "warning: load at non-movable joint, setting load to zero"
        loads[(~jfreedom[:,0]) & (~jfreedom[:,1])] = (0,0)

    # element freedom table - for each element say to which joint indices it belongs
    eft = N.column_stack((
        jdof_indices[elements[:,0],0], jdof_indices[elements[:,0], 1], 
        jdof_indices[elements[:,1],0], jdof_indices[elements[:,1], 1] ))

    # calculate element (edge) metrics
    ls = elem_lengths(jpos, elements)
    masses = elem_masses(ls, r, p)
    weights = masses * gravity
    # add weights of elements to load
    for i, weight in enumerate(weights):
        loads[elements[i]] += (0, (-weight/2.)*weight_factor)
    loads[~jfreedom] = 0
    k = elem_stiffness(ls, E, r)
    # normalized edge direction vectors
    dxdy = (jpos[elements[:,1]] - jpos[elements[:,0]]) / ls[:, N.newaxis]
    s = dxdy[:,1] # y component of normalized edge - sin
    c = dxdy[:,0] # x component of normalized edge - cos
    # TODO: the next intermediate step may not be needed, could collect complete global stiffness matrix in one step maybe?
    # setup local stiffness matrices
    c2 = c**2
    sc = s*c
    s2 = s**2
    elem_K = N.column_stack((
                 c2,  sc, -c2, -sc,
                 sc,  s2, -sc, -s2,
                -c2, -sc,  c2,  sc,
                -sc, -s2,  sc,  s2)). reshape((-1,4,4))
    elem_K *= k[:, N.newaxis, N.newaxis]
    # global stiffness matrix K
    K = N.zeros((len(jpos)*2, len(jpos)*2), N.double)
    # combine local element stiffness matrices according to element freedom index table (eft)
    # TODO: this possible with numpy/indexing somehow?
    for elem,k in enumerate(elem_K):
        for (i,j),v in N.ndenumerate(k):
            K[eft[elem,i], eft[elem,j]] += v

    # degree of freedom pointer matrix
    pdof = jdof_indices[~jfreedom]
    # build modified stiffness matrix
    # this cancels out rows for which there is no movement freedom and puts a 1 in the diagonal entry of that row
    # Kmod will not be singular unless the truss system is really under rigid motion
    Kmod = K.copy()
    Kmod[pdof,:] = 0
    Kmod[:,pdof] = 0 # TODO: it works without this line. why?
    Kmod[pdof,pdof] = 1

    f = loads.ravel().astype(N.double)
    time_preparation = time.time() - before
    before = time.time()
    # solve Kmod * u = f
    try:
        u = N.linalg.solve(Kmod, f)
        solution_found = True
    except N.linalg.LinAlgError, e:
        solution_found = False
        status = 'solution failed (%s)' % e.message
        u = N.zeros_like(f)
        jpos_d = jpos.copy()
        strains = N.zeros(len(elements), N.double)
        stresses = N.zeros(len(elements), N.double)

    time_solution = time.time() - before

    before = time.time()
    if solution_found:
        jpos_d = jpos + u[jdof_indices]
        new_edge = jpos_d[elements[:,1]] - jpos_d[elements[:,0]]
        new_edge_lengths = N.sqrt(new_edge[:,0]**2 + new_edge[:,1]**2)
        # strain (german: dehnung) = change in length / original length
        strains = new_edge_lengths / ls - 1
        # stress (german: spannung)
        stresses = E * strains

        # detect rigid motion of joints simply by testing if very large displacements occured
        if N.abs(u).max() > N.abs(jpos).max() * 1000:
            status = 'rigid motion'
        else:
            status = 'success'

    # TODO: do i have to check if system is in force equilibrium?
    time_postprocessing = time.time() - before

    return u.reshape((-1,2)), strains, stresses, N.sum(masses), status, \
            (time_preparation, time_solution, time_postprocessing)


def visual_test(jpos, jfreedom, loads, elements, E, r, p, max_stress=3e+8):
    displacements, strains, stresses, mass, status, times = analyse_truss(jpos, jfreedom, loads, elements, E, r, p)
    jpos_d = jpos + displacements*100
    print status
    print times
    print 'displacements: ', displacements
    print 'strains: ', strains
    print 'stresses: ', stresses
    #strains_abs = N.abs(strains)

    from enthought.tvtk.api import tvtk
    from enthought.tvtk.tools import ivtk
    from enthought.pyface.api import GUI
    v = ivtk.viewer(False, False)
    v.scene.z_plus_view()

    pd = tvtk.PolyData()
    pts = jpos[elements].reshape((-1,2))
    pd.points = N.column_stack((pts[:,0], pts[:,1], N.zeros(len(pts))))
    pd.lines = N.r_[:len(elements)*2].reshape((-1,2))
    pd.cell_data.scalars = -strains
    pd.point_data.scalars = N.column_stack((r, r)).ravel()
    #tubes = tvtk.TubeFilter(input=pd, radius=N.sqrt(element_A.max() / N.pi), number_of_sides=16)
    tubes = tvtk.TubeFilter(input=pd, number_of_sides=16, vary_radius='vary_radius_by_absolute_scalar', capping=True)
    #tubes = tvtk.RibbonFilter(input=pd, use_default_normal=True, vary_width=True)
    b = tvtk.Actor(mapper=tvtk.PolyDataMapper(input=tubes.output, scalar_range=(-N.abs(strains).max(), N.abs(strains).max()), scalar_mode='use_cell_data'))
    b.mapper.lookup_table.hue_range = (0, 0.66)
    v.scene.add_actor(b)

    pd1 = tvtk.PolyData()
    pd1.points = N.column_stack((jpos_d[:,0], jpos_d[:,1], N.zeros(len(jpos))))
    pd1.lines = elements
    tubes1 = tvtk.TubeFilter(input=pd1, radius=0.01, number_of_sides=16)
    a = tvtk.Actor(mapper=tvtk.PolyDataMapper(input=tubes1.output))
    a.property.opacity = 0.3
    v.scene.add_actor(a)

    print "strain: ", strains.min(), strains.max()

    v.scene.reset_zoom()
    #ipshell = IPython.Shell.IPShellWX( [ '-wthread' ], user_ns=locals() )
    #ipshell.mainloop()
    GUI().start_event_loop()


if __name__ == '__main__':
    #   j6--j7--j8
    #  / |\/| \/| \
    # /  |/\| /\|  \
    #j1--j2-j3--j4--j5
    #jpos = N.array( [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0], [4.0, 0.0], [0.0, 1.0], [2.0, 1.0], [4.0, 1.0]] )
    #jfreedom = N.array( [[False, False], [True, True], [True, True], [True, True], [True, False], [True, True], [True, True], [True, True]] )
    #elements = N.array( [[0, 5], [0, 1], [1, 5], [1, 2], [1, 6], [2, 5], [2, 3], [2, 6], [2, 7], [3, 6], [3, 7], [3, 4], [4, 7], [5, 6], [6, 7]] )
    #loads = N.array( [[0.0, 0.0], [0.0, -5.0], [0.0, -5.0], [0.0, -5.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]] )
    #element_E = N.zeros( len(elements), N.double)
    #element_E[:] = 20000
    #element_A = N.zeros( len(elements), N.double)
    #element_A[:] = 0.01
    ##element_eas[:] = 2.1e+11 * 0.001 # 210 kN/mm^2 == 2.1e+11 N/m^2, A = 1mm^2 = 0.001m^2
    #visual_test(jpos, jfreedom, loads, elements, element_E, element_A)

    # another test from the textbook
    #jpos = N.array(     [[0.0, 0.0],     [1000.0, 0.0],   [2000.0, 0.0],   [3000.0, 0.0],   [4000.0, 0.0]])
    #jfreedom = N.array( [[False, False], [True, False], [True, False], [True, False], [True, False]] )
    #elements = N.array([[0,1], [1,2], [2,3], [3,4]])
    ##loads = N.array([[0,0], [0,1000], [0,1000], [0,1000], [0,1000]], N.float)
    #loads = N.array([[0,0], [1000,0], [1000,0], [1000,0], [1000,0]], N.float)
    #element_E = N.zeros( len(elements), N.double)
    #element_E[:] = 200000
    #element_A = N.zeros( len(elements), N.double)
    #element_A[:] = 100
    #visual_test(jpos, jfreedom, loads, elements, element_E, element_A)

    #jpos = N.array([[0.0, 0.0],     [2.0, -1.0],   [0.0, -1.0]])
    jpos = N.array([[-.1, 0.0],     [.1,0.0],   [0.0, -1.0]])
    jfreedom = N.array( [[False, False], [False, False], [True, True]] )
    elements = N.array([[1,2], [0,2]])
    #loads = N.array([[0,0], [0,1000], [0,1000], [0,1000], [0,1000]], N.float)
    p = 7700. # unit: kg/m^3
    r = N.array([0.001, 0.001])
    edge = jpos[elements[:,1]] - jpos[elements[:,0]]
    ls = N.sqrt(edge[:,0]**2 + edge[:,1]**2)
    #weights = p * N.pi * r**2 * ls * g # unit: Newton
    loads = N.array([[0,0], [0, 0], [0,-2000]], N.float)
    #loads = N.zeros_like(jpos)
    #for i,weight in enumerate(weights):
    #    loads[elements[i]] += (0, -weight/2.)
    #loads[~jfreedom] = 0
    element_E = N.zeros( len(elements), N.double)
    #element_E = N.array([2.1e+11, 0.1])
    element_E[:] =2.1e+11
    #element_A = N.zeros( len(elements), N.double)
    #element_A[:] = N.pi * r**2
    visual_test(jpos, jfreedom, loads, elements, element_E, r, 7700.0)

    #jpos = N.array([[0.0, 0.0],     [0.0, -2.0]])
    #jfreedom = N.array( [[False, False], [True, True]] )
    #elements = N.array([[0,1]])
    ##loads = N.array([[0,0], [0,1000], [0,1000], [0,1000], [0,1000]], N.float)
    #p = 7700. # unit: kg/m^3
    #r = N.array([0.0005])
    #loads = N.array([[0,0], [0,-250]], N.float)
    #element_E = N.zeros( len(elements), N.double)
    #element_E[:] = 2.1e+11
    ##element_A = N.zeros( len(elements), N.double)
    ##element_A[:] = N.pi * r**2
    #visual_test(jpos, jfreedom, loads, elements, element_E, r, 7700.0)



