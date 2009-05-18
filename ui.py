#! /usr/bin/env python

# Copyright 2009 Thomas Neumann
#
# Redistribution of this file is permitted under
# the terms of the GNU Public License (GPL) version 2.

from enthought.traits.api import *
from enthought.traits.ui.api import *
from enthought.pyface.api import ImageResource
from enthought.traits.ui.tabular_adapter import TabularAdapter
from enthought.traits.ui.menu import OKButton
import model

def general_edit_view(available_element_materials):
    return View(
        Tabbed(
            Group(
                Group(
                    HGroup(
                        Item(name='open_construction', label='Load', show_label=False, enabled_when='interaction_mode == "select"'),
                        Item(name='save_construction', label='Save', show_label=False, enabled_when='interaction_mode == "select"'),
                        Item(name='new_construction', label='Clear', show_label=False, enabled_when='interaction_mode == "select"'),
                    ),
                    Group(
                        Item(name='enter_add_joint_mode', label='Add new Joint', show_label=False, 
                            enabled_when='interaction_mode == "select"'),
                        Item(name='snap_joints', label='Snap joint coordinates', show_label=False, 
                            enabled_when='interaction_mode == "select"'),
                    ),
                    Group(
                        Item(name='new_element_material', editor=EnumEditor(values=_material_select(available_element_materials)), show_label=False),
                        Item(name='build_full_connections', label='Connect all joints', show_label=False, 
                            enabled_when='interaction_mode == "select"'),
                        label='Connect all Joints', show_border=True,
                    ),
                    Group(
                        Item(name='edit_available_materials', show_label=False),
                    ),
                    HGroup(
                        Item(name='width', object='object.construction'),
                        Item(name='height', object='object.construction'),
                        label='Bounds', show_border=True,
                    ),
                    Group(
                        Item(name='joint_force_magnitude', label='Force Magnitude (N)', object='object.construction'),
                        #Item(name='weight_factor', label='Weight Factor', object='object.construction'),
                    ),
                    label='Construction',
                ),
                Group(
                    Item(name='show_simulation', label='simulate'),
                    Item(name='status', label='simulation status', style='readonly', object='object.construction.simulation_result'),
                    Item(name='time', style='readonly', object='object.construction.simulation_result',
                        label='time taken (ms)', format_str='%.5f'),
                    #Item(name='time_preparation_percent', format_str='%.2f',
                    #    style='readonly', object='object.construction.simulation_result', label='time preparation (%)'),
                    #Item(name='time_solution_percent', format_str='%.2f',
                    #    style='readonly', object='object.construction.simulation_result', label='time solution (%)'),
                    #Item(name='time_postprocessing_percent', format_str='%.2f',
                    #    style='readonly', object='object.construction.simulation_result', label='time postprocess (%)'),
                    #Item(name='time_conversion_percent', format_str='%.2f',
                    #    style='readonly', object='object.construction.simulation_result', label='time conversion (%)'),
                    label='Simulation',
                ),
                Group(
                    Item(name='show_scalar', label='visualize'),
                    #Item(name='auto_scalar_range', label='auto max strain', enabled_when='show_scalar == "strain"'),
                    #Item(name='scalar_range_max', label='max strain (%)', enabled_when='not auto_scalar_range',
                    #    editor=RangeEditor(low=0, high=100, mode='slider')),
                    Item(name='displacement_amplify_factor', label='displacement amplification'),
                    Item(name='amplify_radius', label='radius amplification'),
                    label='Visualization',
                ),
                label = 'Editor',
            ),
        ),
        Tabbed(
            Group(
                Group(
                    Group(
                        Item(name='start_evolve', label="evolve!", editor=ButtonEditor(), show_label=False,
                            enabled_when='not interaction_mode == "evolve"'),
                        Item(name='pause_evolve', label="pause evolution", editor=ButtonEditor(), show_label=False,
                            enabled_when='interaction_mode == "evolve"'),
                        Item(name='reset_evolve', label="reset evolution", editor=ButtonEditor(), show_label=False,
                            enabled_when='ga.inited and not interaction_mode == "evolve"'),
                    ),
                    Group(
                        Item(name='num_steps', label='Generations', object='object.ga', style='readonly'),
                        Item(name='num_keep_generations', label='Keep generations', object='object.ga'),
                        Item(name='current_best_raw_fitness', label='best', 
                            object='object.ga', style='readonly'),
                        Item(name='current_mean_raw_fitness', label='mean', 
                            object='object.ga', style='readonly'),
                        Item(name='current_worst_raw_fitness', label='worst', 
                            object='object.ga', style='readonly'),
                    ),
                    label = 'Status & Control',
                ),
                Group(
                    Group(
                        Item(name='show_evolution_progress', label='Render Individual'),
                        Item(name='selected_population', label='Population', editor=EnumEditor(name='selectable_populations'),
                            enabled_when='show_evolution_progress'),
                        Item(name='selected_individual', label='Individual', editor=EnumEditor(name='selectable_individuals'),
                            enabled_when='show_evolution_progress'),
                    ),
                    label = 'Browse',
                ),
                label = 'Evolution',
            ),
        ),
        Tabbed(
            Group(
                Group(
                    Group(
                        Item(name='population_size', object='object.ga'),
                        Item(name='crossover_rate', object='object.ga'),
                        Item(name='mutation_rate', object='object.ga'),
                    ),
                    Group(
                        Item(name='maximum_offset', label='Joint Mutate Offset', object='object.ga.mutation_operator.joint_position_mutation'),
                        label = 'Mutation',
                    ),
                    Group(
                        Item(name='arena_size', object='object.ga.selection_operator'),
                        label = 'Selection',
                    ),
                    Group(
                        Item(name='unfeasible_constraint_penalty', label='Unfeasible penalty', object='object.ga.survival'),
                        HGroup(
                            Item(name='force_variety_topology', label='Enabled', object='object.ga.survival'),
                            Item(name='topology_variety_distance', label='Tolerance', object='object.ga.survival'),
                            label='Enforce topology diversity', show_border=True,
                        ),
                        HGroup(
                            Item(name='force_variety_shape', label='Enabled', object='object.ga.survival'),
                            Item(name='shape_variety_distance', label='Tolerance', object='object.ga.survival'),
                            label='Enforce shape diversity', show_border=True,
                        ),
                        HGroup(
                            Item(name='force_variety_size', label='Enabled', object='object.ga.survival'),
                            label='Enforce size diversity', show_border=True,
                        ),
                        Item(name='num_similar_individuals', label='Num similar allowed', object='object.ga.survival'),
                        label = 'Fitness',
                    ),
                    Group(
                        Item(name='deleted_material_gene_places', label='Material gene deletion places', object='object.ga.world'),
                        label = 'Genes',
                    ),
                ),
                label = 'GA Basics',
            ),
        ),
        Tabbed(
            Group(
                Group(
                    Group(
                        Item(name='joint_offset', label='Joint Init Offset', object='object.ga.init_individual', ),
                        Item(name='copy_probability', object='object.ga.init_individual', ),
                        Item(name='mutate_probability',object='object.ga.init_individual', ),
                        label = 'Initialization',
                    ),
                    Group(
                        Item(name='check_joint_overlap', label='Kill when joint overlap', object='object.ga.validator'),
                        Item(name='joint_overlap_tolerance', label='Joint overlap tolerance', object='object.ga.validator'),
                        Item(name='check_simulatable', label='Kill when not simulatable', object='object.ga.validator'),
                        Item(name='kill_stability_constraint', label='Kill when stress constraint', object='object.ga.validator'),
                        Item(name='kill_displacement_constraint', label='Kill when displacement constraint', object='object.ga.validator'),
                        label = 'Validation/Birth',
                    ),
                ),
                label = 'GA Advanced',
            ),
        ),
    )


def joint_edit_view(w, h):
    return View(
        Group(
            Group(
                Item(label='X', name='x', object='joint', editor=RangeEditor(low=-w/2., high=w/2., mode='slider')),
                Item(label='Y', name='y', object='joint', editor=RangeEditor(low=-h/2., high=h/2., mode='slider')),
                label='Position', show_border=True,
            ),
            HGroup(
                Item(name='movable_x', label='X', object='joint'),
                Item(name='movable_y', label='Y', object='joint'),
                label='Movable in simulation', show_border=True,
            ),
            Group(
                Item(label='X', name='max_displacement_x', object='joint', editor=RangeEditor(low=-w/2., high=w/2., mode='slider')),
                Item(label='Y', name='max_displacement_y', object='joint', editor=RangeEditor(low=-h/2., high=h/2., mode='slider')),
                label='Maximal Displacement', show_border=True,
            ),
            Group(
                Item(name='force_x', label='X', object='joint', enabled_when='joint.movable_xy'),
                Item(name='force_y', label='Y', object='joint', enabled_when='joint.movable_xy'),
                label='Force', show_border=True,
            ),
            Group(
                Item(name='mutatable_x', object='joint'),
                Item(name='mutatable_y', object='joint'),
                label = 'Genetic Algorithm Properties', show_border=True,
            ),
            Group(
                Item(name='remove_selected_object', label='Remove this Joint', show_label=False),
                Item(name='enter_add_element_mode', show_label=False, label='Connect to other joint'),
                label = 'Actions', show_border=True
            ),
            orientation='vertical', label='Edit Joint'
        )
    )

def element_edit_view(available_materials):
    return View(
        Group(
            Group(
                Item(name='material', editor=EnumEditor(values=_material_select(available_materials)), object='element'),
                Item(name='mutatable_material', object='element'),
                Item(name='deletable', label='Deletable by GA', object='element'),
            ),
            Group(
                Item(name='density', label='Density [kg/m^3]', object='element.material', style='readonly'),
                Item(name='maximum_stress_mpa', label='Maximum Stress [MPa]', object='element.material', style='readonly'),
                Item(name='radius', label='Radius [m]', object='element.material', style='readonly'),
                label = 'Material Info', show_border=True,
            ),
            Group(
                Item(name='remove_selected_object', label='Remove this Element', show_label=False),
            ),
            orientation='vertical', label='Edit Element'
        )
    )

class IndividualSelect(HasTraits):
    def __init__(self, individual, num_in_population):
        self.individual = individual
        self.num_in_population = num_in_population

    def __str__(self):
        return "#%04d (%d, %.2f/%.2f)" % (self.num_in_population, self.individual.age, self.individual.mass, self.individual.raw_fitness)

class PopulationSelect(HasTraits):
    def __init__(self, population):
        self.population = population

    def __str__(self):
        return "#%04d (best: %.2f, mean: %.2f)" % (self.population.number, self.population.best.raw_fitness, self.population.mean_fitness)

def _material_select(available_materials):
    material_select = {}
    for i, mat in enumerate(available_materials):
        material_select[mat] = '%02d:%.2fcm %s' % (i+1, mat.radius*100, mat.name)
    return material_select

class ElementMaterialTabularAdapter(TabularAdapter):
    columns = [('Name', 'name'),
            ('Radius (m)', 'radius'),
            ('Elasticity module','E'), 
            ('Maximum stress', 'maximum_stress'),
            ('Density', 'density') ]

    default_value = Property()
    def _get_default_value(self):
        return model.ElementMaterial()

def edit_available_elements_view():
    return View(
            Item( name='available_element_materials',
                  style='custom', show_label=False, 
                  editor = TabularEditor(adapter = ElementMaterialTabularAdapter(),
                      operations=['delete','append','edit']),),
            height=360, width=600, title='Edit available materials', buttons=[OKButton],
            )

