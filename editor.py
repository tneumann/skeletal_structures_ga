#! /usr/bin/env python

# Copyright 2009 Thomas Neumann
#
# Redistribution of this file is permitted under
# the terms of the GNU Public License (GPL) version 2.

from enthought.traits.api import *
from enthought.tvtk.api import tvtk
import enthought.pyface.api as pyface
from enthought.tvtk.pyface.api import Scene
import enthought.traits.ui.api as tui
from enthought.chaco.api import Plot, ArrayPlotData
import numpy as N
import wx
import pickle
import sys

import physics
import ga as GA
import ga_truss as GAT
import data
import model
import ui

#def fitness_for_construction(ga, construction):
#    w = GAT.World(construction = construction)
#    return ga.fitness_function(w, w.individual_from_construction(construction))

def pts2dto3d(pts):
    if len(pts) == 0:
        return None
    return N.column_stack((pts[:,0], pts[:,1], N.zeros(len(pts), pts.dtype)))


class ConstructionEditor(pyface.SplitApplicationWindow):

    def __init__(self, **kwargs):
        self._handle_by_joint = {}
        super(ConstructionEditor, self).__init__(**kwargs)

#############################################################
# visualization related
    scene = Instance(Scene) # vtk scene
    show_simulation = Bool(True)
    show_scalar = Enum('stress', 'strain')
    scalar_range_max = Float(50.)
    auto_scalar_range = Bool(True)
    displacement_amplify_factor = Float(1000.)
    amplify_radius = Float(1.5)
    show_cross_sectional_areas_labels = Bool(True)

    _on_init = Event() # internal event that is triggered when scene has been setup

    _bg_plane_actor = None
    _axis = None
    _bg_plane_picker = None
    _bg_plane_width = None
    _bg_plane_height = None

    _reset_zoom_needed = Bool(False)

    def reset_zoom(self):
        if self.scene:
            self.scene.reset_zoom()

    # schedule a render when model changes
    @on_trait_change('construction.[joints,elements].+, construction.[joint_force_magnitude,width,height,weight_factor], selected_object, displacement_amplify_factor')
    def render_if_not_interacting(self):
        if self.scene and not self.scene._interacting:
            self.render_interacting()

    def render_interacting(self):
        if self.scene:
            old_interacting = self.scene._interacting
            self.scene._interacting = True
            self.scene.render()
            self.scene._interacting = old_interacting

    def _setup_scene(self):
        # scalar bar for strain
        lut_strain = tvtk.LookupTable(hue_range=(0.66, 0.0))
        lut_strain.build()
        self._scalar_bar_strain = tvtk.ScalarBarActor(lookup_table=lut_strain, orientation='horizontal',
                text_position='succeed_scalar_bar', maximum_number_of_colors=256, number_of_labels=9,
                position=(0.1, 0.01), position2=(0.8, 0.08), title='element strain (%)', visibility=False)
        self.scene.add_actor(self._scalar_bar_strain)
        # scalar bar for stress
        # lookup table from green to yellow, and last 2 values dark red
        lut_stress = tvtk.LookupTable(hue_range=(0.33, 0.1), number_of_table_values=256)
        lut_stress.build()
        lut_stress.set_table_value(254, (1.0, 0.4, 0.0, 1.0))
        lut_stress.set_table_value(255, (1.0, 0.0, 0.0, 1.0))
        self._scalar_bar_stress = tvtk.ScalarBarActor(lookup_table=lut_stress, orientation='horizontal',
                text_position='succeed_scalar_bar', maximum_number_of_colors=256, number_of_labels=9,
                position=(0.1, 0.01), position2=(0.8, 0.08), title='element stress', visibility=False)
        self.scene.add_actor(self._scalar_bar_stress)
        # setup elements visualization
        self._elements_polydata = tvtk.PolyData()
        self._tubes = tvtk.TubeFilter(input=self._elements_polydata, number_of_sides=6, 
                vary_radius='vary_radius_by_absolute_scalar', radius_factor=1.5)
        mapper = tvtk.PolyDataMapper(input=self._tubes.output, lookup_table=lut_strain, 
                interpolate_scalars_before_mapping=True, scalar_mode='use_cell_data')
        self._elements_actor = tvtk.Actor(mapper=mapper)
        self.scene.add_actor(self._elements_actor)
        # show elements in deformed state as wireframe
        self._deformed_elements_polydata = tvtk.PolyData()
        self._deformed_elements_actor = tvtk.Actor(mapper=tvtk.PolyDataMapper(input=self._deformed_elements_polydata))
        self._deformed_elements_actor.property.set(opacity = 0.2, representation = 'wireframe')
        self.scene.add_actor(self._deformed_elements_actor)
        # highlight one element via a ribbon outline
        self._hl_element_ribbons = tvtk.RibbonFilter(input=tvtk.PolyData(), use_default_normal=True, width=1.0)
        self._hl_element_actor = tvtk.Actor(mapper=tvtk.PolyDataMapper(input=self._hl_element_ribbons.output), visibility=False)
        self._hl_element_actor.property.set(ambient=1, ambient_color=(1,1,1), diffuse=0)
        self.scene.add_actor(self._hl_element_actor)
        # cross sectional radius labels
        self._elements_label_polydata = tvtk.PolyData()
        self._label_cellcenters = tvtk.CellCenters(input=self._elements_label_polydata)
        self._label_visps = tvtk.SelectVisiblePoints(renderer=self.scene.renderer, input=self._label_cellcenters.output, tolerance=10000)
        self._label_actor = tvtk.Actor2D(mapper=tvtk.Dynamic2DLabelMapper(input=self._label_visps.output, label_mode='label_scalars'))
        self._label_actor.mapper.label_text_property.set(bold=True, italic=False, justification='centered', font_size=14)
        #self.scene.add_actor(self._label_actor)
        # force glyphs (use arrows for that)
        self._force_glyphs = tvtk.Glyph3D(scale_mode='scale_by_vector', vector_mode='use_vector', 
                color_mode='color_by_vector', scale_factor=1.)
        self._force_glyphs.set_source(0, tvtk.ArrowSource(shaft_radius=0.04, tip_resolution=8, tip_radius=0.2).output)
        self._force_polydata = tvtk.PolyData()
        self._force_glyphs.set_input(0, self._force_polydata)
        self._force_actor = tvtk.Actor(mapper=tvtk.PolyDataMapper(input=self._force_glyphs.output, scalar_range=(0,10)))
        self._force_actor.mapper.lookup_table.hue_range = (0.33, 0.0)
        self.scene.add_actor(self._force_actor)
        # current status display
        self._text_actor = tvtk.TextActor(position=(0.5,0.95))
        self._text_actor.position_coordinate.coordinate_system = 'normalized_display'
        self._text_actor.text_property.set(font_size=14, justification='center')
        self.scene.add_actor(self._text_actor)
        # a nice gradient background
        self.scene.renderer.set(background2=(0.28, 0.28, 0.28), background=(0.01, 0.01, 0.02), gradient_background=True)
        # setup events
        self.interaction_mode = 'select'
        self.scene.interactor.add_observer('MouseMoveEvent', self._mouse_move)
        self.scene.interactor.add_observer('LeftButtonPressEvent', self._mouse_press)
        self.scene.renderer.add_observer('StartEvent', self._before_render)
        self._on_init = True
        self._reset_zoom_needed = True

    def bw_mode(self):
        self.scene.renderer.set(background = (1,1,1), gradient_background = False)
        self._deformed_elements_actor.visibility=False
        self._elements_actor.mapper.scalar_visibility=False
        self._elements_actor.property.set(ambient=1, ambient_color=(.3,.3,.3), diffuse_color=(0,0,0))
        self._axis.visibility = False
        #self._tubes.set(vary_radius='vary_radius_off', radius=0.1)
        self._bg_plane_actor.visibility = False
        self._text_actor.property.color = (0,0,0)
        self._label_actor.mapper.label_text_property.set(color=(0,0,0))

    def pic_mode(self):
        self.scene.renderer.set(background = (1,1,1), gradient_background = False)
        self._deformed_elements_actor.visibility=False
        self._axis.visibility = False
        self._bg_plane_actor.visibility = False
        self._text_actor.property.color = (0,0,0)
        self._label_actor.visibility = False
        self._scalar_bar_stress.lookup_table.set(hue_range=(0,0), value_range=(0,1))
        self._scalar_bar_stress.lookup_table.force_build()
        self.scene.reset_zoom()

    def _before_render(self, *args):
        self._redraw_background_plane()
        self._redraw_joint_handles()
        self._redraw_joint_forces()
        self._redraw_elements()
        self._redraw_element_labels()
        self._redraw_joints()
        self._redraw_caption()
        self._setup_elements_picker()
        self._show_scalar_bar_when_simulating()
        if self._reset_zoom_needed:
            self.reset_zoom()
            self._reset_zoom_needed = False

    def _redraw_background_plane(self):
        if self.construction and self.scene:
            if self._bg_plane_width == self.construction.width and self._bg_plane_height == self.construction.height:
                return
            if self._bg_plane_actor:
                self.scene.remove_actor(self._bg_plane_actor)
            if self._axis:
                self.scene.remove_actor(self._axis)
            w, h = self.construction.width, self.construction.height
            plane = tvtk.PlaneSource(x_resolution=int(w), y_resolution=int(h))
            scalation = tvtk.Transform()
            scalation.scale((w,h,1))
            scale_plane = tvtk.TransformPolyDataFilter(transform=scalation, input=plane.output)
            self._bg_plane_actor = tvtk.Actor(mapper=tvtk.PolyDataMapper(input=scale_plane.output))
            self._bg_plane_actor.property.set(representation='wireframe', line_stipple_pattern=0xF0F0, opacity=0.15)
            self.scene.add_actor(self._bg_plane_actor)
            self._axis = tvtk.CubeAxesActor2D(camera=self.scene.camera, z_axis_visibility=False, corner_offset=0, bounds=[-w/2, w/2, -h/2, h/2, 0, 0])
            self.scene.add_actor(self._axis)
            self._bg_plane_picker = tvtk.CellPicker(tolerance=0, pick_from_list=True)
            self._bg_plane_picker.pick_list.append(self._bg_plane_actor)
            self._bg_plane_width, self._bg_plane_height = self.construction.width, self.construction.height

    def _show_scalar_bar_when_simulating(self):
        if self.construction:
            essential = self.show_simulation and len(self.construction.elements) > 0
            self._scalar_bar_strain.visibility = essential and self.show_scalar == 'strain'
            self._scalar_bar_stress.visibility = essential and self.show_scalar == 'stress'

    def _representation_model_for_joint_handle(self, joint):
        if joint.movable_x and joint.movable_y:
            model = data.joint_model_movable
        elif not joint.movable_x and not joint.movable_y:
            model = data.joint_model_unmovable
        elif joint.movable_x and not joint.movable_y:
            model = data.joint_model_movable_x
        elif not joint.movable_x and joint.movable_y:
            model = data.joint_model_movable_y
        return model

    def _redraw_joint_handles(self):
        if self.construction:
            #print "rebuild_joint_handles"
            for joint in self.construction.joints:
                if not joint in self._handle_by_joint:
                    self._handle_by_joint[joint] = self._make_handle_for_joint(joint)

    def _redraw_joint_forces(self):
        #print "redraw_joint_forces"
        if self.construction:
            points = []
            vectors = []
            for joint in self.construction.joints:
                if joint.force_x != 0 or joint.force_y != 0:
                    points.append((joint.x, joint.y, 0))
                    vectors.append((joint.force_x, joint.force_y, 0))
            if len(points) > 0:
                self._force_polydata.points = points
                self._force_polydata.point_data.vectors = vectors
            else:
                self._force_polydata.points = None
                self._force_polydata.point_data.vectors = None

    def _redraw_elements(self):
        if self.construction:
            jp = self.construction.joint_positions[self.construction.element_index_table].reshape((-1,2))
            self._elements_polydata.points = pts2dto3d(jp)
            self._elements_polydata.lines = N.r_[:len(self.construction.elements)*2].reshape((-1,2)) # => [[0,1],[2,3],[4,5],...]
            radii = N.array([e.material.radius for e in self.construction.elements])
            self._elements_polydata.point_data.scalars = N.column_stack((radii, radii)).ravel() * self.amplify_radius
            # set scalars of elements when enabled by the user
            if self.show_simulation and self.construction.simulated and len(self.construction.elements) > 0:
                # show element stress
                if self.show_scalar == 'stress':
                    self._elements_actor.mapper.lookup_table = self._scalar_bar_stress.lookup_table
                    breakage = (N.abs(self.construction.simulation_result.element_stresses) / self.construction.max_element_stress_array) * 100. # percentual breakage
                    self._elements_polydata.cell_data.scalars = breakage
                    self._scalar_bar_stress.lookup_table.table_range = \
                            self._elements_actor.mapper.scalar_range = (0, 100)
                # show element strain
                elif self.show_scalar == 'strain':
                    self._elements_actor.mapper.lookup_table = self._scalar_bar_strain.lookup_table
                    strains = self.construction.simulation_result.element_strains
                    self._elements_polydata.cell_data.scalars = strains
                    # set strain scalar range
                    if self.auto_scalar_range:
                        max_strain = N.abs(strains).max()
                        scalar_range = (-max_strain, max_strain)
                    else:
                        scalar_range = (-self.scalar_range_max, self.scalar_range_max)
                    self._scalar_bar_strain.lookup_table.table_range = \
                        self._elements_actor.mapper.scalar_range = scalar_range
            else:
                self._elements_polydata.cell_data.scalars = None
            # deformed wireframe model:
            if self.show_simulation and self.construction.simulated and len(self.construction.joints) > 0 \
                    and self.construction.simulation_result.okay:
                # amplify displacements so they become visible
                jp = self.construction.joint_positions + \
                        self.construction.simulation_result.joint_displacements * self.displacement_amplify_factor
                self._deformed_elements_polydata.points = pts2dto3d(jp)
                self._deformed_elements_polydata.lines = self.construction.element_index_table
                self._deformed_elements_actor.visibility = True
            else:
                self._deformed_elements_actor.visibility = False
            # redraw selected element highlight
            if isinstance(self.selected_object, model.Element):
                self._hl_element_ribbons.input = tvtk.PolyData(\
                        points=[list(self.selected_object.joint1.position)+[0], list(self.selected_object.joint2.position)+[0]], lines=[[0,1]])
                self._hl_element_ribbons.width = max(self.selected_object.material.radius, 0.05) * 2.0

    def _redraw_element_labels(self):
        if self.construction:
            eit = self.construction.element_index_table
            jp = self.construction.joint_positions
            self._elements_label_polydata.points = pts2dto3d((jp[eit[:,0]] + jp[eit[:,1]]) / 2)
            self._elements_label_polydata.point_data.scalars = N.round(N.array([e.material.radius for e in self.construction.elements]) * 100, 3)
            self._label_visps.input = self._elements_label_polydata
            self._label_actor.mapper.input = self._label_visps.output

    def _redraw_caption(self):
        if self.construction:
            if len(self.construction.joints) > 0 and len(self.construction.elements) > 0:
                # calculate fitness for shown construction
                #self.construction_fitness = fitness_for_construction(self.ga, self.construction)
                if self.construction.simulation_result.okay:
                    stability = self.construction.simulation_result.stability
                    if stability < 0:
                        stability_txt = 'STATICALLY UNSTABLE (%.0f%%)' % (stability * 100)
                        self._text_actor.text_property.color = (1.0, 0.0, 0.0)
                    else:
                        self._text_actor.text_property.color = (1.0, 1.0, 1.0)
                        stability_txt = 'Stability %.0f%%' % (stability * 100)
                    self._text_actor.input = 'Weight: %.2fkg, %s' % (self.construction.simulation_result.construction_mass, stability_txt)
                else:
                    self._text_actor.input = 'Simulation Error (%s)' % self.construction.simulation_result.status
                self._text_actor.visibility = True
            else:
                self._text_actor.visibility = False

    def _redraw_joints(self):
        if self.scene and self.construction:
            #print "redraw joints"
            # a joint might have been removed - the handle for that should be removed
            joint_handle_pairs_to_remove = []
            for joint, handle in self._handle_by_joint.iteritems():
                if not joint in self.construction.joints:
                    joint_handle_pairs_to_remove.append((joint, handle))
            for joint, handle in joint_handle_pairs_to_remove:
                del self._handle_by_joint[joint]
                handle.off()
            for joint in self.construction.joints:
                handle = self._handle_by_joint[joint]
                handle.representation.handle = self._representation_model_for_joint_handle(joint)
                handle.representation.world_position = (joint.x, joint.y, 0) 



#############################################################
# gui related
    gui = Instance(pyface.GUI)
    property_view_contents = Any()
    title = String('Plane Truss Design via Genetic Algorithms')
    ratio = Float(1.0)
    direction = Str('vertical')
    construction_fitness = Float(0.0)

    def _create_splitter(self, parent):
        s = super(ConstructionEditor, self)._create_splitter(parent)
        self._setup_scene()
        return s

    def _create_lhs(self, parent):
        self.scene = Scene(parent)
        self.scene.interactor.interactor_style = tvtk.InteractorStyleImage()
        #self.scene.camera.parallel_projection = True
        self.scene.z_plus_view()
        return self.scene.control

    def _create_rhs(self, parent):
        self.properties_panel = wx.Panel(parent)
        parent.SetMinimumPaneSize(0)
        parent.SetSashSize(0)
        parent.SetSize((0,100))
        self._splitter_control = parent
        return self.properties_panel

    def _property_view_contents_changed(self, old, new):
        if old:
            old.dispose()
        if new:
            w, h = new.control.GetSize()
            ps = self._splitter_control.GetMinimumPaneSize()
            self._splitter_control.SetMinimumPaneSize(max(ps,w))
            new.control.SetSize((max(ps,w), h))
        else:
            self._splitter_control.SetMinimumPaneSize(0)

    def close(self):
        if self.scene is not None:
            self.scene.close()
        super(ConstructionEditor, self).close()

    @on_trait_change('show_evolution_progress, ga.on_step')
    def _update_evolve_gui(self):
        if self.show_evolution_progress and self.ga.current_population:
            self.selectable_individuals = [ui.IndividualSelect(individual, i+1) for i,individual in enumerate(self.ga.current_population.individuals)]
            self.selected_individual = self.selectable_individuals[0]

#############################################################
# plots
    plot_num_generations = Enum(-1, 10, 50, 500, 5000)

    fitness_plot = Instance(Plot)
    fitness_plot_data = Instance(ArrayPlotData)
    #feasible_plot = Instance(Plot)
    #feasible_plot_data = Instance(ArrayPlotData)
    #stability_plot = Instance(Plot)
    #stability_plot_data = Instance(ArrayPlotData)

    @on_trait_change('_on_init')
    def _init_plots(self):
        # fitness
        self.fitness_plot_data = ArrayPlotData(generation=N.zeros(1), best=N.zeros(1), average=N.zeros(1))
        self.fitness_plot = Plot(self.fitness_plot_data)
        self.fitness_plot.legend.visible = True
        self.fitness_plot.set(padding_top = 5, padding_right = 5, padding_bottom = 20, padding_left = 40)
        self.fitness_plot.plot(('generation', 'best'), color='green', line_width=2, name='best')
        self.fitness_plot.plot(('generation', 'average'), color='black', name='avg')
        ## stability
        #self.stability_plot_data = ArrayPlotData(generation=N.zeros(1), min=N.zeros(1), max=N.zeros(1), average=N.zeros(1))
        #self.stability_plot = Plot(self.stability_plot_data, height=100)
        #self.stability_plot.legend.visible = True
        #self.stability_plot.set(padding_top = 15, padding_right = 5, padding_bottom = 20, padding_left = 40, title='Stability')
        #self.stability_plot.plot(('generation', 'min'), color='red', line_width=2, name='min')
        #self.stability_plot.plot(('generation', 'average'), color='black', name='avg')
        #self.stability_plot.plot(('generation', 'max'), color='green', line_width=2, name='max')
        ## feasible
        #self.feasible_plot_data = ArrayPlotData(generation=N.zeros(1), num=N.zeros(1))
        #self.feasible_plot = Plot(self.feasible_plot_data)
        #self.feasible_plot.set(padding_top = 15, padding_right = 5, padding_bottom = 20, padding_left = 40, title = 'Unfeasible Individuals')
        #self.feasible_plot.plot(('generation', 'num'), color='red', line_width=2)

    @on_trait_change('ga.on_init')
    def _reset_plots(self):
        self._fitness_best_history = []
        self._fitness_avg_history = []
        #self._stability_max_history = []
        #self._stability_min_history = []
        #self._stability_avg_history = []
        #self._num_feasible_history = []

    @on_trait_change('show_evolution_progress, ga.on_step, ga.on_init')
    def _update_plots(self):
        mm = self.plot_num_generations # just an alias so that slicing expressions do not become too big
        gens = N.r_[:self.ga.num_steps][-mm:]
        # fitness plot
        if self.ga.current_population:
            self._fitness_best_history.append(self.ga.current_population.best.raw_fitness)
            self._fitness_avg_history.append(N.array([i.raw_fitness for i in self.ga.current_population.individuals if i.feasible and not i.got_penalized]).mean())
        self.fitness_plot_data.set_data('generation', gens)
        self.fitness_plot_data.set_data('best', self._fitness_best_history[-mm:])
        self.fitness_plot_data.set_data('average', self._fitness_avg_history[-mm:])
        ## stability plot
        #if self.ga.current_population:
        #    stabilities = N.array([i.stability for i in self.ga.current_population.individuals], N.float) * 100
        #    self._stability_min_history.append(stabilities.min())
        #    self._stability_max_history.append(stabilities.max())
        #    self._stability_avg_history.append(stabilities.mean())
        #self.stability_plot_data.set_data('generation', gens)
        #self.stability_plot_data.set_data('min', self._stability_min_history[-mm:])
        #self.stability_plot_data.set_data('max', self._stability_max_history[-mm:])
        #self.stability_plot_data.set_data('average', self._stability_avg_history[-mm:])
        ## feasible plot
        #if self.ga.current_population:
        #    self._num_feasible_history.append(len([i for i in self.ga.current_population.individuals if not i.feasible or i.got_penalized]))
        #self.feasible_plot_data.set_data('generation', gens)
        #self.feasible_plot_data.set_data('num', self._num_feasible_history[-mm:])


#############################################################
# ga related
    ga = Instance(GA.GeneticAlgorithm)
    fitness_function = Callable()

    start_evolve = Button()
    pause_evolve = Button()
    reset_evolve = Button()
    evolving = Bool(False)

    selected_individual = Instance(ui.IndividualSelect)
    selectable_individuals = List(ui.IndividualSelect)

    show_evolution_progress = Bool(True)
    ga_step = Event

    #def _internal_on_ga_step(self):
    #    self.ga_step = True

    #def _ga_changed(self, old, new):
    #    if old != new and new:
    #        new.on_trait_event(self._internal_on_ga_step, 'on_step', dispatch='fast_ui')

    def _start_evolve_fired(self):
        if not self.ga.inited:
            self.ga.world.construction = self.construction
            self.construction = self.construction.clone_traits(copy='deep')
            self.ga.init_evolution()
            self.ga_step = True
        self.interaction_mode = 'evolve'
        self.gui.invoke_later(self._evolve_it)

    def _fitness_function_changed(self):
        self.ga.fitness_function = self.fitness_function
        self.pause_evolve = True

    @on_trait_change('selected_individual, show_evolution_progress')
    def _show_current_individual(self):
        if self.show_evolution_progress:
            self.construction = self.ga.world.construction_from_individual(self.selected_individual.individual)

    def _evolve_it(self):
        #for i in xrange(5):
        self.ga.evolution_step()
        self.gui.process_events()
        if self.interaction_mode == 'evolve':
            self.gui.invoke_later(self._evolve_it)

    def _pause_evolve_fired(self):
        self.interaction_mode = 'select'

    def _reset_evolve_fired(self):
        if self.ga.inited and self.ga.world:
            self.construction = self.ga.world.construction
            self.ga.reset_evolution()

#############################################################
# model
    construction = Instance(model.Construction)

    open_construction = Button()
    save_construction = Button()
    new_construction = Button()

    def _open_construction_fired(self):
        file_dialog = pyface.FileDialog(action='open',
                wildcard='Constructions (*.con)|*.con|')
        if file_dialog.open() == pyface.OK:
            self.construction, self.ga = pickle.load(open(file_dialog.path))

    def _save_construction_fired(self):
        file_dialog = pyface.FileDialog(action='save as',
                wildcard='Constructions (*.con)|*.con|')
        if file_dialog.open() == pyface.OK:
            pickle.dump((self.construction, self.ga), open(file_dialog.path, 'w'))

    def _new_construction_fired(self):
        self.construction = model.Construction(available_element_materials=data.steels, element_deleted_material=data.air)

    edit_available_materials = Button()

    def _edit_available_materials_fired(self):
        self.construction.edit_traits( view=ui.edit_available_elements_view() )

#############################################################
# editing
    enter_add_joint_mode = Button()
    enter_add_element_mode = Button()
    remove_selected_object = Button()
    build_full_connections = Button()
    snap_joints = Button()
    new_element_material = Instance(model.ElementMaterial)

    selected_object = Either(Instance(model.Joint), Instance(model.Element))
    interaction_mode = Enum('select', 'add_joint', 'add_element', 'evolve')

    _last_hovered_element = None

    def _enter_add_joint_mode_fired(self):
        self.selected_object = None
        self.interaction_mode = 'add_joint'

    def _enter_add_element_mode_fired(self):
        self.interaction_mode = 'add_element'

    def _remove_selected_object_fired(self):
        if isinstance(self.selected_object, model.Joint):
            self.construction.joints.remove(self.selected_object)
            self.selected_object = None
        elif isinstance(self.selected_object, model.Element):
            self.construction.elements.remove(self.selected_object)
            self.selected_object = None

    def _build_full_connections_fired(self):
        self.construction.build_full_connections(self.new_element_material)

    def _snap_joints_fired(self):
        for joint in self.construction.joints:
            joint.x = N.round(joint.x)
            joint.y = N.round(joint.y)

    @on_trait_change('selected_object, _on_init')
    def _edit_selected_object_properties(self):
        # no object selected, show general ui
        if self.selected_object == None:
            self.property_view_contents = \
                    self.edit_traits(kind='subpanel', parent=self.properties_panel, 
                            view=ui.general_edit_view([self.construction.element_deleted_material] + self.construction.available_element_materials))
        # joint selected
        elif isinstance(self.selected_object, model.Joint):
            w, h = self.construction.width, self.construction.height
            self.property_view_contents = \
                    self.edit_traits(view=ui.joint_edit_view(w, h), kind='subpanel', parent=self.properties_panel,
                        context={'joint': self.selected_object, 'object': self})
        # element selected
        elif isinstance(self.selected_object, model.Element):
            self.property_view_contents = \
                    self.edit_traits( view = ui.element_edit_view([self.construction.element_deleted_material] + self.construction.available_element_materials), 
                        kind = 'subpanel', parent = self.properties_panel,
                        context = {'element': self.selected_object, 'object': self} )
                        

    def _pick_element(self, x, y):
        self._element_picker.pick((float(x), float(y), 0.0), self.scene.renderer)
        if self._element_picker.cell_id > -1:
            return self.construction.elements[self._element_picker.cell_id]
        else:
            return None

    def _interaction_mode_changed(self):
        # need to disable the handle widgets?
        for h in self._handle_by_joint.values():
            h.process_events = self.interaction_mode in ['select', 'add_element']
        # cursor for adding
        if self.interaction_mode == 'add_joint':
            self.scene.render_window.current_cursor = 10

    def _mouse_press(self, *args):
        x, y = self.scene.interactor.last_event_position
        if self.interaction_mode == 'select':
            element = self._pick_element(x, y)
            if element:
                self.selected_object = element
            else:
                self.selected_object = None
        elif self.interaction_mode == 'add_joint':
            # did we hit the bounds plane?
            self._bg_plane_picker.pick((x,y,0), self.scene.renderer)
            if self._bg_plane_picker.cell_id > -1:
                wx, wy, wz = self._bg_plane_picker.pick_position
                self.construction.joints.append(model.Joint(x=wx, y=wy))
                self.interaction_mode = 'select'

    def _mouse_move(self, *args):
        if self.interaction_mode == 'select':
            x, y = self.scene.interactor.event_position
            element = self._pick_element(x, y)
            if element:
                self.scene.render_window.current_cursor = 9
                self._last_hovered_element = element
            elif self._last_hovered_element:
                self.scene.render_window.current_cursor = 0
                self._last_hovered_element = None

    def _get_handle_by_joint(self, joint):
        try:
            return self._handle_by_joint[joint]
        except KeyError:
            return None

    def _selected_object_changed(self, old, new):
        if old and isinstance(old, model.Joint):
            handle = self._get_handle_by_joint(old)
            if handle:
                handle.representation.property.set(ambient=0, diffuse=1)
        self._hl_element_actor.visibility = False
        if new:
            if isinstance(new, model.Joint):
                self._get_handle_by_joint(new).representation.property.set(ambient=1, diffuse=0.3)
                if self.interaction_mode == 'add_element' and old != new and isinstance(old, model.Joint):
                        elem = model.Element(joint1=old, joint2=new, material=self.new_element_material)
                        if self.construction.element_addable(elem):
                            self.construction.elements.append(elem)
                        else:
                            pyface.MessageDialog(message='Element can not be added because it already exists', 
                                    severity='warning', title='Adding an Element').open()
                        self.interaction_mode = 'select'
            if isinstance(new, model.Element):
                self._hl_element_actor.visibility = True

    def _construct_handle_widget(self):
        handle = tvtk.HandleWidget(allow_handle_resize=False)
        return handle

    def _make_handle_for_joint(self, joint):
        handle = self._construct_handle_widget()
        representation = tvtk.PolygonalHandleRepresentation3D(handle=self._representation_model_for_joint_handle(joint))
        representation.property.set(ambient=0, diffuse=1, diffuse_color=(0.9, 0.9, 0.9), ambient_color=(1,0,0))
        handle.set_representation(representation)
        handle.interactor = self.scene.interactor
        handle.representation.world_position = (joint.position[0], joint.position[1], 0)
        handle.enabled = True
        # bind movements so that widget movement moves the joint
        def widget_move(*args):
            x,y = handle.representation.world_position[:2]
            joint.x = x
            joint.y = y
        handle.add_observer('InteractionEvent', widget_move)
        # bind selection
        def widget_select(*args):
            self.selected_object = joint
        handle.add_observer('StartInteractionEvent', widget_select)
        self._interaction_mode_changed()
        handle.on()
        return handle

    @on_trait_change('construction.elements, construction:joints:position, _on_init')
    def _setup_elements_picker(self):
        if self.scene and self.construction:
            #print "setup elements picker"
            pd = tvtk.PolyData(points=pts2dto3d(self.construction.joint_positions), lines=self.construction.element_index_table)
            self._pick_elements_actor = tvtk.Actor(mapper=tvtk.PolyDataMapper(input=pd))
            self._element_picker = tvtk.CellPicker(pick_from_list=True, tolerance=0.005)
            self._element_picker.pick_list.append(self._pick_elements_actor)


