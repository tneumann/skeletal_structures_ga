#! /usr/bin/env python

# Copyright 2009 Thomas Neumann
#
# Redistribution of this file is permitted under
# the terms of the GNU Public License (GPL) version 2.

import enthought.pyface.api as pyface
# out of some strange reason we need to load chaco very early (before vtk)
from enthought.chaco.api import Plot, ArrayPlotData

import sys

import data
import ga_truss as GAT
import model
import editor

# initialize the editor
gui = pyface.GUI()
app = editor.ConstructionEditor(
        construction = model.Construction(available_element_materials=data.rajan_materials, element_deleted_material=data.air), 
        ga = GAT.default_genetic_algorithm, 
        new_element_material = data.rajan_materials[0], 
        size = (1000,680), gui = gui )
app.open()
app.scene.render()
app.scene.reset_zoom()

if len(sys.argv) == 1:
    gui.start_event_loop()
