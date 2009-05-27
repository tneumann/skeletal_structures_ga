Designing Skeletal Structures with Genetic Algorithms
=====================================================

This tool was programmed for an assignement at [my university](http://www.htw-dresden.de) for the course "Genetic Algorithms". It is a GUI program where one can construct 2D skeletal structures which then can be optimized by a genetic algorithm. It contains a full visualization of the stability of the structure and an interactive editor. The genetic algorithm can also be used without GUI from the command line and can be controlled via Python scripts as seen in example_crane.py.

Requirements
------------

To run the program, you need:

* Python >= 2.5
* NumPy
* Traits and TraitsGUI (the wxpython backend)
* VTK >= 5.2 with Python bindings, and TVTK (comes with MayaVi)

For the screenshot module (not needed by the GUI), matplotlib is required.

Running the program
-------------------

Start "app.py" and load an example from the example folder. See the [project page](http://www.htw-dresden.de/~s52567/ga/) for more details plus the paper I wrote that explains the algorithms that were used.

Licence
-------

The code is licenced under the GNU Public License (GPL) version 2.
