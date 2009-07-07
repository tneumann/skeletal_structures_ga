from bbfreeze import Freezer
import shutil

distdir = "dist"
freeze = Freezer(distdir=distdir,
                 includes=(
                     "enthought.tvtk.pyface.ui.wx.decorated_scene",
                     "enthought.tvtk.tvtk",
                     "vtk",
                     ))

freeze.addScript("app.py", gui_only=True)
freeze()

shutil.copytree('./examples', distdir + '/examples')

