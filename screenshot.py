# little helper module i used to make the screenshots & figures for the paper
from enthought.traits.api import *

class ScreenshotEvolution(HasTraits):
    app = Any()
    def __init__(self, app):
        HasTraits.__init__(self)
        self.best_fitness = 999999999999
        self.app = app

    @on_trait_change('app.ga.on_step')
    def take_screenshot(self):
        if self.app.ga.current_population.best.raw_fitness < self.best_fitness:
            self.app.construction = self.app.ga.world.construction_from_individual(self.app.ga.current_population.best)
            screenshot(self.app, magnification=1, path='/tmp/gen_%05d.png' % self.app.ga.num_steps)
            self.best_fitness = self.app.ga.current_population.best.raw_fitness


def screenshot_all(app, path):
    for i,individual in enumerate(app.selectable_individuals):
        app.selected_individual = individual
        app.gui.process_events()
        screenshot(app, path+'/'+str(i) + '.png')
    app.construction = app.ga.world.construction
    app.gui.process_events()
    screenshot(app, path+'/start.png')

def screenshot(app, path='/tmp/sc.png', magnification=2):
    #app.bw_mode()
    app._text_actor.text_property.font_size *= magnification
    app._label_actor.mapper.label_text_property.font_size = 14*magnification
    app.scene.magnification = magnification
    app.scene.save(path)
    app._label_actor.mapper.label_text_property.font_size = 14
    app._text_actor.text_property.font_size /= magnification


def _drawfig(con, radius_amplify=50, size_factor=0.5, path=None):
    import pylab as P
    import numpy as N

    P.clf()
    for e in con.elements:
        if e.material != con.element_deleted_material:
            P.plot([e.joint1.x, e.joint2.x], [e.joint1.y, e.joint2.y], color='black', linewidth=e.material.radius*radius_amplify)

    for j in con.joints:
        if j.movable_x and j.movable_y:
            marker = 'o'
        elif j.movable_x and not j.movable_y:
            if j.y > 0:
                marker = 'v'
            else:
                marker = '^'
        elif not j.movable_x and j.movable_y:
            if j.x > 0:
                marker = '<'
            else:
                marker = '>'
        else:
            marker = 's'
        P.plot([j.x],[j.y],marker=marker,c='white')

    P.setp(P.gca(), 'yticklabels', [])
    P.setp(P.gca(), 'xticklabels', [])
    P.setp(P.gca(), 'xticks', [])
    P.setp(P.gca(), 'yticks', [])
    ax = P.axis('scaled')
    #xx, yy = abs(ax[0])+abs(ax[1]), abs(ax[2])+abs(ax[3])
    #P.axis([ax[0]-xx/10., ax[1]+xx/10., ax[2]-yy/10., ax[3]+yy/10.])
    P.axis([-con.width/2., con.width/2., -con.height/2., con.height/2.])
    if path:
        P.savefig(path)
    else:
        P.show()

def _setup_figure(con, size_factor):
    import pylab as P
    P.figure(edgecolor='white', facecolor='white', figsize=(con.width*size_factor, con.height*size_factor))
    P.axes(frameon=False)

def figure(con, radius_amplify=50, path=None, size_factor=0.5):
    _setup_figure(con, size_factor)
    _drawfig(con, radius_amplify=radius_amplify, path=path, size_factor=size_factor)
            
def figure_all(app, radius_amplify=50, path='/tmp', format='eps', size_factor=0.5):
    _setup_figure(app.ga.world.construction, size_factor)
    for i,individual in enumerate(app.ga.current_population.individuals):
        con = app.ga.world.construction_from_individual(individual)

        _drawfig(con, radius_amplify=radius_amplify, path = '%s/%d.%s' % (path, i, format), size_factor=size_factor)
    _drawfig(app.ga.world.construction, radius_amplify=radius_amplify, path = '%s/start.%s' % (path, format), size_factor=size_factor)




