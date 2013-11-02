# Richard Darst, May 2012


import matplotlib.figure
import matplotlib.backends.backend_agg
#from matplotlib.patches import CirclePolygon as Circle
import matplotlib.cm as cm
import matplotlib.colors as colors


def get_fig():
    fig = matplotlib.figure.Figure()
    #ax = fig.add_subplot(111, aspect='equal')
    return fig
def write_fig(fig):
    canvas = matplotlib.backends.backend_agg.FigureCanvasAgg(fig)
    #ax.autoscale_view(tight=True)
    canvas.print_figure(fname, dpi=fig.get_dpi()*dpiScale*escale,
                        bbox_inches='tight')


def get_axes(fname, figsize=(13, 10)):
    if fname:
        import matplotlib.figure
        import matplotlib.backends.backend_agg
        fig = matplotlib.figure.Figure(figsize=figsize)
        canvas = matplotlib.backends.backend_agg.FigureCanvasAgg(fig)
        ax = fig.add_subplot(111)#, aspect='equal')
        return canvas, fig, ax
    else:
        raise
def save_fig(canvas, fig, ax, fname):
    canvas.print_figure(fname, dpi=fig.get_dpi(),
                        bbox_inches='tight')
class Figure(object):
    def __init__(self, fname, figsize=None):
        self.fname = fname
        if fname:
            import matplotlib.figure
            import matplotlib.backends.backend_agg
            self.fig = matplotlib.figure.Figure(figsize=figsize)
            self.canvas = matplotlib.backends.backend_agg.FigureCanvasAgg(self.fig)
            self.ax = self.fig.add_subplot(111)#, aspect='equal')
        else:
            raise
    def save(self):
        self.canvas.print_figure(self.fname, dpi=self.fig.get_dpi(),
                                 bbox_inches='tight')

