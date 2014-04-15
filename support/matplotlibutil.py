# Richard Darst, May 2012

import os
import re

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


def get_axes(fname, figsize=(13, 10), ax_hook=None,
             ret_fig=False, figargs={}):
    """Interface to matplotlib plotting.

    fname: str:
        fname is the string to save to.  It can also have an extension
        of something like FILE.[pdf,png,ps] and it will save a file of
        ALL of these extensions.
    """
    if isinstance(fname, str):
        import matplotlib.figure
        import matplotlib.backends.backend_agg
        fig = matplotlib.figure.Figure(figsize=figsize, dpi=100, **figargs)
        canvas = matplotlib.backends.backend_agg.FigureCanvasAgg(fig)
        if not ret_fig:
            ax = fig.add_subplot(111)#, aspect='equal')
        else:
            ax = fig
        return ax, (fname, canvas, fig, ax, ax_hook)
    if isinstance(fname, matplotlib.axes.Axes):
        ax = fname
        return ax, (fname, ax_hook)
    else:
        raise ValueError("Unknown type of object: %s"%type(fname))
def save_axes(ax, extra):
    fname = extra[0]
    if isinstance(fname, str):
        fname, canvas, fig, ax, ax_hook = extra
        if not os.path.exists(os.path.dirname(fname)):
            os.makedirs(os.path.dirname(fname))
        multi_ext_match = re.match(r'(.*\.)\[([A-Za-z,]+?)\]$', fname)
        if ax_hook: ax_hook(ax, locals())
        if multi_ext_match:
            # Support for multi-extension matching.  If fname matches
            # FILE.[ext1,ext2], then write the figure to ALL of these
            # files.
            base = multi_ext_match.group(1)
            for ext in multi_ext_match.group(2).split(','):
                canvas.print_figure(base+ext, dpi=fig.get_dpi(),
                                    bbox_inches='tight')
        else:
            canvas.print_figure(fname, dpi=fig.get_dpi(),
                                bbox_inches='tight')
    elif isinstance(fname, matplotlib.axes.Axes):
        fname, ax_hook = extra
        if ax_hook: ax_hook(ax, locals())
    else:
        raise ValueError("Unknown type of object: %s (also, how did we even get to this point?)"%type(fname))


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

