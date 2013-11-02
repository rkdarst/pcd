# Richard Darst, September 2013

import pcd.cmty

import networkx


def layout(g):
    #pos = networkx.pydot_layout(g, prog='fdp')
    pos = networkx.pygraphviz_layout(g, prog='neato')
    return pos

def draw(g, fname, cmtys, figsize=(20,20), pos=None):

    hulls = False
    try:
        colors = cmtys.nodecolors()
        colors = [colors.get(n, (0.,0.,0.,1.)) for n in g.nodes()]
    except pcd.cmty.OverlapError:
        colors = 'blue'
        hulls = True


    if pos is None:
        pos = layout(g)



    if any(any(isinstance(v, unicode) for v in ndata.itervalues())
           for n,ndata in g.nodes_iter(data=True)):
        print 'some unicode found'

    import matplotlib.figure
    import matplotlib.backends.backend_agg
    fig = matplotlib.figure.Figure(figsize=figsize)
    canvas = matplotlib.backends.backend_agg.FigureCanvasAgg(fig)
    ax = fig.add_subplot(111, aspect='equal')
    dpiScale = 1



    networkx.draw(g, ax=ax, pos=pos, node_color=colors,
                  #labels = dict((n, 'a') for n in g.nodes()),
                  with_labels=True
                  )

    from matplotlib.patches import Rectangle, Polygon
    import numpy

    if hulls:
        cmtyColormap = cmtys.cmtycolors()

        from support.convexhull import convex_hull
        for c in cmtys.cmtynames():
            ns = cmtys.cmtycontents(c)
            if len(ns) < 1: continue
            points = [ pos[n] for n in ns ]
            points = numpy.asarray(points)
            #if len(points) > 5:
            #    print points
            #    points = convex_hull(points.T)
            points = convex_hull(tuple(p) for p in points)
            # Detect if any points are too long
            #if boxsize is not None:
            if 0:
                pts1 = points[:]
                pts2 = points[1:] + points[0:1]
                d = numpy.subtract(pts1, pts2)
                # We don't sum along axis=-1, instead we do it per-axis
                d /= boxsize # inplace
                #print d
                if numpy.abs(d).max() > .5:
                    continue
            # Do actual plotting
            p = Polygon(points, alpha=.25, color=cmtyColormap[c],
                        zorder=-2,
                        axes=ax)
            ax.add_patch(p)

    canvas.print_figure(fname, dpi=fig.get_dpi()*dpiScale,
                        bbox_inches='tight')
