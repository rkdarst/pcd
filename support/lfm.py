"""Inefficient implementation of LFM CD method.


This module is a inefficient and incomplete implementation of the LFM
CD method.  It was made for the sole purpose of making a movie
animating community detection.  It ignores things such as removing
nodes from communities.  When executed as a script, this module makes
a sequence of imaged representing a community detection process.

Reference:

    Detecting the overlapping and hierarchical community structure in
    complex networks Andrea Lancichinetti et al 2009 New J. Phys. 11
    033015
"""



class LFM(object):
    alpha = 1

    def __init__(self, g, alpha=1):
        self.g = g
        self.alpha = alpha
        self.unassigned = set(g.nodes())


    def score(self, nodes):
        k_in = 2*self.g.subgraph(nodes).number_of_edges()
        k_total = sum(self.g.degree(n) for n in nodes)

        return k_in / float((k_total)**self.alpha)

    def seed_set(self, node):
        self.seed = node
        self.cmty = set((node, ))
        self.to_test = set()
        self.tested = set()
        self.unassigned.remove(node)

        self.to_test.update(g.neighbors(node))
        self.cur_score = self.score(self.cmty)


    def test(self):
        if not self.to_test:
            return None
        while True:
            n = random.sample(self.to_test, 1)[0]
            # This is a hack to make a good visualizaiton.  It assumes
            # we know the communities already, and makes sure the
            # first nodes we add match the true community.
            if len(self.cmty) < 3 and n//6 != next(iter(self.cmty))//6:
                continue
            break

        self.to_test.remove(n)

        new_score = self.score(self.cmty | set((n,)))

        if new_score > self.score(self.cmty):
            print '    adding: score %s > %s'%(new_score, self.score(self.cmty))
            self.cmty.add(n)
            self.unassigned.discard(n)
            for neigh in self.g.neighbors(n):
                if neigh in self.cmty:
                    continue
                self.to_test.add(neigh)
            self.tested = set()
        else:
            print '    not adding: score %s <= %s'%(new_score, self.score(self.cmty))
            self.tested.add(n)

        return n


if __name__ == "__main__":
    # This was converted to a movie using:
    # ffmpeg -r 4 -f image2 -i path/to/%05d.png -vcodec mjpeg -r 4
    # -vcodec mpeg4 -y -f avi LFM2.avi

    import networkx as nx
    #g = nx.read_gml('/home/darstr1/pymod/pcd/karate.gml')
    #print len(g)
    import pcd.sbmgraphs
    g, c = pcd.sbmgraphs.sbm2(pin=.9, pout=.1, q=3, n=6)

    import random, numpy, scipy
    random.seed(13); numpy.random.seed(14); scipy.random.seed(15)
    runner = LFM(g, alpha=1.5,)


    import matplotlib ; matplotlib.use('Agg')
    import pylab
    pos = nx.spring_layout(g)
    def savefig(col):
        global idx
        pylab.clf()
        #from pcd.support.matplotlibutil import get_axes, save_axes
        #ax, extra = get_axes('/home/darstr1/proj/tmp-img/%05d.png'%0)
        #print ax
        nx.draw_networkx(g, pos=pos,
                         node_list=g.nodes(),
                         node_color=[col[n] for n in g.nodes()],
                         with_labels=False)
        #save_axes(ax, extra)

        fig = pylab.gcf()
        ax = pylab.gca()
        dpi = fig.dpi
        bbox = ax.collections[0].get_window_extent(fig)
        from matplotlib.transforms import TransformedBbox, Affine2D
        bbox2 = TransformedBbox(bbox, Affine2D().scale(1. / dpi))
        bbox._points[0] -= .15 * dpi
        bbox._points[1] += .15 * dpi

        ax = pylab.gca()
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        pylab.savefig('/home/darstr1/proj/tmp-imgs/%05d.png'%idx,
                      bbox_inches=bbox2)
        idx += 1
    idx = 0
    cid = 0
    cmty_colors = ['blue', 'green', 'purple', 'brown', 'cyan', 'pink']
    cmtys = { }

    colors = dict((n, 'red') for n in g)

    savefig(colors)

    runner.seed_set(0)
    colors[0] = cmty_colors[len(cmtys)]

    for n in runner.to_test:
        colors[n] = 'orange'
    savefig(colors)


    savefig(colors)
    print colors


    while True:
        restart = False
        n = runner.test()
        if n is None:
            if len(runner.unassigned) == 0:
                break
            else:
                restart = True
                cmtys[len(cmtys)] = runner.cmty
                new_seed = random.sample(runner.unassigned, 1)[0]
                runner.seed_set(new_seed)
                #continue
        print n

        if not restart:
            oldc = colors[n]
            colors[n] = 'black'
            savefig(colors)
        #colors[n] = oldc
        #savefig(colors)

        # set community colors

        colors = dict((n, 'red') for n in g)

        for cid in cmtys:
            for n in cmtys[cid]:
                colors[n] = cmty_colors[cid]
        for n in runner.cmty:
            colors[n] = cmty_colors[len(cmtys)]
        for n in runner.to_test:
            colors[n] = 'orange'
        savefig(colors)



    print runner.cmty, len(runner.cmty)
