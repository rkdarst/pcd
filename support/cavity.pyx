
import cython
import math
cimport libc.math
import networkx
import numpy
cimport numpy
import random

#from exceptions import Exception
#cdef extern class NanError(Exception):
#    pass
NanError = ValueError

def norm(x):
    return numpy.sqrt(numpy.sum(numpy.multiply(x, x)))
def shuffled(x):
    x = list(x)
    random.shuffle(x)
    return x

cdef class Cavity(object):
    cdef public int N
    cdef public int q
    cdef public double avg_degree
    cdef public object adj
    cdef public object adj_inv

    cdef public object psi
    cdef public object psi1
    cdef public object h
    cdef public object c
    cdef public object n

    #cdef double[:, :] psi_d
    cdef object psi_d
    cdef double[:, :] psi1_d
    cdef double[:] h_d
    cdef double[:, :] c_d
    cdef double[:] n_d

    cdef public double verbosity
    cdef public object nodeMap
    cdef public object nodeMap_inv

    cdef initial
    cdef int[:] initial_d

    #conv_infer = 10  # multiplied by self.N if conv_infer_intensive==True
    cdef public double conv_infer  # multiplied by self.N if conv_infer_intensive==True
    cdef public int conv_infer_intensive
    cdef public int t_max
    cdef public double conv_learn  # multiplied by self.q if conv_learn_intensive==True
    cdef public int conv_learn_intensive
    cdef public int t_max_learn


    def __init__(self, g, q, **kwargs):
        self.q = q
        N = self.N = len(g)
        #assert set(g.nodes()) == set(range(N))
        self.nodeMap_inv = dict(enumerate(sorted(g.nodes())))
        self.nodeMap = dict((x, i) for (i,x) in self.nodeMap_inv.iteritems())
        self.avg_degree = 2*g.number_of_edges()/float(g.number_of_nodes())

        self.adj = [ numpy.asarray(tuple(self.nodeMap[ngh] for ngh in adj_list), dtype=int)
                     for adj_list in g.adjacency_list() ]
        self.adj_inv = [ dict((j, j_index) for j_index, j in enumerate(self.adj[i]) )
                         for i in range(self.N)]
        self.psi = [ numpy.zeros((len(x), q)) for x in self.adj ]
        self.psi1 = numpy.zeros((N, q))
        self.h = numpy.zeros(q)
        self.c = numpy.zeros((q, q), dtype=numpy.double)
        self.n = numpy.zeros(q)
        self.verbosity = 1

        self.initial = None
        self.init(**kwargs)
        self.t_learn = 0

        #self.psi_d = self.psi
        #self.psi_d = [ double [:,:] x for x in psi ]
        self.psi1_d = self.psi1
        self.h_d = self.h
        self.c_d = self.c
        self.n_d = self.n

        self.conv_infer = 10./200.  # multiplied by self.N if conv_infer_intensive==True
        self.conv_infer_intensive = False
        self.t_max = 100
        self.conv_learn = 0.005  # multiplied by self.q if conv_learn_intensive==True
        self.conv_learn_intensive = True
        self.t_max_learn = 10



    def init(self, epsilon=None, p_in=None, p_out=None, sizes=None, initial=None):
        """Initialize the algorithm.

        p_in: known p_in.  self.c[diagonal] is set to self.N*p_in

        epsilon: set c_out/c_in ratio.  Lower should be better, but
        should not be zero or else it can cause numerical
        instabilities.  Try .01
        """
        if epsilon is not None:
            c_avg = self.avg_degree*self.q  # XXX is *self.q correct?
            c_in = c_avg/(1.+epsilon)
            c_out = c_avg*epsilon/(1.+epsilon)
        else:
            if p_out is not None:
                c_out = p_out*self.N
            else:
                c_out = .01
            if p_in is not None:
                c_in = p_in*self.N
            else:
                c_in = min(10., .5*self.N/(self.q))
        # Do the actual setting of c_matrix.
        self.c[:] = c_out
        for i in range(self.q):
            self.c[i,i] = c_in

        # If we are given an initial state, then set it
        if initial:
            keys = set(initial[i] for i in range(self.N))
            if len(keys) > self.q:
                raise ValueError("Initial community set has more communities than q.")
            cmty_map = dict((x, i) for i,x in enumerate(keys))  # map keys to range(q)
            #assert keys == set(range(self.q)), "Initial community set is not range(0, q)."

            self.initial = numpy.zeros(self.N, dtype=numpy.int32)
            self.initial_d = self.initial
            for i in range(self.N):
                self.initial[i] = cmty_map[initial[self.nodeMap[i]]]
            for q in range(self.q):
                self.n[q] = sum(1 for i in range(self.N) if self.initial[i]==q)/float(self.N) # XXX make more efficient
        elif sizes:
            for q in range(self.q):
                self.n[q] = sizes[q]/float(sum(sizes))
        else:
            self.n[:] = 1./self.q
        if self.verbosity >= 1.5:
            print self.c


    def initialize(self):
        """Create initial psi vectors (initial community assignments).

        Initializes self.psi.  Then initializes self.psi1 and self.h
        from those values.

        If self.initial is not None, then use that to seed initial communities."""
        cdef int[:] initial_d
        # initialize psi
        if self.initial is None:
            for i in range(self.N):
                for neigh_idx in range(len(self.adj[i])):
                    vec = [ ]
                    for q in range(self.q):
                        vec.append(random.uniform(0, 1))
                    norm = float(sum(vec))
                    vec = [ x/norm for x in vec ]
                    self.psi[i][neigh_idx,:] = vec
        else:
            # self.initial is not None.  self.initial set in __init__
            # to a map from nodes to range(q), so should all be
            # properly initialized.
            initial_d = self.initial_d
            for i in range(self.N):
                vec = [0.] * self.q
                vec[initial_d[i]] = 1.
                for neigh_idx in range(len(self.adj[i])):
                    self.psi[i][neigh_idx,:] = vec

        self.update_psi1()
        self.update_h()

    @cython.cdivision(True)
    @cython.boundscheck(False)
    cdef update_h(self):
        cdef double[:,:] c = self.c_d
        cdef double[:,:] psi1 = self.psi1_d
        cdef double[:] h = self.h_d
        cdef double hh
        for q in range(self.q):
            hh = 0.0
            for k in range(self.N):
                for q2 in range(self.q):
                    hh += c[q2, q] * psi1[k,q2]
            h[q] = hh/self.N
    #@cython.cdivision(True)
    #@cython.boundscheck(False)
    cdef update_psi1(self):
        cdef double[:,:] c = self.c_d
        psi = self.psi
        #cdef double[:,:] psi = self.psi_d
        cdef double[:,:] psi1 = self.psi1_d
        cdef double[:] n = self.n_d
        cdef double[:] h
        if True:
            h = self.h_d
        else:
            # normalize h:
            h = self.h.copy()
            avg_ = numpy.mean(h)
            for q in range(self.q):
                h[q] -= avg_
        conv = 0.0
        #c = self.c
        cdef double psi1_
        cdef double sum_
        cdef double norm_
        for k in range(self.N):
            for tk in range(self.q):
                #Zk = self.Zi(k)
                #if Zk == 0:
                #    print "Zk(%d)==0"
                #    raise
                psi1_ = n[tk] * libc.math.exp(-h[tk])
                for j in self.adj[k]:
                    k_index = self.adj_inv[j][k] #list(self.adj[j]).index(k)  # XXX slow
                    #print 'psi[j][k_index]:', self.psi[j][k_index]

                    sum_ = 0
                    for tj in range(self.q):
                        # Using the /N version of this
                        #sum_ += c[tj, tk] * psi[j][k_index, tj]
                        sum_ += c[tj, tk] * psi[j][k_index, tj] #/ self.N
                    #print 'sum:', sum_
                    psi1_ *= sum_
                psi1[k, tk] = psi1_
            # do normalization
            norm_ = numpy.sum(psi1[k, :])
            #print 'n:', self.n
            #print 'c:', self.c
            #print 'h:', self.h
            #print 'psi_new', tuple(psi1[k,:])
            #print 'psi1 norm:', norm_
            #raw_input('psi1>')
            for q in range(self.q):
                psi1[k, q] = psi1[k, q] / norm_
    cdef update_psi1_one(self, int k):
        cdef double[:,:] c = self.c_d
        psi = self.psi
        #cdef double[:,:] psi = self.psi_d
        cdef double[:,:] psi__j
        cdef double[:,:] psi1 = self.psi1_d
        cdef double[:] n = self.n_d
        cdef double[:] h = self.h_d
        adj = self.adj
        adj_inv = self.adj_inv
        cdef double psi1_
        cdef double sum_
        cdef double norm_
        #psi1_new = [ ]
        cdef int j, k_index
        cdef int tk, tj
        cdef int q = self.q
        psi1_new = numpy.zeros(q) # dtype=numpy.double
        for tk in range(q):
            psi1_ = n[tk] * libc.math.exp(-h[tk])
            #print '  psi1_, start:', psi1_
            for j in adj[k]:
                k_index = adj_inv[j][k] #list(self.adj[j]).index(k)  # XXX slow
                psi__j = psi[j]

                sum_ = 0
                for tj in range(q):
                    #sum_ += c[tj, tk] * psi[j][k_index, tj]
                    sum_ += c[tj, tk] * psi__j[k_index, tj]
                    #print '      int:', c[tj, tk] * psi[j][k_index, tj], c[tj, tk], psi[j][k_index, tj]
                #print '    sum_:', sum_
                psi1_ *= sum_
            #print '  psi1_', psi1_
            #psi1[k, tk] = psi1_
            psi1_new[tk] = psi1_
        # Normalize
        #norm_ = numpy.sum(psi1[k, :])
        norm_ = psi1_new.sum()
        #print 'update_psi1_one norm:', norm_
        if norm_ == 0.0:
            print 'h:', self.h
            print 'n:', self.n
            print 'c:', self.c
            print 'norm_:', norm_
            print 'psi_new:', tuple(psi1[k, :])
            raise ValueError("vector norms to zero: exhdfi")
        #raw_input('psi1_one>')
        for q in range(self.q):
            #psi1[k, q] = psi1[k, q] / norm_
            psi1[k, q] = psi1_new[q] / norm_

    cdef double Zij(self, i, j_index, j, i_index) except? -1.:
        cdef double[:,:] c = self.c_d
        psi = self.psi
        #cdef double[:,:] psi = self.psi_d
        #cdef double[:,:] psi1 = self.psi1_d
        cdef double[:] n = self.n_d
        #cdef double[:] h = self.h_d
        cdef double sum_ = 0
        if i_index == -1:
            i_index = self.adj_inv[j][i] #list(self.adj[j]).index(i)  # XXX slow
        for a in range(self.q):
            for b in range(a):   # goes to a-1
                sum_ += c[a, b] * \
                       (psi[i][j_index,a]*psi[j][i_index,b] +
                        psi[i][j_index,b]*psi[j][i_index,a])
        for a in range(self.q):
            sum_ += c[a, a] * psi[i][j_index,a]*psi[j][i_index,a]
        return sum_
    cdef double Zi(self, i) except? -1.:
        cdef double[:,:] c = self.c_d
        psi = self.psi
        #cdef double[:,:] psi = self.psi_d
        #cdef double[:,:] psi1 = self.psi1_d
        cdef double[:] n = self.n_d
        cdef double[:] h = self.h_d
        cdef double sum_ = 0.0
        cdef double innersum
        cdef double prod
        for ti in range(self.q):
            prod = n[ti] * libc.math.exp(-h[ti])
            for j in self.adj[i]:
                i_index = self.adj_inv[j][i] #list(self.adj[j]).index(i)  # XXX slow
                innersum = 0.0
                for tj in range(self.q):
                    innersum += c[tj, ti]*psi[j][i_index, tj] #XXXfixmaybe?
                prod *= innersum
            sum_ += prod
        return sum_
    cpdef extern double c_total(self):
        """Calculate the average degree from self.c and self.n."""
        cdef double[:,:] c = self.c
        cdef double[:] n = self.n
        cdef double c_total = 0
        # undirected version:
        #for a in range(self.q):
        #    for b in range(self.q):
        #        c_total += c[a,b]*n[a]*n[b]
        # directed version (two loops):
        cdef int a, b
        for a in range(self.q):
            for b in range(a+1, self.q):
                c_total += c[a,b]*n[a]*n[b]
        for a in range(self.q):
            c_total += c[a,a]*n[a]*n[a]/2.
        return c_total
    def f(self):
        """Free energy of this partition"""
        cdef double f = 0
        for i in range(self.N):
            f += - math.log(self.Zi(i))/self.N
            # second term
            for j_index, j in enumerate(self.adj[i]):
                if j > i:    # We want to iterate through each edge only once.
                    continue
                # pass i_index = -1, calculate it in the function.
                f += math.log(self.Zij(i, j_index, j, -1))/self.N
        f += - self.c_total()/2.
        #f += - self.c_total()/2.
        #f += self.c_total()/2.
        #f += - self.avg_degree/2.
        print 'c_avg, c_total:', self.avg_degree, self.c_total()
        return f


    #def update_psi(self):
    #    cdef double[:,:] c = self.c_d
    #    psi = self.psi
    #    cdef double[:,:] psi1 = self.psi1_d
    #    cdef double[:] n = self.n_d
    #    cdef double[:] h = self.h_d
    #    order = [ (i, neigh, #q
    #               )
    #              #for q in range(self.q)
    #              for neigh in self.adj[n]
    #              for i in range(self.N) ]
    #    random.shuffle(order)
    #    cdef double conv = 0.0
    #    cdef double sum_
    #    cdef double psi_
    #    raise NotImplementedError("Incomplete.")
    #    for i, j, in order:
    #        psi_old = psi[i][neigh_idx,:].copy()
    #        for ti in range(q):
    #
    #            psi_ = 1.0   # XXX normalizing factor
    #            psi_ *= n[ti] * libc.math.exp(-h[ti])
    #            for k in self.adj[i]:
    #                if k == j: continue
    #                sum_ = 0.0
    #                for tk in range(self.q):
    #                    sum_ += c[tk, ti] * psi[k][i, tk]
    #                psi_ *= sum_
    cdef update_psi_loop(self):
        cdef int q = self.q
        cdef double[:,:] c = self.c_d
        psi = self.psi
        cdef double[:,:] psi1 = self.psi1_d
        cdef double[:] n = self.n_d
        cdef double[:] h = self.h_d
        cdef double conv = 0.0
        cdef double psi_
        cdef int i_index
        cdef double sum_
        cdef double norm_
        adj = self.adj
        adj_inv = self.adj_inv
        cdef int i, j, k

        cdef double [:,:] psi__k

        # This is a shared temporary variable for all rounds, each
        # round replaces the previous contents.
        psi_new = numpy.zeros(q, dtype=numpy.double)

        order = [ (i, neigh_idx)
                  for i in range(self.N)
                  for neigh_idx in range(len(self.adj[i]))
                  ]
        random.shuffle(order)
        for round, (i, neigh_idx) in enumerate(order):
            if round%100000 == 0:
                print round, i, neigh_idx
            psi_old = psi[i][neigh_idx,:].copy()
            j = adj[i][neigh_idx]

            # Update it  (Eq. 1)
            for ti in range(q):
                psi_ = 1.0   # XXX normalizing factoer?
                psi_ *= n[ti] * libc.math.exp(-h[ti])
#                print '  psi start:', psi_
                for k in adj[i]:
                    if k == j: continue
                    i_index = adj_inv[k][i]  # XXX slow
                    psi__k = psi[k]
                    sum_ = 0.0
                    for tk in range(q):
                        #sum_ += c[tk, ti] * psi[k][i_index, tk]
                        sum_ += c[tk, ti] * psi__k[i_index, tk]
#                        print '      int:', c[tk, ti] * psi[k][i_index, tk], c[tk, ti], psi[k][i_index, tk]
#                    print '    sum:', sum_
                    psi_ *= sum_
#                print '  psi:', psi_
                #psi_new.append(psi_)
                psi_new[ti] = psi_
            #print psi_new
            norm_ = psi_new.sum()
            #psi_new = [ x/norm_ for x in psi_new ]
            psi_new /= norm_
#            print "update_psi_loop norm:", norm_, tuple(psi_new)
            if numpy.isnan(psi_new[0]):
                #print self.psi
                print self.h
                print self.n
                print self.c
                print norm_
                print psi_new
                raise ValueError("NaNs detected: hdmbddh")

            psi[i][neigh_idx,:] = psi_new

            conv += norm(psi_old - psi_new)

            self.update_psi1_one(self.adj[i][neigh_idx])
            self.update_h()  # XXX fixme - make this more efficient
        return conv

    cdef update_n(self):
        #print 'n:', self.n
        for a in range(self.q):
            sum = 0.0
            for i in range(self.N):
                sum += self.psi1[i, a]
            self.n[a] = sum/float(self.N)
        #print 'n:', self.n
        if numpy.isnan(self.n).any():
            #print self.psi1
            print self.N
            print self.n
            print self.c
            raise ValueError("NaNs in self.n: %s"%str(self.n).replace('\n', ' '))
    cdef update_c(self):
        #print self.c
        for a in range(self.q):
            for b in range(self.q):
                # equation begins here
                sum = 0.0
                for i in range(self.N):
                    for j_index, j in enumerate(self.adj[i]):
                        # Want to iterate over edges, so only one-way
                        if j > i: continue
                        i_index = self.adj_inv[j][i]  # XXX slow
                        neum = self.c[a, b] * \
                               (self.psi[i][j_index,a]*self.psi[j][i_index,b] +
                                self.psi[i][j_index,b]*self.psi[j][i_index,a])
                        deno =  self.Zij(i,j_index, j, i_index)
                        #print "a,b, i,j, n,d:", a,b,i,j,neum, deno
                        sum += neum/deno
                #print "N, n[a], n[b]:", self.N, self.n[a], self.n[b]
                self.c[a,b] = sum / float(self.N * self.n[a] * self.n[b])
        #print self.c

    def run_inference(self):
        #print self.adj
        #print self.psi
        #print self.psi1
        #print self.detected_communities()
        #raise
        self.initialize()
        if self.verbosity >= 2:
            print self.detected_communities()
        #if self.verbosity >= 1:
        #    print self.n
        #    print self.c
        conv = float('nan')
        for t in range(self.t_max):
            #if self.verbosity >= 1:
            #    print "BP-inteference iteration %d (last conv was %f)"%(t, conv)
            t += 1
            #conv = 0
            #self.update_psi()
            conv = self.update_psi_loop()
            if self.verbosity >= 1:
                print "%3d %6.2f"%(t, conv),
                if self.verbosity >=2:
                    print self.detected_communities()
                    print 'h:', self.h
                else: print
            if self.conv_infer_intensive:
                conv_criteria = self.conv_infer * self.N
            else:
                conv_criteria = self.conv_infer
            if conv < conv_criteria:
                break
            #for x in self.psi:
            #    print x
            #print self.psi
        if self.verbosity >= 1.49:
            print "BP-inteference done (%d iterations, last conv was %f)"%(t, conv)
            print self.f()
    cpdef t_learn
    def run_learn(self):
        conv = float('nan')
        for i in range(self.t_max_learn):
            self.t_learn = i+1
            #if self.verbosity >= 1:
            #    print "BP-learn iteration %d (last conv was %f)"%(i,conv)
            self.run_inference()
            n_old = self.n
            p_old = self.p
            p_cmp = p_old*numpy.outer(n_old, n_old)
            self.update_n()
            self.update_c()
            #print self.psi
            if self.verbosity >= 1.49:
                print 'n:', self.n
                print 'c:', self.c
            p_cmp_new = self.p * numpy.outer(self.n, self.n)
            #print 'scaled p:'
            #print p_cmp
            #print p_cmp_new
            #conv = numpy.sum(numpy.abs(self.n - n_old)) + numpy.sum(numpy.abs(self.p-p_old))
            conv = numpy.sum(numpy.abs(self.cmtysizes()*(self.n - n_old))) \
                   + numpy.sum(numpy.abs(p_cmp-p_cmp_new))
            #print conv, numpy.sum(numpy.abs(self.cmtysizes()*(self.n-n_old))), numpy.sum(numpy.abs(p_cmp-p_cmp_new))
            if self.conv_learn_intensive:
                conv_criteria = self.conv_learn * self.q
            else:
                conv_criteria = self.conv_learn
            if self.verbosity >= 1:
                print "BP-learn done(%d iters, last conv was %f, f=%f)"%(i,conv, self.f())
            print conv, conv_criteria
            print p_cmp
            print p_cmp_new
            print n_old
            print self.n
            print self.q
            #print self.psi1
            if conv < conv_criteria:
                if self.verbosity >= 2:
                    print self.detected_communities()
                break
            #raw_input('(conv: %s)>'%conv);
    cpdef detected_communities(self):
        cmtys = [ ]
        for i in range(self.N):
            cmtys.append(numpy.argmax(self.psi1[i]))
        return cmtys
    def detected_communities_iter(self):
        for i in range(self.N):
            yield numpy.argmax(self.psi1[i])
    def cmtysizes(self):
        sizes = numpy.zeros(self.q)
        for x in self.detected_communities():
            sizes[x] += 1
        return sizes
    cpdef cmtynodes(self):
        cmtynodes = { }
        for i, c in enumerate(self.detected_communities()):
            if c not in cmtynodes:
                cmtynodes[c] = set()
            cmtynodes[c].add(self.nodeMap_inv[i])
        return cmtynodes
    @property
    def p(self):
        return self.c / self.N



# cython cavity.pyx -a --embed=main
# gcc -I/usr/include/python2.6/ -rdynamic -fPIC -pie -lpython2.6 cavity.c -o cavity.so

def main():
    from pcd import sbmgraphs
    import pcd.graphs
    n = 20
    p = .8
    pOut = .2
    q = 2
    n = 20;  q=2;  p=10./n; pOut=5./n
    #n = 25000;  q=4;  p=8.5/n; pOut=2.5/n   # 4-groups test in ε=.3 (detectable)
    #n = 25000;  q=4;  p=7.3/n; pOut=2.9/n   # 4-groups test in ε=.4 (decelle limit)
    #g1 = networkx.complete_graph(n)
    #g2 = networkx.complete_graph(n)
    #g = networkx.disjoint_union(g1, g2)
    #g.add_edge(0, n)
    #g1 = networkx.gnp_random_graph(n, p)
    #g2 = networkx.gnp_random_graph(n, p)
    #g = networkx.disjoint_union(g1, g2)
    #g.add_edge(0, n)
    n=10000; q=4; p = 12.307692307692307/n; pOut = 1.2307692307692308/n

    #g = sbmgraphs.sbm(U=range(n*q), q=q, n=n, p=p, pOut=pOut,
    #                  edges_constructor=sbmgraphs.add_edges_exact_sparse)

    #g = networkx.Graph(); g.add_edges_from(((0,1), (1,2), (2,0)))
    g = pcd.graphs.karate_club() ; q=2
    #g = pcd.graphs.polbooks() ; q=3

    c = Cavity(g, q)
    c.conv_learn = 1e-10
    c.conv_infer = 1e-10
    c.verbosity = 1.5
    #print c.adj[0]
    #print c.psi[0]
    #print c.h[0]
    #print c.c

    c.init(.1)
    c.c[:] = [[8.96, 1.29], [1.29, 7.87],]
    c.n[:] = [ .525, .475 ]
    c.t_max_learn = 100
    #c.init(p_in=p, p_out=pOut,
    #       initial=dict((n, list(g.node[n]['cmtys'])[0]) for n in g.nodes())
    #       )
    #c.c[:] = 1.
    #for i in range(len(c.c)):
    #    c.c[i,i] = 10.
    #c.n[:] = 1./q
    #c.init()
    #c.run_inference()
    c.run_learn()

if __name__ == "__main__":
    main()
