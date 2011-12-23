#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <glib.h>

#include "cmodels.h"
#include "SFMT.h"

#define DEBUGLISTS (0)
#ifndef DEBUG
#define DEBUG (0)
#endif

int test(Graph_t G) {
  /* Print some test debugging information
   */
  printf("nNodes: %d\n", G->N);
  int i;
  // Print the communities
  printf("cmtys: ");
  for (i=0 ; i<G->N ; i++) {
    printf("%d ", G->cmty[i]);
    G->cmty[i] += 10;
  }
  printf("\n");

  // Print the community lists
  for (i=0 ; i<G->N ; i++) {
    //printf("cmtyi: %x %x %x\n", G->cmtyii, G->cmtyi, G->cmtyi[i]);
    //    printf("cmtyi: %x\n", G->cmtyii);

    //printf("cmty %d: %d %d\n", i, G->cmtyi[i][0], G->cmtyi[i][1]);
  }

  // Print the interaction maxtrix
  for (i=0 ; i<G->N ; i++) {
    int j;
    printf("  ");
    for (j=0 ; j<G->N ; j++) {
#ifdef IMATRIX_T_INT
      printf("%2d ", G->imatrix[i*G->N+j]);
#else
      printf("%5.2f ", G->imatrix[i*G->N+j]);
#endif
    }
    printf("\n");
  }

  return (0);
}

void gThreadInit() {
  //assert(g_thread_supported());
  //assert(!g_thread_get_initialized());
  g_thread_init(NULL);
  //assert(g_thread_get_initialized());
}

int isInCmty(Graph_t G, int c, int n) {
  /* assert(G->oneToOne); */
  /* return (G->cmty[n] == c); */
  if (G->oneToOne)
    return (G->cmty[n] == c);
  else {
    int in_cmty1 = g_hash_table_lookup_extended(G->cmtyListHash[c],
    					       GINT_TO_POINTER(n),
    					       //&n,
    					       NULL, NULL);
    /* int in_cmty1 = cmtyListIsInCmty(G, c, n); */
    if (DEBUG) {
      int in_cmty2 = cmtyListIsInCmty(G, c, n);
      if (DEBUGLISTS) printf("        isInCmty %d %d %d %d\n", c, n,
    			     in_cmty1, in_cmty2);
      assert(!in_cmty1 == !in_cmty2);
    }
    return(in_cmty1);
  }
}


/*
 * Functions which deal with the community lists:
 * - G->cmty[n]   (direct community mapping)
 * - G->cmtyl[c][i]
 * - G->cmtyN[c]
 * - G->Ncmty
 */
inline int cmtyListIsInCmty(Graph_t G, int c, int n) {
  /* Return 1 if particle n is in community c.
   *
   * This is primarily useful for systems that have overlaps. If there
   * are no overlaps (G->oneToOne==1), then G->cmty[n]==0 is a faster
   * test.
   */
  int i;
  for (i=0 ; i<G->cmtyN[c] ; i++) {
    if (G->cmtyl[c][i] == n)
      return (1);
  }
  return (0);
}
inline void cmtyListAddOverlap(Graph_t G, int c, int n) {
  /* Add particle n to community c
   *
   * Adapted for systems that have overlaps.
   */
  if (DEBUG)
    assert(!isInCmty(G, c, n));
  int position = G->cmtyN[c];
  G->cmtyl[c][position] = n;
  G->cmtyN[c]++;
  if (DEBUGLISTS)
    printf("cLAO: %2d %2d %2d %2d %2d\n", c, n, position,G->cmtyl[c][position],
			 G->cmty[n]);
  g_hash_table_insert(G->cmtyListHash[c],
		      GINT_TO_POINTER(n),
		      GINT_TO_POINTER(position)
		      );
  /*G->cmty[n] = c;*/
  if (c >= G->Ncmty)
    G->Ncmty = c+1;
}
inline void cmtyListAdd(Graph_t G, int c, int n) {
  /* Add particle n to community c
   */
  cmtyListAddOverlap(G, c, n);
  // The difference for systems without overlaps is we keep G->cmty[n]
  // up to date
  G->cmty[n] = c;
}
inline void cmtyListRemoveOverlap(Graph_t G, int c, int n) {
  /* Remove particle n from community c
   *
   * Adapted for systems that have overlaps
   */
  void *position_p=(void *) -2;
  int found = g_hash_table_lookup_extended(G->cmtyListHash[c],
					   GINT_TO_POINTER(n),
					   NULL, &position_p);
  int position = GPOINTER_TO_INT(position_p);

  if (DEBUGLISTS)
    printf("cLRO: %2d %2d %2d %2d %2d\n", c, n, position,G->cmtyl[c][position],
	   G->cmty[n]);
  assert (found);
  assert (G->cmtyl[c][position] == n);

  /* // Find where it is in the lists */
  /* for (position=0 ; position<G->cmtyN[c] ; position++) { */
  /*   if (G->cmtyl[c][position] == n) */
  /*     break; */
  /* } */
  if (DEBUG && G->cmtyl[c][position] != n)
    printf("****** wrong particle: c=%d n=%d pos=%d G->cmty[c][pos]=%d\n",
	   c, n, position, G->cmtyl[c][position]);
  // Remove the particle
  if (position != G->cmtyN[c]-1) {
    int m = G->cmtyl[c][G->cmtyN[c]-1];
    G->cmtyl[c][position] = m;
    g_hash_table_insert(G->cmtyListHash[c],
			GINT_TO_POINTER(m),
			GINT_TO_POINTER(position));
    g_hash_table_remove(G->cmtyListHash[c],
			GINT_TO_POINTER(n));
  }
  else {
    //no op, just decrement counter
    g_hash_table_remove(G->cmtyListHash[c],
			GINT_TO_POINTER(n));
  }
  G->cmtyN[c]-- ;
  // If we just removed the greatest-numbered community
  if (c == G->Ncmty-1  &&  G->cmtyN[c] == 0 ) {
    // Altar Ncmty too if we just removed the greatest-numbered community
    int i;
    for (i=G->Ncmty-1 ; i>=0 ; i--) {
      if (G->cmtyN[i] == 0)
	G->Ncmty--;
      else
	break;
    }
  }
}
inline void cmtyListRemove(Graph_t G, int c, int n) {
  /* Remove particle n from community c
   */
  cmtyListRemoveOverlap(G, c, n);
  // The difference for systems without overlaps is we keep G->cmty[n]
  // up to date
  G->cmty[n] = NO_CMTY;
}
inline void cmtyMove(Graph_t G, int n, int cOld, int cNew) {
  /* Move node n from cOld to cNew.
   * This function assumes: the * particle's primary community is
   * cOld, and it is not already in * cNew.
   */
  if (DEBUG) assert (G->cmty[n] == cOld);
  // cOld = G->cmty[n];
  cmtyListRemoveOverlap(G, cOld, n);
  cmtyListAddOverlap(G, cNew, n);
  G->cmty[n] = cNew;
}
inline void cmtyMoveSafe(Graph_t G, int n, int cOld, int cNew) {
  /* Move node n from cOld to cNew.

   * Differes from cmtyMove in that this function can handle several
   * corner c
   */
  if (DEBUG) assert (isInCmty(G, cOld, n));
  //if (isInCmty(G, cOld, n))
  cmtyListRemoveOverlap(G, cOld, n);
  if (! isInCmty(G, cNew, n))
    cmtyListAddOverlap(G, cNew, n);
  G->cmty[n] = cNew;
}
inline void cmtyListInit(Graph_t G) {
  /* Initialize the community lists.
   *
   * Initialize all the G->cmtyl[c][0...i] based on G->cmty[n].
   */
  int c, n;
  // This method is always going to initialize to a oneToOne mapping
  // by definition.
  G->oneToOne = 1;
  // Set all lists lengths to zero
  for (c=0 ; c<G->N ; c++) {
    G->cmtyN[c] = 0;
  }
  // Set Ncmty to zero. cmtyListAdd automatically increments this as needed.
  G->Ncmty = 0;
  // Iterate through particles adding them to community lists
  for (n=0 ; n<G->N ; n++) {
    int c = G->cmty[n];
    G->cmty[n] = NO_CMTY; // Without this line, cmtyListAdd fails errorcheck
    cmtyListAdd(G, c, n);
    if (G->cmty[n] >= G->Ncmty)
      printf("****** community %d is greater than Ncmty=%d\n", G->cmty[n], 
	     G->Ncmty);
  }
  // Reset Ncmty to what it should be
  for(c=G->N-1 ; c>=0 ; c--) {
    if (G->cmtyN[c] > 0) {
      G->Ncmty = c+1;
      break;
    }
  }
  hashInit(G);
}
int cmtyListCheck(Graph_t G) {
  /* Check the community lists for consistency.
   *
   * This function will need revisions once one node can be in
   * multiple groups.
   *
   * Returns the numbers of errors found.
   */
  int errors=0;
  int i, j, n, m, cmty;
  // Check Ncmty is indeed the maximum number of communities.
  for (cmty=0; cmty < G->N; cmty++) {
    if ( G->cmtyN[cmty] > 0   &&  cmty >= G->Ncmty ) {
      printf("cmty %d has nodes in it (G->cmtyN[%d]=%d) but Ncmty=%d\n",
	     cmty, cmty, G->cmtyN[cmty], G->Ncmty);
      errors++;
    }
  }
  // First go cmty value -> list
  for (n=0; n < G->N; n++) {
    cmty = G->cmty[n];
    int is_in_list=0;
    for (i=0; i<G->cmtyN[cmty]; i++) {
      if (G->cmtyl[cmty][i] == n) {
	is_in_list = 1;
	break;
      }
    }
    if (! is_in_list) {
      printf("node %d has G->cmty[%d]=%d but not in G->cmtyl[%d][...]\n",
	     n, n, cmty, cmty);
      errors++;
    }
  }
  // Then go and check list -> raw values.  This will break once you
  // can have one node in multiple communities.
  assert(G->oneToOne);
  for (cmty=0; cmty < G->Ncmty; cmty++) {
    for (i=0; i<G->cmtyN[cmty]; i++) {
      n = G->cmtyl[cmty][i];
      if (G->cmty[n] != cmty) {
	printf("node G->cmtyl[%d][%d]=%d not equal to G->cmty[%d]=%d\n",
	       cmty, i, n, n, G->cmty[n]);
	errors++;
      }
    }
  }
  // Check that no node is in a community more than once
  for (cmty=0; cmty < G->Ncmty; cmty++) {
    // For each community
    for (i=0; i<G->cmtyN[cmty]; i++) {
      // For each particle in the community
      n = G->cmtyl[cmty][i];
      for (j=0; j<G->cmtyN[cmty]; j++) {
	// For each particle in that same community, provided they do
	// are not the same:
	if (i == j)
	  continue;
	m = G->cmtyl[cmty][j];
	if (n == m)
	  errors++;
      }
    }
  }
  // check that cmtyN is correct for each list.
  for (cmty=0; cmty < G->Ncmty; cmty++) {
    int n_found = 0;
    for (n=0 ; n<G->N ; n++) {
      if (G->cmty[n] == cmty)
	n_found ++;
    }
    if (n_found != G->cmtyN[cmty])
      errors++;
  }
  return (errors);
}

void hashCreate(Graph_t G) {
  assert(G->cmtyListHash[0] == NULL  &&  G->cmtyListHash[1] == NULL);
  int c;
  for (c=0 ; c<G->N ; c++) {
    G->cmtyListHash[c] = g_hash_table_new(g_direct_hash,
					  g_direct_equal
					  );
  }
}
void hashInit(Graph_t G) {
  /* Initialize hashes from community lists */
  int c, i;
  for (c=0 ; c<G->N ; c++) {
    //assert(g_hash_table_size(G->cmtyListHash[c]) == 0);
    g_hash_table_remove_all(G->cmtyListHash[c]);
    for (i=0 ; i<G->cmtyN[c] ; i++) {
      g_hash_table_insert(G->cmtyListHash[c],
			  GINT_TO_POINTER(G->cmtyl[c][i]),
			  GINT_TO_POINTER(i)
			  );
    }
  }
}
void hashDestroy(Graph_t G) {
  int c;
  for (c=0 ; c<G->N ; c++) {
    g_hash_table_destroy(G->cmtyListHash[c]);
    G->cmtyListHash[c] = NULL;
  }
  //free(G->cmtyListHash);
}
void hashPrintKeys(GHashTable *HT) {
  void printkey(void *key, void *value, void *data) {
    value=NULL;
    data=NULL;
    printf("%d ", GPOINTER_TO_INT(key));
  }
  g_hash_table_foreach(HT, printkey, NULL);
  printf("\n");
}






int n_intersect_nodes(int *cmtyl0, int *cmtyl1,
		      int  cmtyN0, int  cmtyN1) {
  /* Return the number of nodes that are in both the first and second
   * communities (set intersection)
   */
  int n_shared=0;
  int i, j, n, m;
  for (i=0 ; i<cmtyN0 ; i++) {
    // For nodes in the first community
    n = cmtyl0[i];
    for (j=0 ; j<cmtyN1 ; j++) {
      // For nodes in the second community
      m = cmtyl1[j];
      // If the nodes are the same, then increment the number of
      // shared particles and continue.
      if (m == n) {
	n_shared++;
	continue;
      }
    }
  }
  return (n_shared);
}
int n_union_nodes(int *cmtyl0, int *cmtyl1,
		  int  cmtyN0, int  cmtyN1) {
  /* Return the number of nodes that are in the first OR second
   * communities (set union)
   */
  // Start with the number of nodes in the second community.
  int n_union=cmtyN1;
  int i, j, n, m;
  for (i=0 ; i<cmtyN0 ; i++) {
    // For nodes in the first community
    n = cmtyl0[i];
    int is_in_c1=0;
    for (j=0 ; j<cmtyN1 ; j++) {
      // For nodes in the second community
      m = cmtyl1[j];
      // Is node in c0 also in c1?
      if (m == n) {
	is_in_c1 = 1;
	break;
      }
    }
    // If node in c0 is *not* in c1, then increment n_union.
    if (! is_in_c1)
      n_union++;
  }
  return (n_union);
}
int cmtyIsSubset(Graph_t G, int csmall, int cbig) {
  /* is csmall a subset of cbig?  This also allows csmall to be equal
   * to cbig (not a strict superset).  */
  GHashTableIter hashIter;
  GHashTable *hashSmall = G->cmtyListHash[csmall];
  GHashTable *hashBig   = G->cmtyListHash[cbig];
  g_hash_table_iter_init(&hashIter, hashSmall);
  void *n_p;
  while (g_hash_table_iter_next(&hashIter, &n_p, NULL)) {
    //int n = GPOINTER_TO_INT(n_p);
    int found = g_hash_table_lookup_extended(hashBig,
					     &n_p,
					     NULL, NULL);
    if (!found)
      return(FALSE);
  }
  return (TRUE);
}


int find_empty_cmty(Graph_t G) {
  /* Find the lowest numbered empty community
   */
  int c;
  for (c=0 ; c<G->N ; c++) {
    if (G->cmtyN[c] == 0) {
      int n;
      // A bit of error-checking
      for (n=0 ; n<G->N ; n++) {
	if (G->cmty[n] == c) {
	  printf("Inconsistent state: c=%d should be empty (n=%d)\n", c, n);
	  exit(53);
	}
      } // end errorchecking
      return(c);
    }
  }
  return(-1);
}

//int is_in_list()


int q(Graph_t G) {
  /* Number of communities in graph G.
   */
  int q=0;
  int c;
  for (c=0 ; c<G->Ncmty ; c++) {
    if (G->cmtyN[c] > 0)
      q++;
  }
  return (q);
}
double entropy(Graph_t G) {
  /* Information entropy of community allocation of graph G
   * Requires graph to be one-to-one.
   */
  assert(G->oneToOne);
  double N = (double)G->N;
  double H = 0;
  int c;
  for (c=0; c<G->Ncmty ; c++) {
    int n = G->cmtyN[c];
    if (n == 0)
      continue;
    H += n/N * log2((double)(n/N));
  }
  return -H;
}

double mutual_information(Graph_t G0, Graph_t G1) {
  assert(G0->N == G1->N);
  int N = G0->N;
  double MI=0.0;

  int c0, c1, n0, n1, n_shared;
  for (c0=0 ; c0 < G0->Ncmty ; c0++) {
    n0 = G0->cmtyN[c0];
    for (c1=0 ; c1 < G1->Ncmty ; c1++) {
      n1 = G1->cmtyN[c1];

      n_shared = n_intersect_nodes(G0->cmtyl[c0], G1->cmtyl[c1],
				   n0, n1);
      if (n_shared == 0)
	continue;
      MI += (n_shared/(float)N) * log2((double)(n_shared*N/((double)n0*n1)));
    }
  }
  return (MI);
}

inline double h(double p) {
  if ((p==0) || (p==1))
    return 0;
  return -p * log2(p);
}
inline double H(Graph_t G, int c) {
  return (  h((       G->cmtyN[c]) / (double)G->N)
	  + h((G->N - G->cmtyN[c]) / (double)G->N)
         );
}
double H2(Graph_t GX, Graph_t GY, int cX, int cY) {
  double N = (double) GX->N;
  // cX, cY members, cX, cY number in cmty
  int *cXm = GX->cmtyl[cX];
  int *cYm = GY->cmtyl[cY];
  int  cX_n = GX->cmtyN[cX];
  int  cY_n = GX->cmtyN[cY];

  //printf("  c %d %d\n", n_intersect_nodes(cXm, cYm, cX_n, cY_n),
  //                 n_union_nodes(cXm, cYm, cX_n, cY_n));
  double hP11 = h((       n_intersect_nodes(cXm, cYm, cX_n, cY_n))/N);
  double hP10 = h((cX_n - n_intersect_nodes(cXm, cYm, cX_n, cY_n))/N);
  double hP01 = h((cY_n - n_intersect_nodes(cXm, cYm, cX_n, cY_n))/N);
  double hP00 = h(( N   - n_union_nodes    (cXm, cYm, cX_n, cY_n))/N);
  if (hP11 + hP00 <= hP01 + hP10)
    return 1/0.;
  double hPY1 = h( (  cY_n) / N);
  double hPY0 = h( (N-cY_n) / N);
  return(hP11+hP00+hP01+hP10 - hPY1 - hPY0);
}
double HX_Ynorm(Graph_t GX, Graph_t GY) {
  double HX_Y_total=0;  // These two are to find the average
  int HX_Y_n=0;
  int cX;
  for (cX=0 ; cX < GX->Ncmty; cX++) {
    if (GX->cmtyN[cX] == 0)
      continue;
    double HX_Yhere=1/0.;
    int cY;
    for (cY=0 ; cY < GY->Ncmty ; cY++) {
      if (GY->cmtyN[cY] == 0)
	continue;
      double H2_this = H2(GX, GY, cX, cY);
      if (H2_this < HX_Yhere)
	HX_Yhere = H2_this;
    }
    if (HX_Yhere == 1/0.)
      HX_Yhere = H(GX, cX);
    double _ = H(GX, cX);
    if (_ == 0) {
      HX_Y_total += 0;
      HX_Y_n += 1;
    }
    else {
      HX_Y_total += HX_Yhere / _;
      HX_Y_n += 1;
    }

  }
  return(HX_Y_total / HX_Y_n);
}





double energy_naive(Graph_t G, double gamma) {
  /* Naive energy loop, looping over all pairs of particles.  SLOW.
   */
  assert(0);
  assert(G->hasFull);
  imatrix_t attractions=0;
  imatrix_t repulsions =0;
  int n;
  for (n=0 ; n<G->N ; n++) {
    int m;
    for (m=0 ; m<G->N ; m++) {
      if (m == n)
	continue;
      if (G->cmty[m] == G->cmty[n]) {
	if (G->rmatrix == NULL) {
	  imatrix_t interaction = G->imatrix[n*G->N + m];
	  if (interaction > 0)
	    repulsions  += interaction;
	  else
	    attractions += interaction;
	}
	else {
	  attractions += G->imatrix[n*G->N + m];
	  repulsions += G->rmatrix[n*G->N + m];
	}
      }
    }
  }
  return(.5 * (attractions + gamma*repulsions));
}



double energy(Graph_t G, double gamma) {
  /* Calculate energy using community lists.  Much faster.
   */
  if (!G->hasFull) return energy_sparse(G, gamma);
  assert(G->hasFull);
  imatrix_t attractions=0;
  imatrix_t repulsions =0;
  int c, n;
  //cmtyListInit(G);

  for (c=0 ; c<G->Ncmty ; c++) {
    // for communities c
    int i, j, m;
    for (i=0 ; i<G->cmtyN[c] ; i++) {
      // Do symmetric: both directions.
      for (j=0 ; j<G->cmtyN[c] ; j++) {
	if (i == j)
	  continue;
	n = G->cmtyl[c][i];
	m = G->cmtyl[c][j];
	assert(n != m);
	if (G->rmatrix == NULL) {
	  imatrix_t interaction = G->imatrix[n*G->N + m];
	  if (interaction > 0)
	    repulsions  += interaction;
	  else
	    attractions += interaction;
	}
	else {
	  attractions += G->imatrix[n*G->N + m];
	  repulsions += G->rmatrix[n*G->N + m];
	}
      }
    }
  }
  return(.5 * (attractions + gamma*repulsions));
}


double energy_sparse(Graph_t G, double gamma) {
  /* Calculate energy using community lists.  Much faster.
   */
  assert(0); // Not implemented yet.
  assert(G->hasSparse);
  imatrix_t attractions=0;
  imatrix_t repulsions =0;
  int c, n;

  int nDefined=0;
  int j;
  for (n=0 ; n<G->N ; n++) {
    int m = G->simatrixId[n*G->simatrixLen + j];
    if (G->cmty[m] != c)
      continue;
    nDefined += 1;
  }

  for (c=0 ; c<G->Ncmty ; c++) {
    // for communities c
    int i, j, m;
    for (i=0 ; i<G->cmtyN[c] ; i++) {
      // Do symmetric: both directions.
      for (j=0 ; j<G->cmtyN[c] ; j++) {
	if (i == j)
	  continue;
	n = G->cmtyl[c][i];
	m = G->cmtyl[c][j];
	assert(n != m);
	if (G->rmatrix == NULL) {
	  imatrix_t interaction = G->imatrix[n*G->N + m];
	  if (interaction > 0)
	    repulsions  += interaction;
	  else
	    attractions += interaction;
	}
	else {
	  attractions += G->imatrix[n*G->N + m];
	  repulsions += G->rmatrix[n*G->N + m];
	}
      }
    }
  }
  return(.5 * (attractions + gamma*repulsions));
}




double energy_cmty(Graph_t G, double gamma, int c) {
  /* Calculate the energy of only one community `c`.
   */
  assert(G->hasFull);
  imatrix_t attractions=0;
  imatrix_t repulsions =0;
  int n;

  // for communities c
  int i, j, m;
  for (i=0 ; i<G->cmtyN[c] ; i++) {
    // Do symmetric: both directions.  Someday, this could matter for
    // directed graphs, right now it is irrelevent.
    for (j=0 ; j<G->cmtyN[c] ; j++) {
  	if (i == j)
  	  continue;
  	n = G->cmtyl[c][i];
  	m = G->cmtyl[c][j];
	assert(n != m);
	if (G->rmatrix == NULL) {
	  imatrix_t interaction = G->imatrix[n*G->N + m];
	  if (interaction > 0)
	    repulsions  += interaction;
	  else
	    attractions += interaction;
	}
	else {
	  attractions += G->imatrix[n*G->N + m];
	  repulsions += G->rmatrix[n*G->N + m];
	}
    }
  }
  return(.5 * (attractions + gamma*repulsions));
}

double energy_cmty_n(Graph_t G, double gamma, int c, int n) {
  /* Calculate the energy of only one community `c`, if it had node n
   * in it.  Node n does not have to actually be in that community.
   */
  assert(G->hasFull);
  imatrix_t attractions=0;
  imatrix_t repulsions =0;

  // for communities c
  int j, m;
  for (j=0 ; j<G->cmtyN[c] ; j++) {
    m = G->cmtyl[c][j];
    if (m == n)
      continue;
    if (G->rmatrix == NULL) {
      imatrix_t interaction = G->imatrix[n*G->N + m];
      if (interaction > 0)
	repulsions  += interaction;
      else
	attractions += interaction;
    }
    else {
      attractions += G->imatrix[n*G->N + m];
      repulsions += G->rmatrix[n*G->N + m];
    }
  }
  return(.5 * (attractions + gamma*repulsions));
}

double energy_cmty_n_sparse(Graph_t G, double gamma, int c, int n) {
  /* Calculate the energy of only one community `c`, if it had node n
   * in it.  Node n does not have to actually be in that community.
   */
  assert(G->hasSparse);
  imatrix_t attractions=0;
  imatrix_t repulsions =0;

  int j;
  int nUnDefined = G->cmtyN[c];
  // For each adjoining particle in the list:
  for (j=0 ; j<G->simatrixN[n] ; j++) {
    //int m = G->simatrixId[n*G->simatrixLen + j];
    int m = G->simatrixIdl[n][j];
    if (m == n)
      continue;
    if (! isInCmty(G, c, m))
      continue;
    /* if (G->cmty[m] != c) */
    /*   continue; */
    nUnDefined -= 1;

    if (G->srmatrix == NULL) {
      imatrix_t interaction = G->simatrix[n*G->simatrixLen + j];
      if (interaction > 0)
	repulsions  += interaction;
      else
	attractions += interaction;
    }
    else {
      attractions += G->simatrix[n*G->simatrixLen + j];
      repulsions += G->srmatrix[n*G->simatrixLen + j];
    }
  }
  // -1 comes from self interaction not being counted as undefined.
  if (isInCmty(G, c, n))
    nUnDefined -= 1;
  attractions += nUnDefined * G->simatrixDefault;
  repulsions  += nUnDefined * G->srmatrixDefault;

  double E = .5 * (attractions + gamma*repulsions);
  /* printf(" e_c_n_s %d %d(%d,%d)\n", c, n, G->cmtyN[c], nUnDefined); */
  //if (DEBUG && E != energy_cmty_n(G, gamma, c, n)) {
  if (DEBUG) {
    double E2 = energy_cmty_n(G, gamma, c, n);
    if ( fabs(E-E2)>.01  && fabs(E-E2)/E > .0001) {
      printf(" e_c_n_s %d %d(%d,%d) (%d) %f %f\n", c, n, G->cmtyN[c],nUnDefined,
	     isInCmty(G, c, n),
	     E, E2);
      //assert(!cmtyListCheck(G));
      assert(0);
    }
  }
  return(E);
}
double energy_cmty_n_which(Graph_t G, double gamma, int c, int n) {
  if (G->hasSparse) return energy_cmty_n_sparse(G, gamma, c, n);
  return energy_cmty_n(G, gamma, c, n);
}



double energy_cmty_cmty(Graph_t G, double gamma, int c1, int c2) {
  /* Total energy of interaction between two communities.
   */
  assert(G->hasFull);
  double E=0;
  int i1, n1;
  for (i1=0 ; i1 < G->cmtyN[c1] ; i1++) {
    n1 = G->cmtyl[c1][i1];
    E += energy_cmty_n(G, gamma, c2, n1);
  }
  return (E);
}
double energy_n(Graph_t G, double gamma, int n) {
  /* Energy of particle n in its own community.
   */
  assert(G->hasFull);
  assert(G->oneToOne);
  int c = G->cmty[n];
  return energy_cmty_n(G, gamma, c, n);
}


int minimize_naive(Graph_t G, double gamma) {
  /* OBSELETE minimization routine.
   */
  assert(0);
  assert(G->hasFull);
  int changes=0;
  int i, n;
  // Loop over particles
  for (i=0 ; i<G->N ; i++) {
    n = G->randomOrder[i];
    int oldcmty  = G->cmty[n];
    int bestcmty = G->cmty[n];
    double Ebest = energy(G, gamma);
    /* printf("Partile %d, old community %d\n", i, oldcmty); */
    cmtyListRemove(G, oldcmty, n);


    int newcmty;
    for (newcmty=0 ; newcmty<G->N ; newcmty++) {

      // Try partiicle in each new cmty
      if (newcmty == oldcmty)
	continue;
      if (G->cmtyN[newcmty] == 0) {
	continue;
      }

      double Enew;
      //G->cmty[n] = newcmty;
      cmtyListAdd(G, newcmty, n);
      Enew = energy(G, gamma);
      cmtyListRemove(G, newcmty, n);

      if (Enew < Ebest) {
      /* printf("  Better option for particle %d: %d %d %d\n", */
      /* 	    i, oldcmty, bestcmty, newcmty); */
      bestcmty = newcmty;
      Ebest = Enew;
      }
    }
    //G->cmty[n] = bestcmty;
    cmtyListAdd(G, bestcmty, n);
    if (oldcmty != bestcmty) {
      // Change community
      changes += 1;
      /* printf("particle %4d: cmty change %4d->%4d\n",  */
      /*        i, oldcmty, bestcmty); */
      //G->cmtyN[oldcmty]  --;
      //G->cmtyN[bestcmty] ++;
    }
  }
return (changes);
}

int minimize(Graph_t G, double gamma) {
  /* Core minimization routine.  Do one sweep, moving each particle
   * (in order of G->randomOrder into the community that most lowers the
   * energy.
   */
  if (G->hasSparse) return(minimize_sparse(G, gamma));
  assert(G->hasFull);
  int changes=0;
  int nindex, n;
  // Loop over particles
  for (nindex=0 ; nindex<G->N ; nindex++) {
    n = G->randomOrder[nindex];
    // Keep a record of the running best community to move to.
    // Default to current community (no moving).
    double deltaEbest = 0.0;
    int bestcmty = G->cmty[n];

    // Store our old community and energy change when we remove a
    // particle from the old community.  We see if (energy from
    // removing from old community + energy from adding to new
    // community) is less than deltaEbest to see where we should move.
    int oldcmty  = G->cmty[n];
    double deltaEoldCmty = - energy_cmty_n(G, gamma, oldcmty, n);


    // Try particle in each new cmty.  Accept the new community
    // that has the lowest new energy.
    // There are various ways of doing this inner loop:

    // Method 1 (all other communities) //
    /* int newcmty; */
    /* for (newcmty=0 ; newcmty<G->Ncmty ; newcmty++) { */

    // Method 2 (only interacting cmtys, fixed order) //
    int m;
    for (m=0 ; m<G->N ; m++) {
      /* if (G->imatrix[G->N*n+m] != 500) */
      /* 	printf("  %d %d %f\n", n, m, G->imatrix[G->N*n+m]); */

      if (G->imatrix[n*G->N + m] > 0)
    	continue;
      int newcmty = G->cmty[m];

    // Method 3 (only interacting cmtys, random order) //
    /* int mindex, m; */
    /* for (mindex=0 ; mindex<G->N ; mindex++) { */
    /*   m = G->randomOrder2[mindex]; */
    /*   if (G->imatrix[n*G->N + m] > 0) */
    /* 	continue; */
    /*   int newcmty = G->cmty[m]; */

      if (newcmty == oldcmty)
	continue;
      if (G->cmtyN[newcmty] == 0) {
	continue;
      }

      double deltaEnewCmty = energy_cmty_n(G, gamma, newcmty, n);

      // Our conditional on if we want to move to this new place.  If
      // we do, update our bestcmty and deltaEbest to say so.
      if (deltaEoldCmty + deltaEnewCmty < deltaEbest) {
	bestcmty = newcmty;
	deltaEbest = deltaEoldCmty + deltaEnewCmty;
      }
    }
    // Is it better to move a particle into an _empty_ community?
    if (deltaEoldCmty < deltaEbest) {
      bestcmty = find_empty_cmty(G);
      // deltaEbest = deltaEoldCmty;  // Not needed (not used after this)
      assert(bestcmty != -1);
    }
    if (oldcmty != bestcmty) {
      cmtyMove(G, n, oldcmty, bestcmty);
      changes += 1;
    }
  }
  return (changes);
}


int minimize_sparse(Graph_t G, double gamma) {
  /* Core minimization routine.  Do one sweep, moving each particle
   * (in order of G->randomOrder into the community that most lowers the
   * energy.
   */
  assert(G->hasSparse);
  int changes=0;
  int nindex, n;

  /* double E_avg = energy(G, gamma) / G->N; */

  // Loop over particles
  for (nindex=0 ; nindex<G->N ; nindex++) {
    n = G->randomOrder[nindex];

    /* float E_particle = energy_n(G, gamma, n); */
    /* if (E_particle < E_avg) */
    /*   continue; */

    // Keep a record of the running best community to move to.
    // Default to current community (no moving).
    double deltaEbest = 0.0;
    int bestcmty = G->cmty[n];

    // Store our old community and energy change when we remove a
    // particle from the old community.  We see if (energy from
    // removing from old community + energy from adding to new
    // community) is less than deltaEbest to see where we should move.
    int oldcmty  = G->cmty[n];
    double deltaEoldCmty = - energy_cmty_n_sparse(G, gamma, oldcmty, n);

    // Try particle in each new cmty.  Accept the new community
    // that has the lowest new energy.
    // There are various ways of doing this inner loop:

    SetClear(G->seenList);
    int j;
    for(j=0 ; j<G->simatrixN[n] ; j++) {
      if (G->simatrix[G->simatrixLen*n + j] >= 0)
      	continue;
      int m = G->simatrixId[G->simatrixLen*n + j];
      int newcmty = G->cmty[m];


      if (newcmty == oldcmty)
	continue;
      if (G->cmtyN[newcmty] == 0) {
	continue;
      }
      if (SetContains(G->seenList, newcmty))
      	continue;
      SetAdd(G->seenList, newcmty);

      double deltaEnewCmty = energy_cmty_n_sparse(G, gamma, newcmty, n);

      // Our conditional on if we want to move to this new place.  If
      // we do, update our bestcmty and deltaEbest to say so.
      if (deltaEoldCmty + deltaEnewCmty < deltaEbest) {
	bestcmty = newcmty;
	deltaEbest = deltaEoldCmty + deltaEnewCmty;
      }
    }
    // Is it better to move a particle into an _empty_ community?
    if (deltaEoldCmty < deltaEbest) {
      bestcmty = find_empty_cmty(G);
      // deltaEbest = deltaEoldCmty;  // Not needed (not used after this)
      assert(bestcmty != -1);
    }
    if (oldcmty != bestcmty) {
      cmtyMove(G, n, oldcmty, bestcmty);
      changes += 1;
    }
  }
  return (changes);
}


int overlapMinimize_add(Graph_t G, double gamma) {
  /* Do a minimization attempt, but adding particles to new
   * overlapping communities.  However, we can't remove particles from
   * their original communities.
   *
   * One principle in this function is that G->cmty[n] lists the
   * _original_ community each particle was in.
   *
   * This makes a round of adding particles to communities.
   */
  int changes = 0;
  G->oneToOne = 0;
  assert(G->hasSparse);

  /* int nindex, n; */
  /* // Loop over particles */
  /* for (nindex=0 ; nindex<G->N ; nindex++) { */
  /*   n = G->randomOrder[nindex]; */

    /* // Method 1 (all other communities) // */
    /* int c; */
    /* for (c=0 ; c<G->Ncmty ; c++) { */
    /*   if (G->cmtyN[c] == 0) */
    /* 	continue; */

    /* // Method 2 (all other communities, random order) // */
    /* int cindex, c; */
    /* for (cindex=0 ; cindex<G->N ; cindex++) { */
    /*   c = G->randomOrder2[cindex]; */
    /*   if (c >= G->Ncmty) */
    /* 	continue; */
    /*   if (G->cmtyN[c] == 0) */
    /* 	continue; */

    /* // Method 3 (all neighboring communities using sparse matrix) // */
    /* LListClear(G->seenList); */
    /* int mindex; */
    /* for (mindex=0 ; mindex<G->simatrixN[mindex] ; mindex++) { */
    /*   int m = G->simatrixId[mindex]; */
    /*   int c = G->cmty[m]; */
    /*   if (LListContains(G->seenList, c)) */
    /* 	continue; */
    /*   LListAdd(G->seenList, c); */

    /* // Method 4 (all neighboring cmtys using sparse matrix, correctly) // */
    /* LListClear(G->seenList); */
    /* int mindex; */
    /* for (mindex=0 ; mindex<G->simatrixN[mindex] ; mindex++) { */
    /*   int m = G->simatrixId[mindex]; */
    /*   int c; */
    /* for (c=0 ; c<G->Ncmty ; c++) { */
    /*   if (! isInCmty(G, c, m)) */
    /* 	continue; */
    /*   if (LListContains(G->seenList, c)) */
    /* 	continue; */
    /*   LListAdd(G->seenList, c); */


  // Method Beta - loop over communities first.
  int c;
  for (c=0 ; c<G->Ncmty ; c++) {
    if (G->cmtyN[c] == 0)
      continue;
    SetClear(G->seenList);
    int mindex;
    for (mindex=0 ; mindex<G->cmtyN[c] ; mindex++) {
      int m;
      m = G->cmtyl[c][mindex];
    int nindex;
    for (nindex=0 ; nindex<G->simatrixN[m] ; nindex++) {
      if (G->simatrix[m*G->simatrixLen + nindex] >= 0)
	continue;
      int n;
      n = G->simatrixIdl[m][nindex];
      if (SetContains(G->seenList, n))
      	continue;
      SetAdd(G->seenList, c);



      // If this is in the original comminty, can skip it.
      if (c == G->cmty[n])
	continue;
      // Skip if already in this community.
      if (
	  isInCmty(G, c, n)
	  //cmtyListIsInCmty(G, c, n)
	  //g_hash_table_lookup_extended(G->cmtyListHash[c], (void *)(gint64)n,NULL, NULL)
	  )
	continue;

      // Should this be added to community?
      double deltaE = energy_cmty_n_which(G, gamma, c, n);
      if (deltaE < 0) {
	if (DEBUGLISTS) printf("    Adding in overlapmin_add %d %d\n", c, n);
	cmtyListAddOverlap(G, c, n);
	changes++;
      }

    }  // Extra brace
    }
  }
  return (changes);
}
int overlapMinimize_remove(Graph_t G, double gamma) {
  /* Do a minimization attempt, but adding particles to new
   * overlapping communities.  However, we can't remove particles from
   * their original communities.
   *
   * One principle in this function is that G->cmty[n] lists the
   * _original_ community each particle was in.
   *
   * This makes a round of removing particles from each community.
   */
  int changes = 0;
  G->oneToOne = 0;

  int c;
  // Loop over communities
  for (c=0 ; c<G->Ncmty ; c++) {

    // Loop over particles within that community
    int i, n;
    for (i=0 ; i<G->cmtyN[c] ; i++) {
      n = G->cmtyl[c][i];

      // If this is in the original comminty, can skip the tests.
      if (c == G->cmty[n])
	continue;

      // Should this be removed from the community?
      double deltaE = energy_cmty_n_which(G, gamma, c, n);
      if (deltaE > 0) {
	if (DEBUGLISTS) printf("    Removing in overlapmin_rem %d %d\n", c, n);
	cmtyListRemoveOverlap(G, c, n);
	changes++;
      }
    }
  }
  return (changes);
}

int anneal(Graph_t G, double gamma, double beta,
	   int steps, double deltabeta) {
  int changes=0;
  int step;

  for (step=0 ; step<steps ; step++) {
    beta += deltabeta;


    // Random node
    int n = G->N * genrand_real2();
    // Random community
    int c = G->N * genrand_real2();
    int c_old = G->cmty[n];
    // Random chance of trying a new communty
    if (genrand_real2() < .0001 ) {
      c = find_empty_cmty(G);
      if (c == -1)
	continue;
    }

    // If already in this community, skip it.
    if (c == c_old)
      continue;
    float deltaE =   energy_cmty_n(G, gamma, c, n)
                   - energy_cmty_n(G, gamma, c_old, n);

    double x = exp(- beta * deltaE);
    double ran = genrand_real2();
    if (ran < x) {
      // accept
      cmtyMove(G, n, c_old, c);
      changes += 1;
    }
    else {
      // reject
    }
  }
  return (changes);
}



int combine_cmtys(Graph_t G, double gamma) {
  /* Attempt to merge communities to get a lower energy assignment.
   * Pairwise attempt to merge all.
   */
  if (G->hasSparse) return(combine_cmtys_sparse(G, gamma));
  assert(G->oneToOne);
  assert(G->hasFull);
  int changes = 0;
  // Move particles from c2 into c1
  int i1, i2, c1, c2;
  //printf("gamma: %f\n", gamma);
  /* for (c1=0 ; c1<G->Ncmty-1 ; c1++) { */
  /*   if (G->cmtyN[c1] == 0) */
  /*     continue; */
  /*   int bestcmty = c1; */
  /*   double deltaEbest = 0; */
  /*   for (c2=c1+1 ; c2<G->Ncmty ; c2++) { */
  /*     if (G->cmtyN[c2] == 0) */
  /* 	continue; */
  for (i1=0 ; i1<G->N ; i1++) {
    c1 = G->randomOrder[i1];
    if (G->cmtyN[c1] == 0)
      continue;

    int bestcmty = c1;
    double deltaEbest = 0;

    for (i2=0 ; i2<G->N ; i2++) {
      c2 = G->randomOrder2[i2];
      if (c1 <= c2)
  	continue;
      if (G->cmtyN[c2] == 0)
  	continue;
      /* //double Eold = energy(G, gamma); */
      /* double Eold = energy_cmty(G, gamma, c1) + energy_cmty(G, gamma, c2);*/


      // Calculate change of energy if we moved c1 into c2
      double deltaE = 0;
      int j1, n1;
      for (j1=0 ; j1 < G->cmtyN[c1] ; j1++) {
      	n1 = G->cmtyl[c1][j1];
      	deltaE += energy_cmty_n(G, gamma, c2, n1);
      }
      for (j1=0 ; j1 < G->cmtyN[c2] ; j1++) {
      	n1 = G->cmtyl[c2][j1];
      	deltaE += energy_cmty_n(G, gamma, c1, n1);
      }

      // Do we accept this change?
      if (deltaE < deltaEbest) {
	bestcmty = c2;
	deltaEbest = deltaE;
      }
    }

    // No change
    if (c1 == bestcmty)
      continue;

    int c2oldN = G->cmtyN[bestcmty];
    //Move all from bestcmty into c1
    int i;
    for (i=c2oldN-1 ; i>=0 ; i--) {
      cmtyMove(G, G->cmtyl[bestcmty][i], bestcmty, c1);
    }

    changes += 1;


  }
  return (changes);
}


int combine_cmtys_sparse(Graph_t G, double gamma) {
  /* Attempt to merge communities to get a lower energy assignment.
   * Pairwise attempt to merge all.
   */
  assert(G->oneToOne);
  assert(G->hasSparse);
  int changes = 0;
  // Move particles from c2 into c1
  //int i1, i2, c1, c2;
  int c1;
  //printf("gamma: %f\n", gamma);
  for (c1=0 ; c1<G->Ncmty-1 ; c1++) {
    if (G->cmtyN[c1] == 0)
      continue;
    int bestcmty = c1;
    double deltaEbest = 0;


  /* for (i1=0 ; i1<G->N ; i1++) { */
  /*   c1 = G->randomOrder[i1]; */
  /*   if (G->cmtyN[c1] == 0) */
  /*     continue; */

  /*   int bestcmty = c1; */
  /*   double deltaEbest = 0; */

  /*   for (i2=0 ; i2<G->N ; i2++) { */
  /*     c2 = G->randomOrder2[i2]; */
  /*     if (c1 <= c2) */
  /*   	continue; */
  /*     if (G->cmtyN[c2] == 0) */
  /*   	continue; */
      /* //double Eold = energy(G, gamma); */
      /* double Eold = energy_cmty(G, gamma, c1) + energy_cmty(G, gamma, c2);*/


    /* for (c2=c1+1 ; c2<G->Ncmty ; c2++) { */
    /*   if (G->cmtyN[c2] == 0) */
    /* 	continue; */


    SetClear(G->seenList);
    int ii;
    for (ii=0 ; ii<G->cmtyN[c1] ; ii++ ) {
      int n = G->cmtyl[c1][ii];
      //printf("n %d %d %d\n", c1, ii, n);
    int j;
    for(j=0 ; j<G->simatrixN[n] ; j++) {
      if (G->simatrix[G->simatrixLen*n + j] >= 0)
      	continue;
      int m = G->simatrixId[G->simatrixLen*n + j];
      int c2 = G->cmty[m];
      //printf("  m %d %d %d %d\n", c1, j, m, c2);
      if (c1==c2 || SetContains(G->seenList, c2)) {
    	//printf("already in: %d\n", c2);
      	continue;
      }
      SetAdd(G->seenList, c2);
      //printf("adding: %d\n", c2);



      // Calculate change of energy if we moved c1 into c2
      double deltaE = 0;
      int j1, n1;
      for (j1=0 ; j1 < G->cmtyN[c1] ; j1++) {
      	n1 = G->cmtyl[c1][j1];
      	deltaE += energy_cmty_n_sparse(G, gamma, c2, n1);
      }
      for (j1=0 ; j1 < G->cmtyN[c2] ; j1++) {
      	n1 = G->cmtyl[c2][j1];
      	deltaE += energy_cmty_n_sparse(G, gamma, c1, n1);
      }
      //printf("  e %d %d %f %d %f\n", c1, c2, deltaE, bestcmty, deltaEbest);
      // Do we accept this change?
      if (deltaE < deltaEbest) {
	//printf("better one\n");
	bestcmty = c2;
	deltaEbest = deltaE;
      }
    }
    }
    //printf("loop done\n");

    // No change
    if (c1 == bestcmty)
      continue;
    //printf("  c %d %d combining\n", c1, c2);

    int c2oldN = G->cmtyN[bestcmty];
    //Move all from bestcmty into c1
    int i;
    for (i=c2oldN-1 ; i>=0 ; i--) {
      cmtyMove(G, G->cmtyl[bestcmty][i], bestcmty, c1);
    }

    changes += 1;

  }
  return (changes);
}



int remap_cmtys(Graph_t G) {
  /* Moves all the communities to the lowest numbered continuous
     segments.
   */
  int changes=0;
  int c;
  for (c=G->Ncmty-1 ; c>=0 ; c--) {
    // For each community, starting with the highest.
    if (G->cmtyN[c] == 0)
      // Ignore already-empty communities
      continue;
    int cNew = find_empty_cmty(G);
    if (cNew >= c)
      break;
    if (cNew == -1)
      break;

    // Do the actual remapping:
    int i;
    for(i=G->cmtyN[c]-1 ; i>=0 ; i--) {
      cmtyMove(G, G->cmtyl[c][i], c, cNew);
    }
    changes++;
  }
  return(changes);
}


