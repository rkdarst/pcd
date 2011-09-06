#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "cmodels.h"
#include "SFMT.h"

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
      printf("%2d ", G->imatrix[i*G->N+j]);
    }
    printf("\n");
  }

  return (0);
}


void cmtyListAdd(Graph_t G, int c, int n) {
  /* Add particle n to community c
   */
  G->cmtyl[c][G->cmtyN[c]] = n;
  G->cmtyN[c]++;
  G->cmty[n] = c;
  if (c >= G->Ncmty)
    G->Ncmty = c+1;
}
void cmtyListRemove(Graph_t G, int c, int n) {
  /* Remove particle n from community c
   */
  int i;
  // Find where it is in the lists
  for (i=0 ; i<G->cmtyN[c] ; i++) {
    if (G->cmtyl[c][i] == n)
      break;
  }
  if (G->cmtyl[c][i] != n)
    printf("****** wrong particle: c=%d n=%d i=%d G->cmty[c][i]=%d\n", c, n, i, G->cmtyl[c][i]);
  // Remove the particle
  if (i != G->cmtyN[c]-1) {
    G->cmtyl[c][i] = G->cmtyl[c][G->cmtyN[c]-1];
  }
  else {
    //no op, just decrement counter
  }
  G->cmtyN[c]-- ;
  G->cmty[n] = -1;
  // If we just removed the greatest-numbered community
  if (c == G->Ncmty-1  &&  G->cmtyN[c] == 0 ) {
    // Altar Ncmty to If we just removed the greatest-numbered community
    for (i=G->Ncmty-1 ; i>=0 ; i--) {
      if (G->cmtyN[i] == 0)
	G->Ncmty--;
      else
	break;
    }
  }
}
void cmtyListInit(Graph_t G) {
  /* Initialize the community lists.
   *
   * Set G->cmtyl[c][0...i] based on G->cmty
   */
  int c, n;
  // Set all lists lengths to zero
  for (c=0 ; c<G->N ; c++) {
    G->cmtyN[c] = 0;
  }
  // Set Ncmty to zero. cmtyListAdd automatically increments this as needed.
  G->Ncmty = 0;
  // Iterate through particles adding them to community lists
  for (n=0 ; n<G->N ; n++) {
    cmtyListAdd(G, G->cmty[n], n);
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
  return (errors);
}




int shared_nodes_between_communities(int *cmtyl0, int *cmtyl1,
				     int  cmtyN0, int  cmtyN1) {
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
    H += n/N * log2(n/N);
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

      n_shared = shared_nodes_between_communities(G0->cmtyl[c0], G1->cmtyl[c1],
						  n0, n1);
      if (n_shared == 0)
	continue;
      MI += (n_shared/(float)N) * log2(n_shared*N/((double)n0*n1));
    }
  }
  return (MI);
}



double energy(Graph_t G, double gamma) {
  /* Calculate energy using community lists.  Much faster.
   */
  int attractions=0;
  int repulsions =0;
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
	imatrix_t interaction = G->imatrix[n*G->N + m];
	if (interaction > 0)
	  repulsions  += interaction;
	else
	  attractions += interaction;
      }
    }
  }
  return(.5 * (attractions + gamma*repulsions));
}

double energy_cmty(Graph_t G, double gamma, int c) {
  /* Calculate the energy of only one community `c`.
   */
  int attractions=0;
  int repulsions =0;
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
  	int interaction = G->imatrix[n*G->N + m];
  	if (interaction > 0)
  	  repulsions  += interaction;
  	else
  	  attractions += interaction;
    }
  }
  return(.5 * (attractions + gamma*repulsions));
}

double energy_cmty_n(Graph_t G, double gamma, int c, int n) {
  /* Calculate the energy of only one community `c`, if it had node n
   * in it.  Node n does not have to actually be in that community.
   */
  int attractions=0;
  int repulsions =0;

  // for communities c
  int j, m;
  for (j=0 ; j<G->cmtyN[c] ; j++) {
    m = G->cmtyl[c][j];
    if (m == n)
      continue;
    int interaction = G->imatrix[n*G->N + m];
    if (interaction > 0)
      repulsions  += interaction;
    else
      attractions += interaction;
  }
  return(.5 * (attractions + gamma*repulsions));
}

double energy_cmty_cmty(Graph_t G, double gamma, int c1, int c2) {
  double E=0;
  int i1, n1;
  for (i1=0 ; i1 < G->cmtyN[c1] ; i1++) {
    n1 = G->cmtyl[c1][i1];
    E += energy_cmty_n(G, gamma, c2, n1);
  }
  return (E);
}


int minimize_naive(Graph_t G, double gamma) {
  /* OBSELETE minimization routine.
   */
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
    int newcmty;
    for (newcmty=0 ; newcmty<G->Ncmty ; newcmty++) {

    // Method 2 (only interacting cmtys, fixed order) //
    /* int m; */
    /* for (m=0 ; m<G->N ; m++) { */
    /*   if (G->imatrix[n*G->N + m] > 0) */
    /* 	continue; */
    /*   int newcmty = G->cmty[m]; */

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
      cmtyListRemove(G, oldcmty, n);
      cmtyListAdd(G, bestcmty, n);
      changes += 1;
    }
  }
  return (changes);
}

double energy_naive(Graph_t G, double gamma) {
  /* Naive energy loop, looping over all pairs of particles.  SLOW.
   */
  int attractions=0;
  int repulsions =0;
  int n;
  for (n=0 ; n<G->N ; n++) {
    int m;
    for (m=0 ; m<G->N ; m++) {
      if (m == n)
	continue;
      if (G->cmty[m] == G->cmty[n]) {
	imatrix_t interaction = G->imatrix[n*G->N + m];
	if (interaction > 0)
	  repulsions  += interaction;
	else
	  attractions += interaction;
      }
    }
  }
  return(.5 * (attractions + gamma*repulsions));
}


int combine_cmtys(Graph_t G, double gamma) {
  /* Attempt to merge communities to get a lower energy assignment.
   * Pairwise attempt to merge all.
   */
  int count = 0;
  // Move particles from c2 into c1
  int i1, i2, c1, c2;
  //printf("gamma: %f\n", gamma);
  /* for (c1=0 ; c1<G->Ncmty-1 ; c1++) { */
  /*   for (c2=c1+1 ; c2<G->Ncmty ; c2++) { */
  for (i1=0 ; i1<G->N ; i1++) {
    c1 = G->randomOrder[i1];
    if (G->cmtyN[c1] == 0)
      continue;
    for (i2=0 ; i2<G->N ; i2++) {
      c2 = G->randomOrder2[i2];
      if (c1 <= c2)
	continue;
      if (G->cmtyN[c2] == 0)
	continue;
      //double Eold = energy(G, gamma);
      double Eold = energy_cmty(G, gamma, c1) + energy_cmty(G, gamma, c2);

      int c1oldN = G->cmtyN[c1];
      int c2oldN = G->cmtyN[c2];
      //Move all from c2 into c1
      int i;
      for (i=0 ; i<c2oldN ; i++) {
	G->cmtyl[c1][i+c1oldN] = G->cmtyl[c2][i];
	G->cmty[G->cmtyl[c2][i]] = c1;
      }
      // FIX - breaks for multiple community
      G->cmtyN[c1] = c1oldN + c2oldN;
      G->cmtyN[c2] = 0;

      // Do we accept?
      //double Enew = energy(G, gamma);
      double Enew = energy_cmty(G, gamma, c1);
      if (Enew <= Eold) {
	count++;
	continue;
      }

      // Put it all back
      for (i=0 ; i<c2oldN ; i++) {
	G->cmtyl[c2][i] = G->cmtyl[c1][i+c1oldN];
	G->cmty[G->cmtyl[c2][i]] = c2;
      }
      // FIX - breaks for multiple community
      G->cmtyN[c1] = c1oldN;
      G->cmtyN[c2] = c2oldN;

    }
  }
  return (count);
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
    for(i=0 ; i<G->cmtyN[c] ; i++) {
      G->cmtyl[cNew][i] = G->cmtyl[c][i];
    }
    int n;
    for (n=0 ; n<G->N ; n++) {
      if (G->cmty[n] == c)
	G->cmty[n] = cNew;
    }
    // Change cmtyN and Ncmty
    G->cmtyN[cNew] = G->cmtyN[c];
    G->cmtyN[c] = 0;
    changes++;
    G->Ncmty--;
  }
  // Fix Ncmty to indicate our new max cmty number.
  int i;
  for (i=G->Ncmty-1 ; i>=0 ; i--) {
    if (G->cmtyN[i] == 0)
      G->Ncmty--;
    else
      break;
  }
  return(changes);
}


