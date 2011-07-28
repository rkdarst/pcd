#include <stdio.h>

#include "cmodels.h"

int test(Graph_t G) {
  printf("nNodes: %d\n", G->n);
  int i;
  // Print the communities
  for (i=0 ; i<G->n ; i++) {
    printf("%d ", G->cmty[i]);
    G->cmty[i] += 10;
  }
  printf("\n");

  // Print the community lists
  for (i=0 ; i<G->n ; i++) {
    //printf("cmtyi: %x %x %x\n", G->cmtyii, G->cmtyi, G->cmtyi[i]);
    //    printf("cmtyi: %x\n", G->cmtyii);

    //printf("cmty %d: %d %d\n", i, G->cmtyi[i][0], G->cmtyi[i][1]);
  }

  return (0);
}


void cmtyListAdd(Graph_t G, int c, int n) {
  // Add particle n to community c
  G->cmtyl[c][G->cmtyN[c]] = n;
  G->cmtyN[c]++;
  G->cmty[n] = c;
}
void cmtyListRemove(Graph_t G, int c, int n) {
  // Remove particle n from community c
  int i;
  // Find where it is in the lists
  for (i=0 ; i<G->cmtyN[c] ; i++) {
    if (G->cmtyl[c][i] == n)
      break;
  }
  if (G->cmtyl[c][i] != n)
    printf("****** wrong particle: %d %d %d %d\n", c, n, i, G->cmtyl[c][i]);
  // Remove the particle
  if (i != G->cmtyN[c]-1) {
    G->cmtyl[c][i] = G->cmtyl[c][G->cmtyN[c]-1];
  }
  else {
    //no op, just decrement counter
  }
  G->cmtyN[c]-- ;
  G->cmty[n] = -1;
}
void cmtyListInit(Graph_t G) {
  /* Initialize the community lists.
   */
  int c, n;
  // Set all lists lengths to zero
  for (c=0 ; c<G->n ; c++) {
    G->cmtyN[c] = 0;
  }
  // Iterate through particles adding them to community lists
  for (n=0 ; n<G->n ; n++) {
    cmtyListAdd(G, G->cmty[n], n);
    if (G->cmty[n] >= G->Ncmty)
      printf("****** community %d is greater than Ncmty=%d\n", G->cmty[n], 
	     G->Ncmty);
  }
  // Reset Ncmty to what it should be
  for(c=G->n-1 ; c>=0 ; c--) {
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
   */
  int errors=0;
  int i, n, cmty;
  // Check Ncmty is indeed the maximum number of communities.
  for (cmty=0; cmty < G->n; cmty++) {
    if ( G->cmtyN[cmty] > 0   &&  cmty >= G->Ncmty ) {
      printf("cmty %d has nodes in it (G->cmtyN[%d]=%d) but Ncmty=%d\n",
	     cmty, cmty, G->cmtyN[cmty], G->Ncmty);
      errors++;
    }
  }
  // First go cmty value -> list
  for (n=0; n < G->n; n++) {
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
  
  
  return (errors);
}





int minimize(Graph_t G, double gamma) {

  int changes=0;
  int n;
  // Loop over particles
  for (n=0 ; n<G->n ; n++) {
    int oldcmty  = G->cmty[n];
    int bestcmty = G->cmty[n];
    double Ebest = energy(G, gamma);
    /* printf("Partile %d, old community %d\n", i, oldcmty); */
    cmtyListRemove(G, oldcmty, n);

    int newcmty;
    for (newcmty=0 ; newcmty<G->n ; newcmty++) {
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
	/* 	 i, oldcmty, bestcmty, newcmty); */
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

double energy1(Graph_t G, double gamma) {
  // Naive energy loop, looping over all pairs of particles.
  int attractions=0;
  int repulsions =0;
  int n;
  for (n=0 ; n<G->n ; n++) {
    int m;
    for (m=0 ; m<G->n ; m++) {
      if (m == n) 
	continue;
      if (G->cmty[m] == G->cmty[n]) {
	int interaction = G->interactions[n*G->n + m];
	if (interaction > 0)
	  repulsions  += interaction;
	else
	  attractions += interaction;
      }
    }
  }
  return(.5 * (attractions + gamma*repulsions));
}



double energy(Graph_t G, double gamma) {
  // Calculate energy using community lists.
  int attractions=0;
  int repulsions =0;
  int c, n;
  //cmtyListInit(G);

  for (c=0 ; c<G->n ; c++) {
    // for communities c
    int i, j, m;
    for (i=0 ; i<G->cmtyN[c] ; i++) {
      // Do symmetric: both directions.
      for (j=0 ; j<G->cmtyN[c] ; j++) {
	if (i == j) 
	  continue;
	n = G->cmtyl[c][i];
	m = G->cmtyl[c][j];
	int interaction = G->interactions[n*G->n + m];
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
  /* Calculate the energy of only one community.
   */
  int attractions=0;
  int repulsions =0;
  int n;

  // for communities c
  int i, j, m;
  for (i=0 ; i<G->cmtyN[c] ; i++) {
    // Do symmetric: both directions.
    for (j=0 ; j<G->cmtyN[c] ; j++) {
  	if (i == j) 
  	  continue;
  	n = G->cmtyl[c][i];
  	m = G->cmtyl[c][j];
  	int interaction = G->interactions[n*G->n + m];
  	if (interaction > 0)
  	  repulsions  += interaction;
  	else
  	  attractions += interaction;
    }
  }
  return(.5 * (attractions + gamma*repulsions));
}

int combine_cmtys(Graph_t G, double gamma) {
  int count = 0;
  // Move particles from c2 into c1
  int c1, c2;
  //printf("gamma: %f\n", gamma);
  for (c1=0 ; c1<G->Ncmty-1 ; c1++) {
    for (c2=c1+1 ; c2<G->Ncmty ; c2++) {
      if (G->cmtyN[c1] == 0 || G->cmtyN[c2] == 0)
	continue;
      double Eold = energy(G, gamma);
      
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
      double Enew = energy(G, gamma);
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

