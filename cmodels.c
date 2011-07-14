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


void addToCmty(Graph_t G, int c, int n) {
  // Add particle n to community c
  G->cmtyi[c][G->cmtyN[c]] = n;
  G->cmtyN[c]++;
  G->cmty[n] = c;
}
void removeFromCmty(Graph_t G, int c, int n) {
  // Remove particle n from community c
  int i;
  // Find where it is in the lists
  for (i=0 ; i<G->cmtyN[c] ; i++) {
    if (G->cmtyi[c][i] == n)
      break;
  }
  if (G->cmtyi[c][i] != n)
    printf("****** wrong particle: %d %d %d %d\n", c, n, i, G->cmtyi[c][i]);
  // Remove the particle
  if (i != G->cmtyN[c]-1) {
    G->cmtyi[c][i] = G->cmtyi[c][G->cmtyN[c]-1];
  }
  else {
    //no op, just decrement counter
  }
  G->cmtyN[c]-- ;
  G->cmty[n] = -1;
}
void initCmty(Graph_t G) {
  int c, n;
  // Set all lists lengths to zero
  for (c=0 ; c<G->n ; c++) {
    G->cmtyN[c] = 0;
  }
  // Iterate through particles adding them to community lists
  for (n=0 ; n<G->n ; n++) {
    add_to_cmty_list(G, G->cmty[n], n);
  }
}




int minimize(Graph_t G, int gamma) {

  int changes=0;
  int i;
  // Loop over particles
  for (i=0 ; i<G->n ; i++) {
    int oldcmty  = G->cmty[i];
    int bestcmty = G->cmty[i];
    double Ebest = energy(G, gamma);
    /* printf("Partile %d, old community %d\n", i, oldcmty); */

    int newcmty;
    for (newcmty=0 ; newcmty<G->n ; newcmty++) {
      // Try partiicle in each new cmty
      if (newcmty == oldcmty)
	continue;
      if (G->cmtyN[newcmty] == 0) {
	continue;
      }
      double Enew;
      G->cmty[i] = newcmty;
      Enew = energy(G, gamma);
      if (Enew < Ebest) {
	/* printf("  Better option for particle %d: %d %d %d\n", */
	/* 	 i, oldcmty, bestcmty, newcmty); */
	bestcmty = newcmty;
	Ebest = Enew;
      }
    }
    G->cmty[i] = bestcmty;
    if (oldcmty != bestcmty) {
      // Change community
      changes += 1;
      /* printf("particle %4d: cmty change %4d->%4d\n",  */
      /*        i, oldcmty, bestcmty); */
      G->cmtyN[oldcmty]  --;
      G->cmtyN[bestcmty] ++;
    }
  }
return (changes);
}

double energy1(Graph_t G, int gamma) {
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



double energy(Graph_t G, int gamma) {
  // Calculate energy using community lists.
  int attractions=0;
  int repulsions =0;
  int c, n;
  make_cmty_list(G);

  for (c=0 ; c<G->n ; c++) {
    // for communities c
    int i, j, m;
    for (i=0 ; i<G->cmtyN[c] ; i++) {
      // Do symmetric: both directions.
      for (j=0 ; j<G->cmtyN[c] ; j++) {
	if (i == j) 
	  continue;
	n = G->cmtyi[c][i];
	m = G->cmtyi[c][j];
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

void combine_communities(Graph_t G, int gamma, int c1, int c2) {
  double Eold = energy(G, gamma);
  // Move particles from c2 into c1
  int i;
  for (i=0 ; i<G->cmtyN[c2] ; i++)
    {
      //Move all of these to community 1
      Eold=0;
    }

}

