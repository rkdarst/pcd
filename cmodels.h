/* Richard Darst, July 2011 */

typedef struct Graph {
  int n;

  int *cmty;
  int *cmtyN;
  int *interactions;

  int **cmtyi;
  int  *cmtyii;

  } *Graph_t;

double energy(Graph_t G, int gamma);
