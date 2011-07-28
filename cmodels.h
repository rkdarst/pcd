/* Richard Darst, July 2011 */

typedef struct Graph {
  int n;
  int Ncmty;

  int *cmty;
  int *interactions;

  int **cmtyl;
  int  *cmtyll;
  int *cmtyN;

  } *Graph_t;

double energy(Graph_t G, double gamma);
