/* Richard Darst, July 2011 */

typedef struct Graph {
  int N;
  int Ncmty;
  int oneToOne;

  int *cmty;
  int *interactions;

  int **cmtyl;
  int  *cmtyll;
  int *cmtyN;
  int *randomOrder;
  int *randomOrder2;

  } *Graph_t;

double energy(Graph_t G, double gamma);
double energy_cmty(Graph_t G, double gamma, int c);
