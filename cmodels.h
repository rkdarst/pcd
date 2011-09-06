/* Richard Darst, July 2011 */

typedef int imatrix_t;

typedef struct Graph {
  int N;
  int Ncmty;
  int oneToOne;

  int *cmty;
  imatrix_t *imatrix;

  int **cmtyl;
  int  *cmtyll;
  int *cmtyN;
  int *randomOrder;
  int *randomOrder2;

  } *Graph_t;

double energy(Graph_t G, double gamma);
double energy_cmty(Graph_t G, double gamma, int c);
double energy_cmty_n(Graph_t G, double gamma, int c, int n);
