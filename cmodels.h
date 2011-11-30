/* Richard Darst, July 2011 */

#include <glib.h>

#include "imatrix_t.h"

typedef IMATRIX_T imatrix_t;

#define NO_CMTY (-1)

typedef struct LList {
  int count;
  int maxcount;
  int *data;
} *LList_t;

typedef struct Graph {
  int N;
  int Ncmty;
  int oneToOne;

  int *cmty;
  imatrix_t *imatrix;
  imatrix_t *rmatrix;

  /* For sparse implementation */
  int hasSparse;
  int hasFull;
  imatrix_t *simatrix;
  imatrix_t *srmatrix;
  int simatrixLen;
  int *simatrixN;
  int *simatrixId;
  int **simatrixIdl;
  imatrix_t simatrixDefault;
  imatrix_t srmatrixDefault;

  int **cmtyl;
  int  *cmtyll;
  int *cmtyN;
  int *randomOrder;
  int *randomOrder2;

  LList_t seenList;
  GHashTable **cmtyListHash;

  } *Graph_t;

inline int cmtyListIsInCmty(Graph_t G, int c, int n);

double energy(Graph_t G, double gamma);
double energy_cmty(Graph_t G, double gamma, int c);
double energy_cmty_n(Graph_t G, double gamma, int c, int n);

double energy_sparse(Graph_t G, double gamma);
int minimize_sparse(Graph_t G, double gamma);
int combine_cmtys_sparse(Graph_t G, double gamma);




inline void LListAdd(LList_t lst, int value) {
  if (lst->count >= lst->maxcount) {
    return ; // We are already full, silently ignore adding it.
  }
  lst->data[lst->count] = value ;
  lst->count ++;
}
inline int LListContains(LList_t lst, int value) {
  int i = lst->count;
  //printf("Doing LlistLookup: n:%d\n", llist->n);
  while (i--) {
    //printf(" ..LlistLookup: %d %d\n", llist->n, i);
    if (lst->data[i] == value)
      return 1;
  }
  return 0;
}
void LListClear(LList_t lst) {
  lst->count = 0;
}
