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
  int hasPrimaryCmty;
  int const_q;

  int *cmty;
  imatrix_t *imatrix;
  imatrix_t *rmatrix;

  /* For sparse implementation */
  int hasSparse;
  int hasFull;
  imatrix_t *simatrix;
  int simatrixLen;
  int *simatrixN;
  int *simatrixId;
  int **simatrixIdl;
  imatrix_t simatrixDefault;
  imatrix_t *srmatrix;
  /* int srmatrixLen; */
  /* int *srmatrixN; */
  /* int *srmatrixId; */
  /* int **srmatrixIdl; */
  imatrix_t srmatrixDefault;
  imatrix_t srmatrixDefaultOnlyDefined;

  //int **cmtyl;
  //int  *cmtyll;
  int *cmtyN;
  int *randomOrder;
  int *randomOrder2;
  int *tmp;

  //LList_t seenList;
  GHashTable *seenList;
  GHashTable **cmtyListHash;

  } *Graph_t;

void hashInit(Graph_t G);
inline int cmtyListIsInCmty(Graph_t G, int c, int n);

double energy(Graph_t G, double gamma);
double energy_cmty(Graph_t G, double gamma, int c);
double energy_cmty_n(Graph_t G, double gamma, int c, int n);
double energy_cmty_n_sparse(Graph_t G, double gamma, int c, int n);

double energy_sparse(Graph_t G, double gamma);
int greedy_sparse(Graph_t G, double gamma);
int combine_sparse(Graph_t G, double gamma);




/* inline void SetAdd(LList_t lst, int value) { */
/*   if (lst->count >= lst->maxcount) { */
/*     return ; // We are already full, silently ignore adding it. */
/*   } */
/*   lst->data[lst->count] = value ; */
/*   lst->count ++; */
/* } */
/* inline int SetContains(LList_t lst, int value) { */
/*   int i = lst->count; */
/*   //printf("Doing LlistLookup: n:%d\n", llist->n); */
/*   while (i--) { */
/*     //printf(" ..LlistLookup: %d %d\n", llist->n, i); */
/*     if (lst->data[i] == value) */
/*       return 1; */
/*   } */
/*   return 0; */
/* } */
/* void SetClear(LList_t lst) { */
/*   lst->count = 0; */
/* } */
inline void SetAdd(GHashTable *HT, int value) {
  g_hash_table_insert(HT,
		      GINT_TO_POINTER(value),
		      NULL
		      );
}
inline int SetContains(GHashTable *HT, int value) {

  int found = g_hash_table_lookup_extended(HT, GINT_TO_POINTER(value),
					   NULL, NULL);
  return(found);
}
void SetClear(GHashTable *HT) {
  g_hash_table_remove_all(HT);
}
GHashTable *SetInit() {
  return (g_hash_table_new(g_direct_hash, g_direct_equal));
}
void SetDestroy(GHashTable *HT) {
  g_hash_table_destroy(HT);
}
