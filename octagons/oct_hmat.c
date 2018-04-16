/*
 * oct_hmat.c
 *
 * Half-matrices - Basic management.
 *
 * APRON Library / Octagonal Domain
 *
 * Copyright (C) Antoine Mine' 2006
 *
 */

/* This file is part of the APRON Library, released under LGPL license
   with an exception allowing the redistribution of statically linked
   executables.

   Please read the COPYING file packaged in the distribution.
*/

#include "oct.h"
#include "oct_internal.h"

/* We consider matrices of 2n*2n upper bounds.
   Let us denote by (i,j) the matrix element at line i, column j; the matrix
   induces the following constraints:
     Vj-Vi <= (2i,2j)
     Vj+Vi <= (2i+1,2j)
    -Vj-Vi <= (2i,2j+1)
    -Vj+Vi <= (2i+1,2j+1)

   Actually, this representation is redudant, and so, we manipulate 
   2x2 block lower-triangular matrices. 
   Only elements (i,j) such that j/2 <= i/2 are represented:

       j ->  0 1 2 3 4 5
            ___
        0  |_|_|
        1  |_|_|___
  i ->  2  |_|_|_|_|
        3  |_|_|_|_|___
        4  |_|_|_|_|_|_|
        5  |_|_|_|_|_|_|

                 
                 j
             0     -2x0
            2x0      0
       i
           x0-x1  -x0-x1      0   -2x1
           x0+x1  -x0+x1     2x1    0


   Elements such that j/2 > i/2 are retreived by coherence: (i,j) = (j^1,i^1)
*/


/* alloced but not initialized */

#define CACHESIZE 20000


struct 
{
  bound_t value;
  signed short index; 
} pairs[CACHESIZE];

bound_t values[CACHESIZE];

signed short dbmnext = 0;

void setinfty(dbm* d, size_t index);
void setzero(dbm* d, size_t index);

void printpairs(void) {
  signed short i;
  printf("[\n");
  for(i=0; i<dbmnext; i++) {
    printf("bound: "); bound_print(pairs[i].value);
    printf(" index: %d;\n", pairs[i].index);
  }
  printf("]\n");
}


void printvalues(void) {
  signed short i;
  printf("[\n");
  for(i=0; i<dbmnext; i++) {
    printf("bound: "); bound_print(values[i]); printf(" ");
  }
  printf("]\n");
}

bool check_values_invariant(void) {
  signed short i;
  for(i=0; i<dbmnext;i++) {
    if(bound_cmp(pairs[i].value, values[pairs[i].index]) != 0) return false;
  }
  return true;
}

void initcache(void) {
  assert(dbmnext == 0);
  bound_t zero, one, infty;
  bound_init(zero);
  bound_init(one); 
  bound_init(infty);
  
  bound_set_int(zero, 0);
  bound_set_int(one, 1);
  bound_set_infty(infty, 1);
  
  bound_set(values[0], zero);
  bound_set(values[1], one);
  bound_set(values[2], infty);
  
  pairs[0].index = 0; bound_set(pairs[0].value, zero);
  pairs[1].index = 1; bound_set(pairs[1].value, one);
  pairs[2].index = 2; bound_set(pairs[2].value, infty);
  
  dbmnext = 3;
}


#if defined(DBMCACHE)
inline dbm* hmat_alloc(oct_internal_t* pr, size_t dim)
{
  dbm* d;
  d = (dbm*) malloc(sizeof(dbm));
  unsigned short* r;
  size_t sz = matsize(dim);
  if (!sz) sz = 1; /* make sure we never malloc a O-sized block */
  checked_malloc(r,unsigned short,sz,return NULL;); 
  d->m = r;
  return d;
}

inline void hmat_free(oct_internal_t* pr, dbm* d, size_t dim)
{
  free(d->m);
  free(d);
}

/* all variables are initialized to 0 */
inline dbm* hmat_alloc_zero(oct_internal_t* pr, size_t dim)
{
  size_t i;
  dbm* d = hmat_alloc(pr,dim);
  for (i=0;i<matsize(dim);i++) setzero(d,i);
  return d;
}

/* all variables are initialized to "don't know" */
inline dbm* hmat_alloc_top(oct_internal_t* pr, size_t dim)
{
  size_t i;
  dbm* d = hmat_alloc(pr,dim);
  for (i=0;i<matsize(dim);i++) setinfty(d,i); 
  for (i=0;i<2*dim;i++) setzero(d,matpos(i,i)); 
  return d;
}
#else
inline dbm* hmat_alloc(oct_internal_t* pr, size_t dim)
{
  dbm* d;
  d = (dbm*) malloc(sizeof(dbm));
  bound_t* r;
  size_t sz = matsize(dim);
  if (!sz) sz = 1; /* make sure we never malloc a O-sized block */
  checked_malloc(r,bound_t,sz,return NULL;);
  bound_init_array(r,matsize(dim));
  d->m = r;
  return d;
}

inline void hmat_free(oct_internal_t* pr, dbm* d, size_t dim)
{
  bound_clear_array(d->m,matsize(dim));
  free(d->m);
  free(d);
}

/* all variables are initialized to 0 */
inline dbm* hmat_alloc_zero(oct_internal_t* pr, size_t dim)
{
  size_t i;
  dbm* d = hmat_alloc(pr,dim);
  for (i=0;i<matsize(dim);i++) bound_set_int(*getdbm(d,i),0);
  return d;
}

/* all variables are initialized to "don't know" */
inline dbm* hmat_alloc_top(oct_internal_t* pr, size_t dim)
{
  size_t i;
  dbm* d = hmat_alloc(pr,dim);
  for (i=0;i<matsize(dim);i++) bound_set_infty(*getdbm(d,i),1);
  for (i=0;i<2*dim;i++) bound_set_int(*getdbm(d,matpos(i,i)),0);
  return d;
}

#endif

void hmat_fdump(FILE* stream, oct_internal_t* pr, dbm* d, size_t dim)
{
  size_t i,j;
  for (i=0;i<2*dim;i++) {
    for (j=0;j<=(i|1);j++) {
      if (j) fprintf(stream," ");
      bound_fprint(stream,*getdbm(d, matpos2(i, j)));
    }
    fprintf(stream,"\n");
  }
}
inline dbm* hmat_copy(oct_internal_t* pr, dbm* m, size_t dim)
{
  if (m) {
    dbm* d = hmat_alloc(pr,dim);
    size_t i;
    for(i=0; i<matsize(dim); i++) {
      setdbm(d,i, *getdbm(m,i));
    }
    return d;
  }
  else {
    return NULL;
  }
}




signed short dbm_insert(bound_t new)
{
  signed short lower = 0;
  signed short upper = dbmnext - 1;
  assert(dbmnext > 0);
  
  while (lower <= upper )
    {
      signed short mid = lower + ((upper - lower)/2);
      int cmp = bound_cmp(new, pairs[mid].value); // 1 if new>pairs[mid].value
      if (cmp < 0) upper = mid-1;
      else if (cmp > 0) lower = mid+1;
      else return pairs[mid].index;
    }
  
  if (dbmnext >= CACHESIZE) {
    return lower + 1;
  } else {
    for (unsigned short i = dbmnext; i > lower; i--)
      {
	bound_set(pairs[i].value, pairs[i-1].value);
	pairs[i].index = pairs[i - 1].index;
      }
    
    bound_set(values[dbmnext], new);
    
    bound_set(pairs[lower].value, new);
    pairs[lower].index = dbmnext;
    dbmnext++;
    return pairs[lower].index;
  }
}

#if defined(DBMCACHE)
void setinfty(dbm* d, size_t index) {
  d->m[index] = 2;
}

void setzero(dbm* d, size_t index) {
  d->m[index] = 0;
}
#endif

#if defined(DBMCACHE)
inline void setdbm(dbm* d, size_t k, bound_t new) {
  signed short value_index = dbm_insert(new); 
  d->m[k] = value_index;
}

void setdbminfty(dbm* d, size_t k) {
  d->m[k] = 2;
}

void setdbmzero(dbm* d, size_t k) {
   d->m[k] = 0;
}
#else
inline void setdbm(dbm* d, size_t k, bound_t new) {
  bound_set(d->m[k], new);
}

inline void setdbminfty(dbm* d, size_t k) {
  bound_set_infty(*getdbm(d,k), 1);
}

inline void setdbmzero(dbm* d, size_t k) {
  bound_set_int(*getdbm(d,k), 0);
}
#endif

void setdbmbmin(dbm* d, size_t k, bound_t b) {
  if(bound_cmp(*getdbm(d,k), b) > 0) {
    setdbm(d,k,b);
  }
}

void setdbmmin(dbm* d, size_t k, bound_t b1, bound_t b2) {
  bound_t tmp; bound_init(tmp);
  bound_min(tmp, b1, b2);
  setdbm(d, k, tmp);
}

void setdbmbmax(dbm* d, size_t k, bound_t b) {
  if(bound_cmp(*getdbm(d,k), b) < 0) {
    setdbm(d,k,b);
  }
}

void setdbmmax(dbm* d, size_t k, bound_t b1, bound_t b2) {
  bound_t tmp; bound_init(tmp);
  bound_max(tmp, b1, b2);
  setdbm(d,k,tmp);
}

void setdbmadd(dbm* d, size_t k, bound_t b1, bound_t b2) {
  bound_t tmp; bound_init(tmp);
  bound_add(tmp, b1, b2);
  setdbm(d,k,tmp);
}

void setdbmsub(dbm* d, size_t k, bound_t b1, bound_t b2) {
  bound_t tmp; bound_init(tmp);
  bound_sub(tmp, b1,b2);
  setdbm(d,k,tmp);
}

void setdbmbadd(dbm* d, size_t k, bound_t b) {
  bound_t tmp; bound_init(tmp);
  bound_set(tmp, *getdbm(d,k));
  bound_badd(tmp,b);
  setdbm(d ,k, tmp);
}

void setdbm_mul_2(dbm* d, size_t k, bound_t b) {
  bound_t tmp; bound_init(tmp);
  bound_mul_2(tmp, b);
  setdbm(d,k,tmp);
}


void setdbm_div_2(dbm* d, size_t k, bound_t b) {
  bound_t tmp; bound_init(tmp);
  bound_div_2(tmp, b);
  setdbm(d,k,tmp);
}

#if defined(DBMCACHE)
inline bound_t* getdbm(dbm* d, size_t k) {
   return &values[d->m[k]];
}
#else
inline bound_t* getdbm(dbm* d, size_t k) {
  return (d->m)+k;
}
#endif 

inline void dbm_set_array(dbm* dst, dbm* src, size_t size) {
  size_t k;
  for(k=0; k<size; k++) {
    setdbm(dst, k, *getdbm(src, k));
  }
}

inline void dbm_set_array_from_point(dbm* dst, dbm* src, size_t point, size_t point2,size_t size) {
  size_t k;
  for(k=0; k<size; k++) {
    setdbm(dst, point+k, *getdbm(src, point2+k));
  }
}

inline void dbm_bound_set_array(bound_t* dst, dbm* src, size_t size) {
  size_t k;
  for(k=0; k<size; k++) {
    bound_set(dst[k], *getdbm(src,k));
  }
}

inline void dbm_bound_set_array2(dbm* dst, bound_t* src, size_t size) {
  size_t k;
  for(k=0; k<size; k++) {
    setdbm(dst,k, src[k]);
  }
}

inline size_t dbm_serialized_size_array(dbm* src, size_t size) {
  size_t i, n=0;
  for(i=0; i<size; i++) {
    n+= bound_serialized_size(*getdbm(src, i));
  }
  return n;
}

inline size_t dbm_serialize_array(void* dst, dbm* src, size_t size) {
  size_t i,n=0;
  for (i=0;i<size;i++)
    n += bound_serialize((char*)dst+n,*getdbm(src,i));
  return n;
}

inline size_t dbm_deserialize_array(dbm* dst, const void *src, size_t size) {
  size_t i,n=0;
  for (i=0;i<size;i++)
    n += bound_deserialize(*getdbm(dst,i),(const char*)src+n);
  return n;
}
