#include "oct.h"
#include "oct_internal.h"
#include "seqalgorithms.h"

bool incrclosure_seq(dbm* m, size_t dim, size_t a, size_t b, bound_t d) {
  size_t i,j;
  size_t bara = a^1;
  size_t barb = b^1;


  // redundancy check
  if (bound_cmp(d, *getdbm(m,matpos2(a,b))) >= 0) {
    return false;
  }

  // sat check
  // m_{b,a} + d < 0
  // m_{bara, barb} < 0 (no need by coherence)
  // m_{bara,a} + d + m_{b,barb} + d < 0

  bound_t temp1; bound_init(temp1);
  bound_t temp2; bound_init(temp2);

  bound_add(temp1, *getdbm(m,matpos2(b,a)), d);
  bound_add(temp2, *getdbm(m,matpos2(bara,a)), d);
  bound_badd(temp2, *getdbm(m,matpos2(b,barb)));
  bound_badd(temp2, d);
  
  if ((bound_sgn(temp1) < 0) || (bound_sgn(temp2) < 0)) {
    bound_clear(temp1);
    bound_clear(temp2);
    return true;
  }

  bound_t twod; bound_init(twod);
  bound_mul_2(twod, d);
 
  bound_add(temp1, twod, *getdbm(m,matpos2(bara,a)));
  bound_add(temp2, twod, *getdbm(m,matpos2(b,barb)));

  size_t sz_mat = matsize(dim);
  if(STRONGINCR==1) strong_helper(m,dim,a,b,d,temp1,temp2);
  single_seq_helper(m, m, 0, sz_mat, a, b, d, temp1, temp2);
  bound_clear(temp1);
  bound_clear(temp2);
  bound_clear(twod);
  if(STRONGINCR==1) {
	return false;
  } else {
  	return hmat_s_step(m, dim);
  }
}

/* bool unary_inequality_incrclosure(bound_t* m, size_t dim, size_t a, bound_t d) { */

/*   size_t i,j; */
/*   size_t bara = a^1; */

/*   // redundancy ccheck */
/*   if (bound_cmp(d, m[matpos2(a,bara)]) >= 0) { */
/*     return false; */
/*   } */

/*   bound_t ij, ia, ibara; */

/*   bound_t tmp; bound_init(tmp); */
/*   bound_add(tmp, d, m[matpos2(bara,a)]); */

/*   if (bound_sgn(tmp) < 0) { */
/*     bound_clear(tmp); */
/*     return true; */
/*   } */
  
/*   bound_init(ij); */
/*   bound_init(ia); */


/*   for(i=0; i<2*dim; i++) { */
/*     size_t br = i|1; */
/*     bound_set(ia, m[matpos2(i,a)]); */
/*     bound_badd(ia, d); */
/*     for(j=0; j<=br; j++) { */
/*       bound_add(ij, ia, m[matpos2(bara,j)]); */
/*       bound_bmin(m[matpos2(i,j)], ij); */
/*     } */
/*   } */

/*   bound_clear(ij); */
/*   bound_clear(ia); */
/*   bound_clear(tmp); */

  
/*   return hmat_s_step(m,dim); */
/* } */

/* bool unary_equality_incrclosure(bound_t* m, size_t dim, size_t a, bound_t d, bound_t dprime ) { */
/*   size_t i, j; */
/*   size_t bara = a^1; */

/*   if (bound_cmp(d, m[matpos2(a,bara)]) >= 0 && bound_cmp(dprime, m[matpos2(bara,a)]) >= 0) { */
/*     return false; */
/*   } */

/*   bound_t ij, ia, ibara; */

/*   bound_t tmp, tmp2; */
/*   bound_init(tmp); */
/*   bound_init(tmp2); */

/*   bound_add(tmp, d, m[matpos2(bara,a)]); */
/*   bound_add(tmp2, dprime, m[matpos2(a,bara)]); */
/*   if ((bound_sgn(tmp) < 0) || (bound_sgn(tmp2) < 0)) { */
/*     bound_clear(tmp); */
/*     bound_clear(tmp2); */
/*     return true; */
/*   } */

/*   bound_init(ij); */
/*   bound_init(ia); */
/*   bound_init(ibara); */
/*   for(i=0; i<2*dim; i++) { */
/*     bound_set(ia, m[matpos2(i,a)]); */
/*     bound_set(ibara, m[matpos2(i,bara)]); */
/*     bound_badd(ia, d); */
/*     bound_badd(ibara, dprime); */
/*     for(j=0; j<=(i|1); j++)  { */
/*       bound_add(ij, ia, m[matpos2(bara,j)]); */
/*       bound_set(tmp, ij); */
/*       bound_add(ij, ibara,m[matpos2(a,j)]); */
/*       bound_bmin(tmp, ij); */
/*       bound_bmin(m[matpos2(i,j)], tmp); */
/*     } */
/*   } */

/*   bound_clear(ij); */
/*   bound_clear(ia); */
/*   bound_clear(ibara); */
/*   bound_clear(tmp); */
/*   bound_clear(tmp2); */
/*   return hmat_s_step(m,dim); */
/* } */

