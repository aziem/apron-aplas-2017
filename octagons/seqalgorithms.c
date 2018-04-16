#include "oct.h"
#include "oct_internal.h"

bool single_seq_helper(dbm* m, dbm* old, size_t start, size_t end, size_t a, size_t b, bound_t d, bound_t temp1, bound_t temp2) {
    bound_t temp5; bound_init(temp5);
    bound_t temp3; bound_init(temp3);
    bound_t temp4; bound_init(temp4);
    bound_t temp6; bound_init(temp6);

    bound_t mibari; bound_init(mibari);

    size_t i;
    size_t j;

    i = (size_t) (sqrt((2*start)+1)-1);
    j = start - ((i+1) * (i+1))/2;

    size_t ibarb = matpos2(i,b^1);
    size_t ia = matpos2(i,a);
    bound_set(mibari, *getdbm(old,matpos(i,i^1)));
    bound_set(temp5,d);
    bound_badd(temp5, *getdbm(old, ia));

    bound_set(temp3, temp1);
    bound_badd(temp3, *getdbm(old, ibarb));
    bound_bmin(temp3, temp5);

    bound_add(temp5,*getdbm(old, ibarb), d);

    bound_add(temp4, *getdbm(old, ia), temp2);
    bound_bmin(temp4, temp5);

    for (size_t k = start; k < end; k++)
    {
	if (j < (i^1))
	{
	  bound_add(temp5, temp3, *getdbm(old, matpos2(b,j)));
	  bound_add(temp6, temp4, *getdbm(old,matpos2(a^1,j)));
	  bound_bmin(temp6, temp5);

	  if (STRONGINCR == 1)
	    {
	      bound_add(temp5, mibari, *getdbm(old,matpos(j^1,j)));
	      bound_div_2(temp5, temp5);
	      bound_bmin(temp6, temp5);
	    }
	  
	  if (bound_cmp(*getdbm(m,k), temp6) > 0) {
	    setdbm(m,k,temp6);
	  }
	
	  
	  j = j + 1;
	}
	else
	{
	  if (j == (i^1))
	    {
	      if (STRONGINCR == 0)
		{
		  bound_add(temp5, temp3, *getdbm(old, matpos2(b,j)));
		  bound_add(temp6, temp4, *getdbm(old,matpos2(a^1,j)));
		  bound_bmin(temp6, temp5);
		  
		  if (bound_cmp(*getdbm(m,k), temp6) > 0) {
		    setdbm(m,k,temp6);
		  }
		}
	    }
	  else
	    {
	      bound_add(temp5, temp3, *getdbm(old,matpos2(b,j)));
	      bound_add(temp6, temp4, *getdbm(old,matpos2(a^1,j)));
	      bound_bmin(temp6, temp5);
	      
	      if (STRONGINCR == 1)
		{
		  bound_add(temp5, mibari, *getdbm(old, matpos(j^1,j)));
		  bound_div_2(temp5, temp5);
		  bound_bmin(temp6, temp5);
		}
	      
	      if (bound_cmp(*getdbm(m,k), temp6) > 0) {
		setdbm(m, k, temp6);
	      }
	    }
	  
	  if (j < (i|1))
	    {
	      j = j+1;
	    }
	  else if (k<end-1)
	    {
	      i = i+1;
	      j = 0;
	      size_t ibarb = matpos2(i,b^1);
	      size_t ia = matpos2(i,a);
	      bound_set(mibari, *getdbm(old,matpos(i,i^1)));
	      bound_set(temp5,d);
	      bound_badd(temp5, *getdbm(old,ia));
	      
	      bound_set(temp3, temp1);
	      bound_badd(temp3, *getdbm(old,ibarb));
	      bound_bmin(temp3, temp5);
	      
	      bound_add(temp5,*getdbm(old,ibarb), d);
	      
	      bound_add(temp4, *getdbm(old,ia), temp2);
	      bound_bmin(temp4, temp5);
	    }
	}
    }
    
    bound_clear(temp3);
    bound_clear(temp4);
    bound_clear(temp5);
    bound_clear(temp6);
    bound_clear(mibari);

    return false;
}

void strong_helper(dbm* m, size_t dim,  size_t a, size_t b, bound_t d, bound_t temp1, bound_t temp2){
  bound_t temp3, temp4;
  bound_t ia, ibarb;
  size_t k; 

  bound_init(ia);
  bound_init(ibarb);
  bound_init(temp3);
  bound_init(temp4);

  for(size_t i=0; i<2*dim; i++) {
    bound_set(ia, *getdbm(m,matpos2(i,a)));
    bound_set(ibarb, *getdbm(m,matpos2(i,b^1)));

    bound_add(temp3, ia, d);
    bound_badd(temp3, ibarb);

    bound_mul_2(temp4, ibarb);
    bound_badd(temp4, temp1);
    bound_bmin(temp3, temp4);

    bound_mul_2(temp4, ia);
    bound_badd(temp4, temp2);
    bound_bmin(temp3, temp4);
    k = matpos(i,i^1);

    if(bound_cmp(*getdbm(m,k), temp3) > 0) {
      setdbm(m,k,temp3);
    }
  }

  bound_clear(temp3);
  bound_clear(temp4);
  bound_clear(ibarb);
  bound_clear(ia);
  return;
}

bool double_seq_helper(bound_t* m, bound_t* old, size_t dim, size_t start, size_t end, size_t a, size_t b, bound_t d, bound_t temp1, bound_t temp2) {

  size_t i,j;
  size_t bara = a^1;
  size_t barb = b^1;

  bound_t temp3; bound_init(temp3);
  bound_t temp4; bound_init(temp4);
  bound_t temp5; bound_init(temp5);
  bound_t temp6; bound_init(temp6);
 
  bound_t ia; bound_init(ia);
  bound_t ibarb; bound_init(ibarb);

  size_t sz_mat = matsize(dim);
  assert(start==0);
  for(i=0; i<2*dim; i++) {
    bound_set(ia, m[matpos2(i,a)]);
    bound_set(ibarb, m[matpos2(i,barb)]);

    bound_add(temp5, ia, d);
    bound_add(temp3, ibarb, temp1);
    bound_min(temp3, temp5, temp3);

    bound_add(temp5, ibarb, d);
    bound_add(temp4, ia, temp2);
    bound_min(temp4, temp5, temp4);

    for(j=0; j<=(i|1); j++) {

      size_t ij = matpos2(i,j);
      
      bound_add(temp5, temp3, m[matpos2(b,j)]);
      bound_add(temp6, temp4, m[matpos2(bara,j)]);

      bound_bmin(temp6, temp5);

      if(bound_cmp(m[ij], temp6) > 0) {
	bound_set(m[ij], temp6);
      }
    }
  }

  bound_clear(temp3);
  bound_clear(temp4);
  bound_clear(temp5);
  bound_clear(temp6);
  bound_clear(ia);
  bound_clear(ibarb);
  return false;
}
