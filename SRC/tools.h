#include "./pblas.h"

#ifdef __STDC__
typedef void (*CPYPTR)(int, int, float *, int, float *, int);
#define SLVOID void
#else
typedef void (*CPYPTR)();
#define SLVOID char
#endif

#define ErrPrnt fprintf(stderr, "line %d of file %s\n",__LINE__, __FILE__);

#define Mdescset(desc, m, n, mb, nb, rsrc, csrc, ictxt, lld) \
{ \
   (desc)[DT_] = BLOCK_CYCLIC_2D; \
   (desc)[CTXT_] = (ictxt); \
   (desc)[M_] = (m); \
   (desc)[N_] = (n); \
   (desc)[MB_] = (mb); \
   (desc)[NB_] = (nb); \
   (desc)[RSRC_] = (rsrc); \
   (desc)[CSRC_] = (csrc); \
   (desc)[LLD_] = (lld); \
}

#define MCindxg2p(IG, nb, srcproc, nprocs) \
        ( ((srcproc) + (IG)/(nb)) % nprocs )

typedef struct {double r, i;} DCOMPLEX;
typedef struct {float r, i;} SCOMPLEX;

#define Mmalloc(M_ptr, M_type, M_elt, M_i, M_ctxt) \
{ \
   void pberror_(); \
   (M_ptr) = ( M_type * ) malloc((M_elt)*(sizeof(M_type))); \
   if (!(M_ptr)) \
   { \
      if ((M_elt) > 0) \
      { \
         (M_i) = 1; \
         fprintf(stderr, "Not enough memory on line %d of file %s!!\n", \
                 __LINE__, __FILE__); \
         pberror_(&(M_ctxt), __FILE__, &(M_i)); \
      } \
   } \
}
