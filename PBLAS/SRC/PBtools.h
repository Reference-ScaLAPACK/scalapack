/* ---------------------------------------------------------------------
*
*  -- PBLAS routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     March 12, 2002 
*
*  ---------------------------------------------------------------------
*/
/*
*  This file includes PBLAS tools definitions. All PBLAS routines include
*  this file.
*
* ----------------------------------------------------------------------
*  #define macro constants
*  ---------------------------------------------------------------------
*
*  Descriptor entries for type 1
*/
#define    BLOCK_CYCLIC_2D     1

#define    DTYPE1_             0                   /* Descriptor Type */
#define    CTXT1_              1                     /* BLACS context */
#define    M1_                 2             /* Global Number of Rows */
#define    N1_                 3          /* Global Number of Columns */
#define    MB1_                4                 /* Row Blocking Size */
#define    NB1_                5              /* Column Blocking Size */
#define    RSRC1_              6            /* Starting Processor Row */
#define    CSRC1_              7         /* Starting Processor Column */
#define    LLD1_               8           /* Local Leading Dimension */
#define    DLEN1_              9                 /* Descriptor Length */
/*
*  Descriptor entries for type 2
*/
#define    BLOCK_CYCLIC_2D_INB 2

#define    DTYPE_              0                   /* Descriptor Type */
#define    CTXT_               1                     /* BLACS context */
#define    M_                  2             /* Global Number of Rows */
#define    N_                  3          /* Global Number of Columns */
#define    IMB_                4         /* Initial Row Blocking Size */
#define    INB_                5      /* Initial Column Blocking Size */
#define    MB_                 6                 /* Row Blocking Size */
#define    NB_                 7              /* Column Blocking Size */
#define    RSRC_               8              /* Starting Process Row */
#define    CSRC_               9           /* Starting Process Column */
#define    LLD_                10          /* Local Leading Dimension */
#define    DLEN_               11                /* Descriptor Length */

#define    CPACKING            'P'
#define    CUNPACKING          'U'

#define    PACKING             "P"
#define    UNPACKING           "U"

#define    CGENERAL            'G'
/* #define    CSYMM               'S'  */
#define    CHERM               'H'

#define    GENERAL             "G"
#define    SYMM                "S"
#define    HERM                "H"

#define    ONE                 1.0
#define    TWO                 2.0
#define    ZERO                0.0
                            /* Input error checking related constants */
#define    DESCMULT            100
#define    BIGNUM              10000
/*
*  ---------------------------------------------------------------------
*  #define macro functions
*  ---------------------------------------------------------------------
*/
#define    ABS( a_ )           ( ( (a_) <   0  ) ? -(a_) : (a_) )
#define    MIN( a_, b_ )       ( ( (a_) < (b_) ) ?  (a_) : (b_) )
#define    MAX( a_, b_ )       ( ( (a_) > (b_) ) ?  (a_) : (b_) )

#define    FLOOR(a,b) (((a)>0) ? (((a)/(b))) : (-(((-(a))+(b)-1)/(b))))
#define    CEIL(a,b)           ( ( (a)+(b)-1 ) / (b) )
#define    ICEIL(a,b) (((a)>0) ? ((((a)+(b)-1)/(b))) : (-((-(a))/(b))))

#define    Mupcase(C)          (((C)>96 && (C)<123) ? (C) & 0xDF : (C))
#define    Mlowcase(C)         (((C)>64 && (C)< 91) ? (C) | 32   : (C))
/*
*  The following macros perform common modulo operations;  All functions
*  except MPosMod assume arguments are < d (i.e., arguments are themsel-
*  ves within modulo range).
*/
                                                /* increment with mod */
#define    MModInc(I, d)       if(++(I) == (d)) (I) = 0
                                                /* decrement with mod */
#define    MModDec(I, d)       if(--(I) == -1) (I) = (d)-1
                                                   /* positive modulo */
#define    MPosMod(I, d)       ( (I) - ((I)/(d))*(d) )
                                                   /* add two numbers */
#define    MModAdd(I1, I2, d) \
           ( ( (I1) + (I2) < (d) ) ? (I1) + (I2) : (I1) + (I2) - (d) )
                                                        /* add 1 to # */
#define    MModAdd1(I, d) ( ((I) != (d)-1) ? (I) + 1 : 0 )
                                              /* subtract two numbers */
#define    MModSub(I1, I2, d) \
           ( ( (I1) < (I2) ) ? (d) + (I1) - (I2) : (I1) - (I2) )
                                                      /* sub 1 from # */
#define    MModSub1(I, d) ( ((I)!=0) ? (I)-1 : (d)-1 )
/*
*  DNROC computes maximum number of local rows or columns. This macro is
*  only used to compute the time estimates in the Level 3 PBLAS routines.
*/

#define    DNROC( n_, nb_, p_ ) \
           ((double)(((((n_)+(nb_)-1)/(nb_))+(p_)-1)/(p_))*(double)((nb_)))
/*
*  Mptr returns a pointer to a_( i_, j_ ) for readability reasons and
*  also less silly errors ...
*
*  There was some problems with the previous code which read:
*
*      #define    Mptr( a_, i_, j_, lda_, siz_ ) \
*                    ( (a_) + ( ( (i_)+(j_)*(lda_) )*(siz_) ) )
* 
*  since it can overflow the 32-bit integer "easily".
*  The following code should fix the problem.
*  It uses the "off_t" command.
*
*  Change made by Julien Langou on Sat. September 12, 2009. 
*  Fix provided by John Moyard from CNES.
*
*  JL :April 2011: Change off_t by long long
*  off_t is not supported under Windows
*/
#define    Mptr( a_, i_, j_, lda_, siz_ ) \
              ( (a_) + ( (long long) ( (long long)(i_)+ \
              (long long)(j_)*(long long)(lda_))*(long long)(siz_) ) )
/*
*  Mfirstnb and Mlastnb compute the global size of the first and last
*  block corresponding to the interval i_:i_+n_-1 of global indexes.
*/
#define    Mfirstnb( inbt_, n_, i_, inb_, nb_ ) \
              inbt_ = (inb_) - (i_); \
              if( inbt_ <= 0 ) \
                 inbt_ = ( (-inbt_) / (nb_) + 1 ) * (nb_) + inbt_; \
              inbt_ = MIN( inbt_, (n_) );

#define    Mlastnb( inbt_, n_, i_, inb_, nb_ ) \
              inbt_ = (i_) + (n_) - (inb_); \
              if( inbt_ > 0 ) \
              { \
                 inbt_ = -( ( (nb_)+inbt_-1 )/(nb_)-1 )*(nb_) + inbt_; \
                 inbt_ = MIN( inbt_, (n_) ); \
              } \
              else { inbt_ = (n_); };
/*
*  Does the index interval i_:i_+n_-1 spans more than one process rows
*  or columns ?
*
*  Mspan returns 0 (false) when the data is replicated (srcproc_ < 0) or
*  when there is only one process row or column in the process grid.
*/
#define    Mspan( n_, i_, inb_, nb_, srcproc_, nprocs_ ) \
              ( ( (srcproc_) >= 0 ) && ( ( (nprocs_) > 1 ) && \
              ( ( (i_) < (inb_) ) ? \
                ( (i_) + (n_) > (inb_) ) : \
                ( (i_) + (n_) > (inb_) + \
                  ( ( (i_) - (inb_) ) / (nb_) + 1 ) * (nb_) ) ) ) )
/*
*  Mindxl2g computes the global index ig_ corresponding to the local
*  index il_ in process proc_.
*/
#define    Mindxl2g( ig_, il_, inb_, nb_, proc_, srcproc_, nprocs_ ) \
           { \
              if( ( (srcproc_) >= 0 ) && ( (nprocs_) > 1 ) ) \
              { \
                 if( (proc_) == (srcproc_) ) \
                 { \
                    if( (il_) < (inb_) ) ig_ = (il_); \
                    else                 ig_ = (il_) + \
                       (nb_)*((nprocs_)-1)*( ((il_)-(inb_)) / (nb_) + 1 ); \
                 } \
                 else if( (proc_) < (srcproc_) ) \
                 { \
                    ig_ = (il_) + (inb_) + \
                          (nb_)*(  ((nprocs_)-1)*((il_)/(nb_)) + \
                                   (proc_)-(srcproc_)-1+(nprocs_) ); \
                 } \
                 else \
                 { \
                    ig_ =  (il_) + (inb_) + \
                           (nb_)*( ((nprocs_)-1)*((il_)/(nb_)) + \
                           (proc_)-(srcproc_)-1 ); \
                 } \
              } \
              else \
              { \
                 ig_ = (il_); \
              } \
           }
/*
*  Mindxg2p returns the process coodinate owning the entry globally
*  indexed by ig_.
*/
#define    Mindxg2p( ig_, inb_, nb_, proc_, srcproc_, nprocs_ ) \
           { \
              if( ( (ig_) >= (inb_) ) && ( (srcproc_) >= 0 ) && \
                  ( (nprocs_) > 1 ) ) \
              { \
                 proc_  = (srcproc_) + 1 + ( (ig_)-(inb_) ) / (nb_); \
                 proc_ -= ( proc_ / (nprocs_) ) * (nprocs_); \
              } \
              else \
              { \
                 proc_ = (srcproc_); \
              } \
           }
/*
*  Mnumroc computes the # of local indexes np_ residing in the process
*  of coordinate proc_ corresponding to the interval of global indexes
*  i_:i_+n_-1 assuming that the global index 0 resides in  the process
*  srcproc_, and that the indexes are distributed from  srcproc_ using
*  the parameters inb_, nb_ and nprocs_.
*/
#define    Mnumroc( np_, n_, i_, inb_, nb_, proc_, srcproc_, nprocs_ ) \
           { \
              if( ( (srcproc_) >= 0 ) && ( (nprocs_) > 1 ) ) \
              { \
                 Int inb__, mydist__, n__, nblk__, quot__, src__; \
                 if( ( inb__ = (inb_) - (i_) ) <= 0 ) \
                 { \
                    src__  = (srcproc_) + ( nblk__ = (-inb__) / (nb_) + 1 ); \
                    src__ -= ( src__ / (nprocs_) ) * (nprocs_); \
                    inb__ += nblk__*(nb_); \
                    if( ( n__ = (n_) - inb__ ) <= 0 ) \
                    { if( (proc_) == src__ ) np_ = (n_); else np_ = 0; } \
                    else \
                    { \
                       if( ( mydist__ = (proc_) - src__ ) < 0 ) \
                          mydist__ += (nprocs_); \
                       nblk__    = n__ / (nb_) + 1; \
                       mydist__ -= nblk__ - \
                          ( quot__ = ( nblk__ / (nprocs_) ) ) * (nprocs_); \
                       if( mydist__ < 0 ) \
                       { \
                          if( (proc_) != src__ ) \
                             np_ = (nb_) + (nb_) * quot__; \
                          else \
                             np_ = inb__ + (nb_) * quot__; \
                       } \
                       else if( mydist__ > 0 ) \
                       { \
                          np_ = (nb_) * quot__; \
                       } \
                       else \
                       { \
                          if( (proc_) != src__ ) \
                             np_ = n__ + (nb_) + (nb_) * ( quot__ - nblk__ ); \
                          else \
                             np_ = (n_) +        (nb_) * ( quot__ - nblk__ ); \
                       } \
                    } \
                 } \
                 else \
                 { \
                    if( ( n__ = (n_) - inb__ ) <= 0 ) \
                    { if( (proc_) == (srcproc_) ) np_ = (n_); else np_ = 0; } \
                    else \
                    { \
                       if( ( mydist__ = (proc_) - (srcproc_) ) < 0 ) \
                          mydist__ += (nprocs_); \
                       nblk__    = n__ / (nb_) + 1; \
                       mydist__ -= nblk__ - \
                          ( quot__ = ( nblk__ / (nprocs_) ) ) * (nprocs_); \
                       if( mydist__ < 0 ) \
                       { \
                          if( (proc_) != (srcproc_) ) \
                             np_ = (nb_) + (nb_) * quot__; \
                          else \
                             np_ = inb__ + (nb_) * quot__; \
                       } \
                       else if( mydist__ > 0 ) \
                       { \
                          np_ = (nb_) * quot__; \
                       } \
                       else \
                       { \
                          if( (proc_) != (srcproc_) ) \
                             np_ = n__ + (nb_) + (nb_) * ( quot__ - nblk__ ); \
                          else \
                             np_ = (n_) +        (nb_) * ( quot__ - nblk__ ); \
                       } \
                    } \
                 } \
              } \
              else \
              { \
                 np_ = (n_); \
              } \
           }

#define    Mnpreroc( np_, n_, i_, inb_, nb_, proc_, srcproc_, nprocs_ ) \
           { \
              if( ( (srcproc_) >= 0 ) && ( (nprocs_) > 1 ) ) \
              { \
                 Int inb__, mydist__, n__, nblk__, quot__, rem__, src__; \
                 if( ( inb__ = (inb_) - (i_) ) <= 0 ) \
                 { \
                    src__  = (srcproc_) + ( nblk__ = (-inb__) / (nb_) + 1 ); \
                    src__ -= ( src__ / (nprocs_) ) * (nprocs_); \
                    if( (proc_) != src__ ) \
                    { \
                       inb__ += nblk__*(nb_); \
                       if( ( n__ = (n_) - inb__ ) <= 0 ) { np_ = (n_); } \
                       else \
                       { \
                          if( ( mydist__ = (proc_) - src__ ) < 0 ) \
                             mydist__ += (nprocs_); \
                          nblk__ = n__ / (nb_) + 1; \
                          rem__  = nblk__ - \
                             ( quot__ = ( nblk__ / (nprocs_) ) ) * (nprocs_); \
                          if( mydist__ <= rem__ ) \
                          { \
                             np_ = inb__ - (nb_) + \
                                   ( quot__ + 1 ) * mydist__ * (nb_); \
                          } \
                          else \
                          { \
                             np_ = (n_) + \
                                   ( mydist__ - (nprocs_) ) * quot__ * (nb_); \
                          } \
                       } \
                    } \
                    else \
                    { \
                       np_ = 0; \
                    } \
                 } \
                 else \
                 { \
                    if( (proc_) != (srcproc_) ) \
                    { \
                       if( ( n__ = (n_) - inb__ ) <= 0 ) { np_ = (n_); } \
                       else \
                       { \
                          if( ( mydist__ = (proc_) - (srcproc_) ) < 0 ) \
                             mydist__ += (nprocs_); \
                          nblk__ = n__ / (nb_) + 1; \
                          rem__  = nblk__ - \
                             ( quot__ = ( nblk__ / (nprocs_) ) ) * (nprocs_); \
                          if( mydist__ <= rem__ ) \
                          { \
                             np_ = inb__ - (nb_) + \
                                   ( quot__ + 1 ) * mydist__ * (nb_); \
                          } \
                          else \
                          { \
                             np_ = (n_) + \
                                   ( mydist__ - (nprocs_) ) * quot__ * (nb_); \
                          } \
                       } \
                    } \
                    else \
                    { \
                       np_ = 0; \
                    } \
                 } \
              } \
              else \
              { \
                 np_ = 0; \
              } \
           }

#define    Mnnxtroc( np_, n_, i_, inb_, nb_, proc_, srcproc_, nprocs_ ) \
           { \
              if( ( (srcproc_) >= 0 ) && ( (nprocs_) > 1 ) ) \
              { \
                 Int inb__, mydist__, n__, nblk__, quot__, rem__, src__; \
                 if( ( inb__ = (inb_) - (i_) ) <= 0 ) \
                 { \
                    src__  = (srcproc_) + ( nblk__ = (-inb__) / (nb_) + 1 ); \
                    src__ -= ( src__ / (nprocs_) ) * (nprocs_); \
                    inb__ += nblk__*(nb_); \
                    if( ( n__ = (n_) - inb__ ) <= 0 ) { np_ = 0; } \
                    else \
                    { \
                       if( ( mydist__ = (proc_) - src__ ) < 0 ) \
                          mydist__ += (nprocs_); \
                       nblk__ = n__ / (nb_) + 1; \
                       rem__  = nblk__ - \
                             ( quot__ = ( nblk__ / (nprocs_) ) ) * (nprocs_); \
                       if( mydist__ < rem__ ) \
                       { \
                          np_ = n__ - ( quot__ * mydist__ + \
                                        quot__ + mydist__ ) * (nb_); \
                       } \
                       else \
                       { \
                          np_ = ( (nprocs_) - 1 - mydist__ ) * quot__ * (nb_); \
                       } \
                    } \
                 } \
                 else \
                 { \
                    if( ( n__ = (n_) - inb__ ) <= 0 ) { np_ = 0; } \
                    else \
                    { \
                       if( ( mydist__ = (proc_) - (srcproc_) ) < 0 ) \
                          mydist__ += (nprocs_); \
                       nblk__ = n__ / (nb_) + 1; \
                       rem__  = nblk__ - \
                             ( quot__ = ( nblk__ / (nprocs_) ) ) * (nprocs_); \
                       if( mydist__ < rem__ ) \
                       { \
                          np_ = n__ - ( quot__ * mydist__ + \
                                        quot__ + mydist__ ) * (nb_); \
                       } \
                       else \
                       { \
                          np_ = ( (nprocs_) - 1 - mydist__ ) * quot__ * (nb_); \
                       } \
                    } \
                 } \
              } \
              else \
              { np_ = 0; } \
           }


#define    Minfog2l( i_, j_, desc_, nr_, nc_, r_, c_, ii_, jj_, pr_, pc_ ) \
           { \
              Int quot__, i__, imb__, inb__, j__, mb__, mydist__, \
                  nb__, nblk__, src__; \
              imb__ = desc_[IMB_]; mb__ = desc_[MB_]; pr_ = desc_[RSRC_]; \
              if( ( pr_ >= 0 ) && ( nr_ > 1 ) ) \
              { \
                 if( ( i__ = (i_) - imb__ ) < 0 ) \
                 { ii_ = ( r_ == pr_ ? (i_) : 0 ); } \
                 else \
                 { \
                    src__     = pr_; \
                    pr_      += ( nblk__ = i__ / mb__ + 1 ); \
                    pr_      -= ( pr_ / nr_ ) * nr_; \
                    if( ( mydist__ = r_ - src__ ) < 0 ) mydist__ += nr_; \
                    if( mydist__ >= nblk__ - ( quot__ = nblk__ / nr_ ) * nr_ ) \
                    { \
                       if( r_ != src__ ) ii_ =  mb__; \
                       else              ii_ = imb__; \
                       if( r_ != pr_ ) \
                          ii_ += ( quot__ - 1 ) * mb__; \
                       else \
                          ii_ += i__ + ( quot__ - nblk__ ) * mb__; \
                    } \
                    else \
                    { \
                       if( r_ != src__ ) ii_ =  mb__ + quot__ * mb__; \
                       else              ii_ = imb__ + quot__ * mb__; \
                    } \
                 } \
              } \
              else \
              { \
                 ii_ = (i_); \
              } \
              inb__ = desc_[INB_]; nb__ = desc_[NB_]; pc_ = desc_[CSRC_]; \
              if( ( pc_ >= 0 ) && ( nc_ > 1 ) ) \
              { \
                 if( ( j__ = (j_) - inb__ ) < 0 ) \
                 { jj_ = ( c_ == pc_ ? (j_) : 0 ); } \
                 else \
                 { \
                    src__     = pc_; \
                    pc_      += ( nblk__ = j__ / nb__ + 1 ); \
                    pc_      -= ( pc_ / nc_ ) * nc_; \
                    if( ( mydist__ = c_ - src__ ) < 0 ) mydist__ += nc_; \
                    if( mydist__ >= nblk__ - ( quot__ = nblk__ / nc_ ) * nc_ ) \
                    { \
                       if( c_ != src__ ) jj_ =  nb__; \
                       else              jj_ = inb__; \
                       if( c_ != pc_ ) \
                          jj_ += ( quot__ - 1 ) * nb__; \
                       else \
                          jj_ += j__ + ( quot__ - nblk__ ) * nb__; \
                    } \
                    else \
                    { \
                       if( c_ != src__ ) jj_ =  nb__ + quot__ * nb__; \
                       else              jj_ = inb__ + quot__ * nb__; \
                    } \
                 } \
              } \
              else \
              { \
                 jj_ = (j_); \
              } \
           }

/*
*  The following macros initialize or translate descriptors.
*/
#define    MDescSet( desc, m, n, imb, inb, mb, nb, rsrc, csrc, ictxt, lld ) \
           { \
              (desc)[DTYPE_] = BLOCK_CYCLIC_2D_INB; \
              (desc)[CTXT_ ] = (ictxt); \
              (desc)[M_    ] = (m);     \
              (desc)[N_    ] = (n);     \
              (desc)[IMB_  ] = (imb);   \
              (desc)[INB_  ] = (inb);   \
              (desc)[MB_   ] = (mb);    \
              (desc)[NB_   ] = (nb);    \
              (desc)[RSRC_ ] = (rsrc);  \
              (desc)[CSRC_ ] = (csrc);  \
              (desc)[LLD_  ] = (lld);   \
           }

#define    MDescCopy(DescIn, DescOut) \
           { \
              (DescOut)[DTYPE_] = (DescIn)[DTYPE_];    \
              (DescOut)[M_    ] = (DescIn)[M_    ];    \
              (DescOut)[N_    ] = (DescIn)[N_    ];    \
              (DescOut)[IMB_  ] = (DescIn)[IMB_  ];    \
              (DescOut)[INB_  ] = (DescIn)[INB_  ];    \
              (DescOut)[MB_   ] = (DescIn)[MB_   ];    \
              (DescOut)[NB_   ] = (DescIn)[NB_   ];    \
              (DescOut)[RSRC_ ] = (DescIn)[RSRC_ ];    \
              (DescOut)[CSRC_ ] = (DescIn)[CSRC_ ];    \
              (DescOut)[CTXT_ ] = (DescIn)[CTXT_ ];    \
              (DescOut)[LLD_  ] = (DescIn)[LLD_  ];    \
           }

#define    MDescTrans(DescIn, DescOut) \
           { \
              if ( (DescIn)[DTYPE_] == BLOCK_CYCLIC_2D ) \
              { \
                 (DescOut)[DTYPE_] = BLOCK_CYCLIC_2D_INB; \
                 (DescOut)[M_    ] = (DescIn)[M1_    ];   \
                 (DescOut)[N_    ] = (DescIn)[N1_    ];   \
                 (DescOut)[IMB_  ] = (DescIn)[MB1_   ];   \
                 (DescOut)[INB_  ] = (DescIn)[NB1_   ];   \
                 (DescOut)[MB_   ] = (DescIn)[MB1_   ];   \
                 (DescOut)[NB_   ] = (DescIn)[NB1_   ];   \
                 (DescOut)[RSRC_ ] = (DescIn)[RSRC1_ ];   \
                 (DescOut)[CSRC_ ] = (DescIn)[CSRC1_ ];   \
                 (DescOut)[CTXT_ ] = (DescIn)[CTXT1_ ];   \
                 (DescOut)[LLD_  ] = (DescIn)[LLD1_  ];   \
              } \
              else if ( (DescIn)[DTYPE_] == BLOCK_CYCLIC_2D_INB ) \
              { \
                 (DescOut)[DTYPE_] = BLOCK_CYCLIC_2D_INB; \
                 (DescOut)[M_    ] = (DescIn)[M_    ];    \
                 (DescOut)[N_    ] = (DescIn)[N_    ];    \
                 (DescOut)[IMB_  ] = (DescIn)[IMB_  ];    \
                 (DescOut)[INB_  ] = (DescIn)[INB_  ];    \
                 (DescOut)[MB_   ] = (DescIn)[MB_   ];    \
                 (DescOut)[NB_   ] = (DescIn)[NB_   ];    \
                 (DescOut)[RSRC_ ] = (DescIn)[RSRC_ ];    \
                 (DescOut)[CSRC_ ] = (DescIn)[CSRC_ ];    \
                 (DescOut)[CTXT_ ] = (DescIn)[CTXT_ ];    \
                 (DescOut)[LLD_  ] = (DescIn)[LLD_  ];    \
              } \
              else \
              { \
                 (DescOut)[DTYPE_] = (DescIn)[0]; \
                 (DescOut)[CTXT_ ] = (DescIn)[1]; \
                 (DescOut)[M_    ] = 0;           \
                 (DescOut)[N_    ] = 0;           \
                 (DescOut)[IMB_  ] = 1;           \
                 (DescOut)[INB_  ] = 1;           \
                 (DescOut)[MB_   ] = 1;           \
                 (DescOut)[NB_   ] = 1;           \
                 (DescOut)[RSRC_ ] = 0;           \
                 (DescOut)[CSRC_ ] = 0;           \
                 (DescOut)[LLD_  ] = 1;           \
              } \
           }

#define    MIndxTrans( I, J, i, j ) \
           { \
              i = *I - 1; \
              j = *J - 1; \
           }

#if( _F2C_CALL_ == _F2C_ADD_ )
/*
*  These defines  set  up  the  naming scheme required to have a FORTRAN
*  routine called by a C routine. No redefinition is necessary  to  have
*  the following FORTRAN to C interface:
*
*           FORTRAN DECLARATION            C CALL
*           SUBROUTINE PDFOO(...)          pdfoo_(...)
*
*  This is the PBLAS default.
*/

#endif

#if( _F2C_CALL_ == _F2C_F77ISF2C )
/*
*  These defines  set  up  the  naming scheme required to have a FORTRAN
*  routine called by a C routine for systems where  the FORTRAN compiler
*  is actually f2c (a FORTRAN to C conversion utility).
*
*           FORTRAN DECLARATION            C CALL
*           SUBROUTINE PDFOO(...)          pdfoo__(...)
*/

#endif

#if( _F2C_CALL_ == _F2C_UPCASE )
/*
*  These defines  set  up  the  naming scheme required to have a FORTRAN
*  routine called by a C routine with the following  FORTRAN to C inter-
*  face:
*
*           FORTRAN DECLARATION            C CALL
*           SUBROUTINE PDFOO(...)          PDFOO(...)
*/
#define    immadd_             IMMADD
#define    smmadd_             SMMADD
#define    dmmadd_             DMMADD
#define    cmmadd_             CMMADD
#define    zmmadd_             ZMMADD

#define    immtadd_            IMMTADD
#define    smmtadd_            SMMTADD
#define    dmmtadd_            DMMTADD
#define    cmmtadd_            CMMTADD
#define    zmmtadd_            ZMMTADD

#define    smmcadd_            SMMCADD
#define    dmmcadd_            DMMCADD
#define    cmmcadd_            CMMCADD
#define    zmmcadd_            ZMMCADD

#define    smmtcadd_           SMMTCADD
#define    dmmtcadd_           DMMTCADD
#define    cmmtcadd_           CMMTCADD
#define    zmmtcadd_           ZMMTCADD

#define    immdda_             IMMDDA
#define    smmdda_             SMMDDA
#define    dmmdda_             DMMDDA
#define    cmmdda_             CMMDDA
#define    zmmdda_             ZMMDDA

#define    smmddac_            SMMDDAC
#define    dmmddac_            DMMDDAC
#define    cmmddac_            CMMDDAC
#define    zmmddac_            ZMMDDAC

#define    immddat_            IMMDDAT
#define    smmddat_            SMMDDAT
#define    dmmddat_            DMMDDAT
#define    cmmddat_            CMMDDAT
#define    zmmddat_            ZMMDDAT

#define    smmddact_           SMMDDACT
#define    dmmddact_           DMMDDACT
#define    cmmddact_           CMMDDACT
#define    zmmddact_           ZMMDDACT

#define    sasqrtb_            SASQRTB
#define    dasqrtb_            DASQRTB

#define    sset_               SSET
#define    dset_               DSET
#define    cset_               CSET
#define    zset_               ZSET

#define    svasum_             SVASUM
#define    dvasum_             DVASUM
#define    scvasum_            SCVASUM
#define    dzvasum_            DZVASUM

#define    sascal_             SASCAL
#define    dascal_             DASCAL

#define    scshft_             SCSHFT
#define    dcshft_             DCSHFT
#define    ccshft_             CCSHFT
#define    zcshft_             ZCSHFT

#define    srshft_             SRSHFT
#define    drshft_             DRSHFT
#define    crshft_             CRSHFT
#define    zrshft_             ZRSHFT

#define    svvdot_             SVVDOT
#define    dvvdot_             DVVDOT
#define    cvvdotc_            CVVDOTC
#define    cvvdotu_            CVVDOTU
#define    zvvdotc_            ZVVDOTC
#define    zvvdotu_            ZVVDOTU

#define    stzpad_             STZPAD
#define    dtzpad_             DTZPAD
#define    ctzpad_             CTZPAD
#define    ztzpad_             ZTZPAD

#define    stzpadcpy_          STZPADCPY
#define    dtzpadcpy_          DTZPADCPY
#define    ctzpadcpy_          CTZPADCPY
#define    ztzpadcpy_          ZTZPADCPY

#define    stzscal_            STZSCAL
#define    dtzscal_            DTZSCAL
#define    ctzscal_            CTZSCAL
#define    ztzscal_            ZTZSCAL

#define    chescal_            CHESCAL
#define    zhescal_            ZHESCAL

#define    ctzcnjg_            CTZCNJG
#define    ztzcnjg_            ZTZCNJG

#define    sagemv_             SAGEMV
#define    dagemv_             DAGEMV
#define    cagemv_             CAGEMV
#define    zagemv_             ZAGEMV

#define    sasymv_             SASYMV
#define    dasymv_             DASYMV
#define    casymv_             CASYMV
#define    zasymv_             ZASYMV
#define    cahemv_             CAHEMV
#define    zahemv_             ZAHEMV

#define    satrmv_             SATRMV
#define    datrmv_             DATRMV
#define    catrmv_             CATRMV
#define    zatrmv_             ZATRMV

#define    csymv_              CSYMV
#define    zsymv_              ZSYMV

#define    csyr_               CSYR
#define    zsyr_               ZSYR

#define    csyr2_              CSYR2
#define    zsyr2_              ZSYR2

#endif

#if( _F2C_CALL_ == _F2C_NOCHANGE )
/*
*  These defines  set  up  the  naming scheme required to have a FORTRAN
*  routine called by a C routine with the following  FORTRAN to C inter-
*  face:
*
*           FORTRAN DECLARATION            C CALL
*           SUBROUTINE PDFOO(...)          pdfoo(...)
*/
#define    immadd_             immadd
#define    smmadd_             smmadd
#define    dmmadd_             dmmadd
#define    cmmadd_             cmmadd
#define    zmmadd_             zmmadd

#define    immtadd_            immtadd
#define    smmtadd_            smmtadd
#define    dmmtadd_            dmmtadd
#define    cmmtadd_            cmmtadd
#define    zmmtadd_            zmmtadd

#define    smmcadd_            smmcadd
#define    dmmcadd_            dmmcadd
#define    cmmcadd_            cmmcadd
#define    zmmcadd_            zmmcadd

#define    smmtcadd_           smmtcadd
#define    dmmtcadd_           dmmtcadd
#define    cmmtcadd_           cmmtcadd
#define    zmmtcadd_           zmmtcadd

#define    immdda_             immdda
#define    smmdda_             smmdda
#define    dmmdda_             dmmdda
#define    cmmdda_             cmmdda
#define    zmmdda_             zmmdda

#define    smmddac_            smmddac
#define    dmmddac_            dmmddac
#define    cmmddac_            cmmddac
#define    zmmddac_            zmmddac

#define    immddat_            immddat
#define    smmddat_            smmddat
#define    dmmddat_            dmmddat
#define    cmmddat_            cmmddat
#define    zmmddat_            zmmddat

#define    smmddact_           smmddact
#define    dmmddact_           dmmddact
#define    cmmddact_           cmmddact
#define    zmmddact_           zmmddact

#define    sasqrtb_            sasqrtb
#define    dasqrtb_            dasqrtb

#define    sset_               sset
#define    dset_               dset
#define    cset_               cset
#define    zset_               zset

#define    svasum_             svasum
#define    dvasum_             dvasum
#define    scvasum_            scvasum
#define    dzvasum_            dzvasum

#define    sascal_             sascal
#define    dascal_             dascal

#define    scshft_             scshft
#define    dcshft_             dcshft
#define    ccshft_             ccshft
#define    zcshft_             zcshft

#define    srshft_             srshft
#define    drshft_             drshft
#define    crshft_             crshft
#define    zrshft_             zrshft

#define    svvdot_             svvdot
#define    dvvdot_             dvvdot
#define    cvvdotc_            cvvdotc
#define    cvvdotu_            cvvdotu
#define    zvvdotc_            zvvdotc
#define    zvvdotu_            zvvdotu

#define    stzpad_             stzpad
#define    dtzpad_             dtzpad
#define    ctzpad_             ctzpad
#define    ztzpad_             ztzpad

#define    stzpadcpy_          stzpadcpy
#define    dtzpadcpy_          dtzpadcpy
#define    ctzpadcpy_          ctzpadcpy
#define    ztzpadcpy_          ztzpadcpy

#define    stzscal_            stzscal
#define    dtzscal_            dtzscal
#define    ctzscal_            ctzscal
#define    ztzscal_            ztzscal

#define    chescal_            chescal
#define    zhescal_            zhescal

#define    ctzcnjg_            ctzcnjg
#define    ztzcnjg_            ztzcnjg

#define    sagemv_             sagemv
#define    dagemv_             dagemv
#define    cagemv_             cagemv
#define    zagemv_             zagemv

#define    sasymv_             sasymv
#define    dasymv_             dasymv
#define    casymv_             casymv
#define    zasymv_             zasymv
#define    cahemv_             cahemv
#define    zahemv_             zahemv

#define    satrmv_             satrmv
#define    datrmv_             datrmv
#define    catrmv_             catrmv
#define    zatrmv_             zatrmv

#define    csymv_              csymv
#define    zsymv_              zsymv

#define    csyr_               csyr
#define    zsyr_               zsyr

#define    csyr2_              csyr2
#define    zsyr2_              zsyr2

#endif
/*
*  ---------------------------------------------------------------------
*  Function prototypes
*  ---------------------------------------------------------------------
*/
#ifdef __STDC__

F_VOID_FCT     immadd_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     smmadd_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     dmmadd_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     cmmadd_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     zmmadd_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );

F_VOID_FCT     smmcadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     dmmcadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     cmmcadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     zmmcadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );

F_VOID_FCT     immtadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     smmtadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     dmmtadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     cmmtadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     zmmtadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );

F_VOID_FCT     smmtcadd_       ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     dmmtcadd_       ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     cmmtcadd_       ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     zmmtcadd_       ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );

F_VOID_FCT     immdda_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     smmdda_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     dmmdda_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     cmmdda_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     zmmdda_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );

F_VOID_FCT     smmddac_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     dmmddac_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     cmmddac_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     zmmddac_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );

F_VOID_FCT     immddat_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     smmddat_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     dmmddat_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     cmmddat_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     zmmddat_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );

F_VOID_FCT     smmddact_       ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     dmmddact_       ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     cmmddact_       ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     zmmddact_       ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );

F_VOID_FCT     sasqrtb_        ( float *,   float *,   float * );
F_VOID_FCT     dasqrtb_        ( double *,  double *,  double * );

F_VOID_FCT     sset_           ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     dset_           ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     cset_           ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     zset_           ( Int *,     char *,    char *,
                                 Int * );

F_VOID_FCT     svasum_         ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     dvasum_         ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     scvasum_        ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     dzvasum_        ( Int *,     char *,    char *,
                                 Int * );

F_VOID_FCT     sascal_         ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     dascal_         ( Int *,     char *,    char *,
                                 Int * );

F_VOID_FCT     scshft_         ( Int *,     Int *,     Int *,
                                 char *,    Int * );
F_VOID_FCT     dcshft_         ( Int *,     Int *,     Int *,
                                 char *,    Int * );
F_VOID_FCT     ccshft_         ( Int *,     Int *,     Int *,
                                 char *,    Int * );
F_VOID_FCT     zcshft_         ( Int *,     Int *,     Int *,
                                 char *,    Int * );

F_VOID_FCT     srshft_         ( Int *,     Int *,     Int *,
                                 char *,    Int * );
F_VOID_FCT     drshft_         ( Int *,     Int *,     Int *,
                                 char *,    Int * );
F_VOID_FCT     crshft_         ( Int *,     Int *,     Int *,
                                 char *,    Int * );
F_VOID_FCT     zrshft_         ( Int *,     Int *,     Int *,
                                 char *,    Int * );

F_VOID_FCT     svvdot_         ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     dvvdot_         ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     cvvdotu_        ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     cvvdotc_        ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     zvvdotu_        ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     zvvdotc_        ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );

F_VOID_FCT     stzpad_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 char *,    char *,    Int * );
F_VOID_FCT     dtzpad_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 char *,    char *,    Int * );
F_VOID_FCT     ctzpad_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 char *,    char *,    Int * );
F_VOID_FCT     ztzpad_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 char *,    char *,    Int * );

F_VOID_FCT     stzpadcpy_      ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     dtzpadcpy_      ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     ctzpadcpy_      ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     ztzpadcpy_      ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 Int *,     char *,    Int * );

F_VOID_FCT     stzscal_        ( F_CHAR_T,  Int *,     Int *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     dtzscal_        ( F_CHAR_T,  Int *,     Int *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     ctzscal_        ( F_CHAR_T,  Int *,     Int *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     ztzscal_        ( F_CHAR_T,  Int *,     Int *,
                                 Int *,     char *,    char *,
                                 Int * );

F_VOID_FCT     chescal_        ( F_CHAR_T,  Int *,     Int *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     zhescal_        ( F_CHAR_T,  Int *,     Int *,
                                 Int *,     char *,    char *,
                                 Int * );

F_VOID_FCT     ctzcnjg_        ( F_CHAR_T,  Int *,     Int *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     ztzcnjg_        ( F_CHAR_T,  Int *,     Int *,
                                 Int *,     char *,    char *,
                                 Int * );

F_VOID_FCT     sagemv_         ( F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     dagemv_         ( F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     cagemv_         ( F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     zagemv_         ( F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );

F_VOID_FCT     sasymv_         ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     dasymv_         ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     casymv_         ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     zasymv_         ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     cahemv_         ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     zahemv_         ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );

F_VOID_FCT     satrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     datrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     catrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     zatrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );

F_VOID_FCT     csymv_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     zsymv_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );

F_VOID_FCT     csyr_           ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int * );
F_VOID_FCT     zsyr_           ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int * );

F_VOID_FCT     csyr2_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     zsyr2_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );

void           PB_Ctzsyr       ( PBTYP_T *, char *,    Int,
                                 Int,       Int,       Int,
                                 char *,    char *,    Int,
                                 char *,    Int,       char *,
                                 Int );
void           PB_Ctzher       ( PBTYP_T *, char *,    Int,
                                 Int,       Int,       Int,
                                 char *,    char *,    Int,
                                 char *,    Int,       char *,
                                 Int );
void           PB_Ctzsyr2      ( PBTYP_T *, char *,    Int,
                                 Int,       Int,       Int,
                                 char *,    char *,    Int,
                                 char *,    Int,       char *,
                                 Int,       char *,    Int,
                                 char *,    Int );
void           PB_Ctzher2      ( PBTYP_T *, char *,    Int,
                                 Int,       Int,       Int,
                                 char *,    char *,    Int,
                                 char *,    Int,       char *,
                                 Int,       char *,    Int,
                                 char *,    Int );
void           PB_Ctztrmv      ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       Int,       Int,
                                 char *,    char *,    Int,
                                 char *,    Int,       char *,
                                 Int );
void           PB_Ctzatrmv     ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       Int,       Int,
                                 char *,    char *,    Int,
                                 char *,    Int,       char *,
                                 Int );
void           PB_Ctzsymv      ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       Int,
                                 Int,       char *,    char *,
                                 Int,       char *,    Int,
                                 char *,    Int,       char *,
                                 Int,       char *,    Int );
void           PB_Ctzhemv      ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       Int,
                                 Int,       char *,    char *,
                                 Int,       char *,    Int,
                                 char *,    Int,       char *,
                                 Int,       char *,    Int );
void           PB_Ctzasymv     ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       Int,
                                 Int,       char *,    char *,
                                 Int,       char *,    Int,
                                 char *,    Int,       char *,
                                 Int,       char *,    Int );
void           PB_Ctzahemv     ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       Int,
                                 Int,       char *,    char *,
                                 Int,       char *,    Int,
                                 char *,    Int,       char *,
                                 Int,       char *,    Int );

void           PB_Ctzsyrk      ( PBTYP_T *, char *,    Int,
                                 Int,       Int,       Int,
                                 char *,    char *,    Int,
                                 char *,    Int,       char *,
                                 Int );
void           PB_Ctzherk      ( PBTYP_T *, char *,    Int,
                                 Int,       Int,       Int,
                                 char *,    char *,    Int,
                                 char *,    Int,       char *,
                                 Int );
void           PB_Ctzsyr2k     ( PBTYP_T *, char *,    Int,
                                 Int,       Int,       Int,
                                 char *,    char *,    Int,
                                 char *,    Int,       char *,
                                 Int,       char *,    Int,
                                 char *,    Int );
void           PB_Ctzher2k     ( PBTYP_T *, char *,    Int,
                                 Int,       Int,       Int,
                                 char *,    char *,    Int,
                                 char *,    Int,       char *,
                                 Int,       char *,    Int,
                                 char *,    Int );
void           PB_Ctztrmm      ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       Int,       Int,
                                 char *,    char *,    Int,
                                 char *,    Int,       char *,
                                 Int );
void           PB_Ctzsymm      ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       Int,
                                 Int,       char *,    char *,
                                 Int,       char *,    Int,
                                 char *,    Int,       char *,
                                 Int,       char *,    Int );
void           PB_Ctzhemm      ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       Int,
                                 Int,       char *,    char *,
                                 Int,       char *,    Int,
                                 char *,    Int,       char *,
                                 Int,       char *,    Int );

void           PB_CpswapNN     ( PBTYP_T *, Int,       char *,
                                 Int,       Int,       Int *,
                                 Int,       char *,    Int,
                                 Int,       Int *,     Int );
void           PB_CpswapND     ( PBTYP_T *, Int,       char *,
                                 Int,       Int,       Int *,
                                 Int,       char *,    Int,
                                 Int,       Int *,     Int );
void           PB_Cpdot11      ( PBTYP_T *, Int,       char *,
                                 char *,    Int,       Int,
                                 Int *,     Int,       char *,
                                 Int,       Int,       Int *,
                                 Int,       VVDOT_T );
void           PB_CpdotNN      ( PBTYP_T *, Int,       char *,
                                 char *,    Int,       Int,
                                 Int *,     Int,       char *,
                                 Int,       Int,       Int *,
                                 Int,       VVDOT_T );
void           PB_CpdotND      ( PBTYP_T *, Int,       char *,
                                 char *,    Int,       Int,
                                 Int *,     Int,       char *,
                                 Int,       Int,       Int *,
                                 Int,       VVDOT_T );
void           PB_CpaxpbyNN    ( PBTYP_T *, char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       Int,       Int *,
                                 char *,    char *,    char *,
                                 Int,       Int,       Int *,
                                 char * );
void           PB_CpaxpbyND    ( PBTYP_T *, char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       Int,       Int *,
                                 char *,    char *,    char *,
                                 Int,       Int,       Int *,
                                 char * );
void           PB_CpaxpbyDN    ( PBTYP_T *, char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       Int,       Int *,
                                 char *,    char *,    char *,
                                 Int,       Int,       Int *,
                                 char * );
void           PB_Cpaxpby      ( PBTYP_T *, char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       Int,       Int *,
                                 char *,    char *,    char *,
                                 Int,       Int,       Int *,
                                 char * );

void           PB_Cpsyr        ( PBTYP_T *, char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       char *,    Int,
                                 char *,    Int,       Int,
                                 Int *,     TZSYR_T );
void           PB_Cpsyr2       ( PBTYP_T *, char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       char *,    Int,
                                 char *,    Int,       char *,
                                 Int,       char *,    Int,
                                 Int,       Int *,     TZSYR2_T );
void           PB_Cptrm        ( PBTYP_T *, PBTYP_T *, char *,
                                 char *,    char *,    char *,
                                 Int,       Int,       char *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    Int,
                                 char *,    Int,       TZTRM_T );
void           PB_Cpsym        ( PBTYP_T *, PBTYP_T *, char *,
                                 char *,    Int,       Int,
                                 char *,    char *,    Int,
                                 Int,       Int *,     char *,
                                 Int,       char *,    Int,
                                 char *,    Int,       char *,
                                 Int,       TZSYM_T );
void           PB_Cpgeadd      ( PBTYP_T *, char *,    char *,
                                 char *,    Int,       Int,
                                 char *,    char *,    Int,
                                 Int,       Int *,     char *,
                                 char *,    Int,       Int,
                                 Int * );
void           PB_Cptradd      ( PBTYP_T *, char *,    char *,
                                 char *,    Int,       Int,
                                 char *,    char *,    Int,
                                 Int,       Int *,     char *,
                                 char *,    Int,       Int,
                                 Int * );
void           PB_Cptran       ( PBTYP_T *, char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       Int,       Int *,
                                 char *,    char *,    Int,
                                 Int,       Int * );
void           PB_Cptrsv       ( PBTYP_T *, Int,       char *,
                                 char *,    char *,    Int,
                                 char *,    Int,       Int,
                                 Int *,     char *,    Int,
                                 char *,    Int );
void           PB_Cptrsm       ( PBTYP_T *, Int,       char *,
                                 char *,    char *,    char *,
                                 Int,       Int,       char *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    Int,
                                 char *,    Int );

void           PB_CpgemmAB     ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       Int,       char *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    Int,
                                 Int,       Int *,     char *,
                                 char *,    Int,       Int,
                                 Int * );
void           PB_CpgemmAC     ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       Int,       char *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    Int,
                                 Int,       Int *,     char *,
                                 char *,    Int,       Int,
                                 Int * );
void           PB_CpgemmBC     ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       Int,       char *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    Int,
                                 Int,       Int *,     char *,
                                 char *,    Int,       Int,
                                 Int * );
void           PB_CpsymmAB     ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       Int,       Int *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    char *,
                                 Int,       Int,       Int * );
void           PB_CpsymmBC     ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       Int,       Int *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    char *,
                                 Int,       Int,       Int * );
void           PB_CpsyrkA      ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       Int,       Int *,
                                 char *,    char *,    Int,
                                 Int,       Int * );
void           PB_CpsyrkAC     ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       Int,       Int *,
                                 char *,    char *,    Int,
                                 Int,       Int * );
void           PB_Cpsyr2kA     ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       Int,       Int *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    char *,
                                 Int,       Int,       Int * );
void           PB_Cpsyr2kAC    ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       Int,       Int *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    char *,
                                 Int,       Int,       Int * );
void           PB_CptrmmAB     ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    char *,
                                 Int,       Int,       char *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    Int,
                                 Int,       Int * );
void           PB_CptrmmB      ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    char *,
                                 Int,       Int,       char *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    Int,
                                 Int,       Int * );
void           PB_CptrsmAB     ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    char *,
                                 Int,       Int,       char *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    Int,
                                 Int,       Int * );
void           PB_CptrsmAB0    ( PBTYP_T *, char *,    char *,
                                 char *,    Int,       Int,
                                 char *,    char *,    Int,
                                 Int,       Int *,     char *,
                                 Int,       Int,       Int *,
                                 char * *,  Int *,     Int * );
void           PB_CptrsmAB1    ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       Int,       Int *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    Int * );
void           PB_CptrsmB      ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    char *,
                                 Int,       Int,       char *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    Int,
                                 Int,       Int * );
#else

F_VOID_FCT     immadd_         ();
F_VOID_FCT     smmadd_         ();
F_VOID_FCT     dmmadd_         ();
F_VOID_FCT     cmmadd_         ();
F_VOID_FCT     zmmadd_         ();

F_VOID_FCT     smmcadd_        ();
F_VOID_FCT     dmmcadd_        ();
F_VOID_FCT     cmmcadd_        ();
F_VOID_FCT     zmmcadd_        ();

F_VOID_FCT     immtadd_        ();
F_VOID_FCT     smmtadd_        ();
F_VOID_FCT     dmmtadd_        ();
F_VOID_FCT     cmmtadd_        ();
F_VOID_FCT     zmmtadd_        ();

F_VOID_FCT     smmtcadd_       ();
F_VOID_FCT     dmmtcadd_       ();
F_VOID_FCT     cmmtcadd_       ();
F_VOID_FCT     zmmtcadd_       ();

F_VOID_FCT     immdda_         ();
F_VOID_FCT     smmdda_         ();
F_VOID_FCT     dmmdda_         ();
F_VOID_FCT     cmmdda_         ();
F_VOID_FCT     zmmdda_         ();

F_VOID_FCT     smmddac_        ();
F_VOID_FCT     dmmddac_        ();
F_VOID_FCT     cmmddac_        ();
F_VOID_FCT     zmmddac_        ();

F_VOID_FCT     immddat_        ();
F_VOID_FCT     smmddat_        ();
F_VOID_FCT     dmmddat_        ();
F_VOID_FCT     cmmddat_        ();
F_VOID_FCT     zmmddat_        ();

F_VOID_FCT     smmddact_       ();
F_VOID_FCT     dmmddact_       ();
F_VOID_FCT     cmmddact_       ();
F_VOID_FCT     zmmddact_       ();

F_VOID_FCT     sasqrtb_        ();
F_VOID_FCT     dasqrtb_        ();

F_VOID_FCT     sset_           ();
F_VOID_FCT     dset_           ();
F_VOID_FCT     cset_           ();
F_VOID_FCT     zset_           ();

F_VOID_FCT     svasum_         ();
F_VOID_FCT     dvasum_         ();
F_VOID_FCT     scvasum_        ();
F_VOID_FCT     dzvasum_        ();

F_VOID_FCT     sascal_         ();
F_VOID_FCT     dascal_         ();

F_VOID_FCT     scshft_         ();
F_VOID_FCT     dcshft_         ();
F_VOID_FCT     ccshft_         ();
F_VOID_FCT     zcshft_         ();

F_VOID_FCT     srshft_         ();
F_VOID_FCT     drshft_         ();
F_VOID_FCT     crshft_         ();
F_VOID_FCT     zrshft_         ();

F_VOID_FCT     svvdot_         ();
F_VOID_FCT     dvvdot_         ();
F_VOID_FCT     cvvdotc_        ();
F_VOID_FCT     cvvdotu_        ();
F_VOID_FCT     zvvdotc_        ();
F_VOID_FCT     zvvdotu_        ();

F_VOID_FCT     stzpad_         ();
F_VOID_FCT     dtzpad_         ();
F_VOID_FCT     ctzpad_         ();
F_VOID_FCT     ztzpad_         ();

F_VOID_FCT     stzpadcpy_      ();
F_VOID_FCT     dtzpadcpy_      ();
F_VOID_FCT     ctzpadcpy_      ();
F_VOID_FCT     ztzpadcpy_      ();

F_VOID_FCT     stzscal_        ();
F_VOID_FCT     dtzscal_        ();
F_VOID_FCT     ctzscal_        ();
F_VOID_FCT     ztzscal_        ();

F_VOID_FCT     chescal_        ();
F_VOID_FCT     zhescal_        ();

F_VOID_FCT     ctzcnjg_        ();
F_VOID_FCT     ztzcnjg_        ();

F_VOID_FCT     sagemv_         ();
F_VOID_FCT     dagemv_         ();
F_VOID_FCT     cagemv_         ();
F_VOID_FCT     zagemv_         ();

F_VOID_FCT     sasymv_         ();
F_VOID_FCT     dasymv_         ();
F_VOID_FCT     casymv_         ();
F_VOID_FCT     zasymv_         ();
F_VOID_FCT     cahemv_         ();
F_VOID_FCT     zahemv_         ();

F_VOID_FCT     satrmv_         ();
F_VOID_FCT     datrmv_         ();
F_VOID_FCT     catrmv_         ();
F_VOID_FCT     zatrmv_         ();

F_VOID_FCT     csymv_          ();
F_VOID_FCT     zsymv_          ();

F_VOID_FCT     csyr_           ();
F_VOID_FCT     zsyr_           ();

F_VOID_FCT     csyr2_          ();
F_VOID_FCT     zsyr2_          ();

void           PB_Ctzsyr       ();
void           PB_Ctzher       ();
void           PB_Ctzsyr2      ();
void           PB_Ctzher2      ();
void           PB_Ctztrmv      ();
void           PB_Ctzatrmv     ();
void           PB_Ctzsymv      ();
void           PB_Ctzhemv      ();
void           PB_Ctzasymv     ();
void           PB_Ctzahemv     ();
void           PB_Ctzsyrk      ();
void           PB_Ctzherk      ();
void           PB_Ctzsyr2k     ();
void           PB_Ctzher2k     ();
void           PB_Ctztrmm      ();
void           PB_Ctzsymm      ();
void           PB_Ctzhemm      ();

void           PB_CpswapNN     ();
void           PB_CpswapND     ();
void           PB_Cpdot11      ();
void           PB_CpdotNN      ();
void           PB_CpdotND      ();
void           PB_CpaxpbyNN    ();
void           PB_CpaxpbyND    ();
void           PB_CpaxpbyDN    ();
void           PB_Cpaxpby      ();

void           PB_Cpsyr        ();
void           PB_Cpsyr2       ();
void           PB_Cptrm        ();
void           PB_Cpsym        ();
void           PB_Cpgeadd      ();
void           PB_Cptradd      ();
void           PB_Cptran       ();
void           PB_Cptrsv       ();
void           PB_Cptrsm       ();

void           PB_CpgemmAB     ();
void           PB_CpgemmAC     ();
void           PB_CpgemmBC     ();
void           PB_CpsymmAB     ();
void           PB_CpsymmBC     ();
void           PB_CpsyrkA      ();
void           PB_CpsyrkAC     ();
void           PB_Cpsyr2kA     ();
void           PB_Cpsyr2kAC    ();
void           PB_CptrmmAB     ();
void           PB_CptrmmB      ();
void           PB_CptrsmAB     ();
void           PB_CptrsmAB0    ();
void           PB_CptrsmAB1    ();
void           PB_CptrsmB      ();

#endif
                                                             /* TOOLS */
#ifdef __STDC__

Int            PB_Cgcd         ( Int,       Int );
Int            PB_Clcm         ( Int,       Int );

void           PB_Cdescset     ( Int *,     Int,       Int,
                                 Int,       Int,       Int,
                                 Int,       Int,       Int,
                                 Int,       Int );
void           PB_Cdescribe    ( Int,       Int,       Int,
                                 Int,       Int *,     Int,
                                 Int,       Int,       Int,
                                 Int *,     Int *,     Int *,
                                 Int *,     Int *,     Int *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           PB_CargFtoC     ( Int,       Int,       Int *,
                                 Int *,     Int *,     Int * );
Int            PB_Cfirstnb     ( Int,       Int,       Int,
                                 Int );
Int            PB_Clastnb      ( Int,       Int,       Int,
                                 Int );
Int            PB_Cspan        ( Int,       Int,       Int,
                                 Int,       Int,       Int );

void           PB_Cainfog2l    ( Int,       Int,       Int,
                                 Int,       Int *,     Int,
                                 Int,       Int,       Int,
                                 Int *,     Int *,     Int *,
                                 Int *,     Int *,     Int *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           PB_Cinfog2l     ( Int,       Int,       Int *,
                                 Int,       Int,       Int,
                                 Int,       Int *,     Int *,
                                 Int *,     Int * );
Int            PB_Cg2lrem      ( Int,       Int,       Int,
                                 Int,       Int,       Int );
Int            PB_Cindxg2p     ( Int,       Int,       Int,
                                 Int,       Int,       Int );
Int            PB_Cnumroc      ( Int,       Int,       Int,
                                 Int,       Int,       Int,
                                 Int );
Int            PB_Cnpreroc     ( Int,       Int,       Int,
                                 Int,       Int,       Int,
                                 Int );
Int            PB_Cnnxtroc     ( Int,       Int,       Int,
                                 Int,       Int,       Int,
                                 Int );

void           PB_Cconjg       ( PBTYP_T *, char *,    char * );


void           PB_Cwarn        ( Int,       Int,       char *,
                                 char *,    ... );
void           PB_Cabort       ( Int,       char *,    Int );
void           PB_Cchkmat      ( Int,       char *,    char *,
                                 Int,       Int,       Int,
                                 Int,       Int,       Int,
                                 Int *,     Int,       Int * );
void           PB_Cchkvec      ( Int,       char *,    char *,
                                 Int,       Int ,       Int ,
                                 Int,       Int *,     Int,
                                 Int,       Int * );

char *         PB_Cmalloc      ( Int );
char *         PB_Cgetbuf      ( char *,    Int );

PBTYP_T *      PB_Citypeset    ( void );
PBTYP_T *      PB_Cstypeset    ( void );
PBTYP_T *      PB_Cdtypeset    ( void );
PBTYP_T *      PB_Cctypeset    ( void );
PBTYP_T *      PB_Cztypeset    ( void );

Int            pilaenv_        ( Int *,     F_CHAR_T );
char *         PB_Ctop         ( Int *,     char *,    char *,
                                 char * );

void           PB_CVMinit      ( PB_VM_T *, Int,       Int,
                                 Int,       Int,       Int,
                                 Int,       Int,       Int,
                                 Int,       Int,       Int,
                                 Int );
Int            PB_CVMnpq       ( PB_VM_T * );
void           PB_CVMcontig    ( PB_VM_T *, Int *,     Int *,
                                 Int *,     Int * );
Int            PB_CVMloc       ( PBTYP_T *, PB_VM_T *, char *,
                                 char *,    char *,    char *,
                                 Int,       Int,       char *,
                                 char *,    Int,       char *,
                                 char *,    Int );
Int            PB_CVMswp       ( PBTYP_T *, PB_VM_T *, char *,
                                 char *,    char *,    Int,
                                 char *,    Int,       char *,
                                 Int );
Int            PB_CVMpack      ( PBTYP_T *, PB_VM_T *, char *,
                                 char *,    char *,    char *,
                                 Int,       Int,       char *,
                                 char *,    Int,       char *,
                                 char *,    Int );
void           PB_CVMupdate    ( PB_VM_T *, Int,       Int *,
                                 Int * );

void           PB_Cbinfo       ( Int,       Int,       Int,
                                 Int,       Int,       Int,
                                 Int,       Int,       Int,
                                 Int *,     Int *,     Int *,
                                 Int *,     Int *,     Int *,
                                 Int *,     Int *,     Int *,
                                 Int *,     Int * );

void           PB_Cplaprnt     ( PBTYP_T *, Int,       Int,
                                 char *,    Int,       Int,
                                 Int *,     Int,       Int,
                                 char * );
void           PB_Cplaprn2     ( PBTYP_T *, Int,       Int,
                                 char *,    Int,       Int,
                                 Int *,     Int,       Int,
                                 char *,    Int,       Int );
void           PB_Cprnt        ( char,      Int,       Int,
                                 Int,       char *,    Int,
                                 Int,       char * );

void           PB_Cplapad      ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       char *,
                                 char *,    char *,    Int,
                                 Int,       Int * );
void           PB_Cplapd2      ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       char *,
                                 char *,    char *,    Int,
                                 Int,       Int * );
void           PB_Cplascal     ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       char *,
                                 char *,    Int,       Int,
                                 Int * );
void           PB_Cplasca2     ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       char *,
                                 char *,    Int,       Int,
                                 Int * );
void           PB_Cplacnjg     ( PBTYP_T *, Int,       Int,
                                 char *,    char *,    Int,
                                 Int,       Int * );

void           PB_CInV         ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       Int *,
                                 Int,       char *,    Int,
                                 Int,       Int *,     char *,
                                 char * *,  Int *,     Int * );
void           PB_CInV2        ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       Int *,
                                 Int,       char *,    Int,
                                 Int,       Int *,     char *,
                                 char *,    Int,       Int * );
void           PB_CInOutV      ( PBTYP_T *, char *,    Int,
                                 Int,       Int *,     Int,
                                 char *,    char *,    Int,
                                 Int,       Int *,     char *,
                                 char * *,  char * *,  Int *,
                                 Int *,     Int *,     Int * );
void           PB_CInOutV2     ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       Int,
                                 Int *,     Int,       char *,
                                 Int,       Int,       Int *,
                                 char *,    char * *,  Int *,
                                 Int *,     Int *,     Int * );
void           PB_COutV        ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       Int *,
                                 Int,       char * *,  Int *,
                                 Int *,     Int * );
void           PB_CGatherV     ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int,       Int *,
                                 char *,    char * *,  Int *,
                                 Int * );
void           PB_CScatterV    ( PBTYP_T *, char *,    Int,
                                 Int,       char *,    Int,
                                 Int,       Int *,     char *,
                                 char *,    char *,    Int,
                                 Int,       Int *,     char * );
#else

Int            PB_Cgcd         ();
Int            PB_Clcm         ();

void           PB_Cdescset     ();
void           PB_Cdescribe    ();
void           PB_CargFtoC     ();
Int            PB_Cfirstnb     ();
Int            PB_Clastnb      ();
Int            PB_Cspan        ();

void           PB_Cainfog2l    ();
void           PB_Cinfog2l     ();
Int            PB_Cg2lrem      ();
Int            PB_Cindxg2p     ();
Int            PB_Cnumroc      ();
Int            PB_Cnpreroc     ();
Int            PB_Cnnxtroc     ();

void           PB_Cconjg       ();

void           PB_Cwarn        ();
void           PB_Cabort       ();
void           PB_Cchkmat      ();
void           PB_Cchkvec      ();

char *         PB_Cmalloc      ();
char *         PB_Cgetbuf      ();

PBTYP_T *      PB_Citypeset    ();
PBTYP_T *      PB_Cstypeset    ();
PBTYP_T *      PB_Cdtypeset    ();
PBTYP_T *      PB_Cctypeset    ();
PBTYP_T *      PB_Cztypeset    ();

Int            pilaenv_        ();
char *         PB_Ctop         ();

void           PB_CVMinit      ();
Int            PB_CVMnpq       ();
void           PB_CVMcontig    ();
Int            PB_CVMloc       ();
Int            PB_CVMswp       ();
Int            PB_CVMpack      ();
void           PB_CVMupdate    ();

void           PB_Cbinfo       ();

void           PB_Cplaprnt     ();
void           PB_Cplaprn2     ();
void           PB_Cprnt        ();

void           PB_Cplapad      ();
void           PB_Cplapd2      ();
void           PB_Cplascal     ();
void           PB_Cplasca2     ();
void           PB_Cplacnjg     ();

void           PB_CInV         ();
void           PB_CInV2        ();
void           PB_CInOutV      ();
void           PB_CInOutV2     ();
void           PB_COutV        ();
void           PB_CGatherV     ();
void           PB_CScatterV    ();

#endif
