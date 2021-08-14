#include <stdlib.h>

#ifndef Int
#define Int int
#endif

void Creshape( Int context_in, Int major_in, Int* context_out, Int major_out,
                    Int first_proc, Int nprow_new, Int npcol_new )
/* major in, major out represent whether processors go row major (1) or
column major (2) in the input and output grids */
{

   /** called subprograms **/
   void proc_inc();
   void Cblacs_gridinfo();
   Int Cblacs_pnum();
   void Cblacs_get();
   void Cblacs_gridmap();

   /** variables **/
   Int i, j;
   Int nprow_in, npcol_in, myrow_in, mycol_in;
   Int nprocs_new;
   Int myrow_old, mycol_old, myrow_new, mycol_new;
   Int pnum;
   Int *grid_new;

/********** executable statements ************/

   nprocs_new = nprow_new * npcol_new;

   Cblacs_gridinfo( context_in, &nprow_in, &npcol_in, &myrow_in, &mycol_in );

   /* Quick return if possible */
   if( ( nprow_in == nprow_new ) && ( npcol_in == npcol_new ) &&
       ( first_proc == 0 ) && ( major_in == major_out ) )
   {
      *context_out = context_in;
      return;
   }

   /* allocate space for new process mapping */
   grid_new = (Int *) malloc( nprocs_new * sizeof( Int ) );

   /* set place in old grid to start grabbing processors for new grid */
   myrow_old = 0; mycol_old = 0;
   if ( major_in == 1 ) /* row major */
   {
      myrow_old = first_proc / nprow_in;
      mycol_old = first_proc % nprow_in;
   }
   else                  /* col major */
   {
      myrow_old = first_proc % nprow_in;
      mycol_old = first_proc / nprow_in;
   }

   myrow_new = 0; mycol_new = 0;

   /* Set up array of process numbers for new grid */
   for (i=0; i< nprocs_new; i++ )
   {
      pnum = Cblacs_pnum( context_in, myrow_old, mycol_old );
      grid_new[ (mycol_new * nprow_new) + myrow_new ] = pnum;
      proc_inc( &myrow_old, &mycol_old, nprow_in, npcol_in, major_in );
      proc_inc( &myrow_new, &mycol_new, nprow_new, npcol_new, major_out );
   }

   /* get context */
   Cblacs_get( context_in, 10, context_out );

   /* allocate grid */
   Cblacs_gridmap( context_out, grid_new, nprow_new, nprow_new, npcol_new );

   /* free malloced space */
   free( grid_new );
}

/*************************************************************************/
void reshape( Int* context_in, Int* major_in, Int* context_out, Int* major_out,
                    Int* first_proc, Int* nprow_new, Int* npcol_new )
{
   Creshape( *context_in, *major_in, context_out, *major_out,
                    *first_proc, *nprow_new, *npcol_new );
}
/*************************************************************************/
void RESHAPE( Int* context_in, Int* major_in, Int* context_out, Int* major_out,
                    Int* first_proc, Int* nprow_new, Int* npcol_new )
{
   Creshape( *context_in, *major_in, context_out, *major_out,
                    *first_proc, *nprow_new, *npcol_new );
}
/*************************************************************************/
void reshape_( Int* context_in, Int* major_in, Int* context_out, Int* major_out,
                    Int* first_proc, Int* nprow_new, Int* npcol_new )
{
   Creshape( *context_in, *major_in, context_out, *major_out,
                    *first_proc, *nprow_new, *npcol_new );
}
/*************************************************************************/
void proc_inc( Int* myrow, Int* mycol, Int nprow, Int npcol, Int major )
{
   if( major == 1) /* row major */
   {
      if( *mycol == npcol-1 )
      {
         *mycol = 0;
         if( *myrow == nprow-1 )
         {
            *myrow = 0;
         }
         else
         {
            *myrow = *myrow + 1;
         }
      }
      else
      {
         *mycol = *mycol + 1;
      }
   }
   else            /* col major */
   {
      if( *myrow == nprow-1 )
      {
         *myrow = 0;
         if( *mycol == npcol-1 )
         {
            *mycol = 0;
         }
         else
         {
            *mycol = *mycol + 1;
         }
      }
      else
      {
         *myrow = *myrow + 1;
      }
   }
}

