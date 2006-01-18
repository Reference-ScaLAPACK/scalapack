#include <stdlib.h>

void Creshape( context_in, major_in, context_out, major_out,
                    first_proc, nprow_new, npcol_new )
int context_in, *context_out, first_proc, major_in, major_out, nprow_new, npcol_new;
/* major in, major out represent whether processors go row major (1) or
column major (2) in the input and output grids */
{

   /** called subprograms **/
   void proc_inc();
   void Cblacs_gridinfo();
   int Cblacs_pnum();
   void Cblacs_get();
   void Cblacs_gridmap();

   /** variables **/
   int i, j;
   int nprow_in, npcol_in, myrow_in, mycol_in;
   int nprocs_new;
   int myrow_old, mycol_old, myrow_new, mycol_new;
   int pnum;
   int *grid_new;

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
   grid_new = (int *) malloc( nprocs_new * sizeof( int ) );

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
void reshape( context_in, major_in, context_out, major_out,
                    first_proc, nprow_new, npcol_new )
int *context_in, *context_out, *first_proc, *major_in, *major_out, *nprow_new, *npcol_new;
{
   Creshape( *context_in, *major_in, context_out, *major_out,
                    *first_proc, *nprow_new, *npcol_new );
}
/*************************************************************************/
void RESHAPE( context_in, major_in, context_out, major_out,
                    first_proc, nprow_new, npcol_new )
int *context_in, *context_out, *first_proc, *major_in, *major_out, *nprow_new, *npcol_new;
{
   Creshape( *context_in, *major_in, context_out, *major_out,
                    *first_proc, *nprow_new, *npcol_new );
}
/*************************************************************************/
void reshape_( context_in, major_in, context_out, major_out,
                    first_proc, nprow_new, npcol_new )
int *context_in, *context_out, *first_proc, *major_in, *major_out, *nprow_new, *npcol_new;
{
   Creshape( *context_in, *major_in, context_out, *major_out,
                    *first_proc, *nprow_new, *npcol_new );
}
/*************************************************************************/
void proc_inc( myrow, mycol, nprow, npcol, major )
int *myrow, *mycol, nprow, npcol, major;
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

