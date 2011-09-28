      program tctst
      include 'mpif.h'
      integer f77com, wgrp, f77grp, Iam, i, ierr
      integer irank(2)
      external Ccommcheck
      integer  Ccommcheck

      call mpi_init(ierr)
      call mpi_comm_size(MPI_COMM_WORLD, i, ierr)
      call mpi_comm_rank(MPI_COMM_WORLD, Iam, ierr)
      if (i .lt. 2) then
         print*,'Need at least 2 processes to run test, aborting.'
      else
         if (Iam .eq. 0) then
            print*,'If this routine does not complete successfully,'
            print*,'Do _NOT_ set TRANSCOMM = -DCSameF77'
            print*,'  '
            print*,'  '
         end if
*
*        Form context with 2 members
*
         irank(1) = 0
         irank(2) = 1
         call mpi_comm_group(MPI_COMM_WORLD, wgrp, ierr)
         call mpi_group_incl(wgrp, 2, irank, f77grp, ierr)
         call mpi_comm_create(MPI_COMM_WORLD, f77grp, f77com, ierr)
         call mpi_group_free(f77grp, ierr)
   
         i = Ccommcheck(MPI_COMM_WORLD, f77com)
         if (Iam .eq. 0) then
            if (i .eq. 0) then
               print*,'Do _NOT_ set TRANSCOMM = -DCSameF77'
            else
            print*,'Set TRANSCOMM = -DCSameF77'
            end if
         end if

         if (f77grp .ne. MPI_COMM_NULL) call mpi_comm_free(f77com, ierr)
      end if
      call mpi_finalize(ierr)

      stop
      end
