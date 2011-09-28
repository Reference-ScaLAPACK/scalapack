      program fmpi_sane
*
*     This program checks to make sure that you can run a basic program
*     on your machine using the Fortran77 interface to MPI.
*     Can increase parameter wastesz, if you think size of executable
*     is causing launching problem.
*
      include 'mpif.h'
      integer nproc, wastesz
      parameter (nproc = 4)
      parameter (wastesz = 100)
      integer i, Iam, Np, ierr
      integer mcom, wgrp, mgrp
      integer irank(nproc), stat(MPI_STATUS_SIZE)
      double precision WasteSpc(wastesz)

      call mpi_init(ierr)
      call mpi_comm_size(MPI_COMM_WORLD, Np, ierr)
      if (Np .lt. nproc) then
         print*,'Not enough processes to run sanity check'
         call mpi_abort(MPI_COMM_WORLD, -1, ierr)
      end if
*
*     Access all of WasteSpc
*
      do 10 i = 1, wastesz
         WasteSpc(i) = 0.0D0
10    continue
*
*     Form context with NPROC members
*
      do 20 i = 1, nproc
         irank(i) = i - 1
20    continue
      call mpi_comm_group(MPI_COMM_WORLD, wgrp, ierr)
      call mpi_group_incl(wgrp, nproc, irank, mgrp, ierr)
      call mpi_comm_create(MPI_COMM_WORLD, mgrp, mcom, ierr)
      call mpi_group_free(mgrp, ierr)
*
*     Everyone in new communicator sends a message to his neighbor
*
      if (mcom .ne. MPI_COMM_NULL) then
         call mpi_comm_rank(mcom, Iam, ierr)
*
*        Odd nodes receive first, so we don't hang if MPI_Send is
*        globally blocking
*
         if (mod(Iam, 2) .ne. 0) then
            call mpi_recv(i, 1, MPI_INTEGER, MOD(nproc+Iam-1, nproc), 
     &                    0, mcom, stat, ierr)
            call mpi_send(Iam, 1, MPI_INTEGER, MOD(Iam+1, nproc), 
     &                    0, mcom, ierr)
         else
            call mpi_send(Iam, 1, MPI_INTEGER, MOD(Iam+1, nproc), 
     &                    0, mcom, ierr)
            call mpi_recv(i, 1, MPI_INTEGER, MOD(nproc+Iam-1, nproc), 
     &                    0, mcom, stat, ierr)
         end if
*
*        Make sure we've received the right information
*
         if (i .ne. MOD(nproc+Iam-1, nproc)) then
            print*,'Communication does not seem to work properly!!'
            call mpi_abort(MPI_COMM_WORLD, -1, ierr)
         end if
      end if

      print*,Iam,' F77 MPI sanity test passed.'
      call mpi_finalize(ierr)

      stop
      end
