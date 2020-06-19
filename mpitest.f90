PROGRAM hello
   include 'mpif.h'
   integer rank, size, ierror, tag, status(MPI_STATUS_SIZE)
   
   call MPI_INIT(ierror)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
   write(*,*) 'mpi-process with rank', rank, ' ok.'
   call sleep(10)
   write(*,*) 'mpi-process with rank', rank, ' finalizing...'
   call MPI_FINALIZE(ierror)
END PROGRAM
