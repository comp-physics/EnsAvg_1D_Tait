!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Module for MPI setup
!
!  Last update: January 22, 2010
!  Author: Keita Ando
!  Department of Mechanical Engineering
!  Division of Engineering and Applied Science
!  California Institute of Technology, Pasadena CA 91125
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mpi_setup

  USE m_globalvar
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  INTEGER :: mpi_rank
  INTEGER :: mpi_size
  INTEGER :: mpi_err
  INTEGER :: mpi_mn, mpi_mx
  INTEGER :: i_nopad_mn, i_nopad_mx 
  INTEGER :: iup, idown
  INTEGER :: istatus(MPI_STATUS_SIZE)
  INTEGER :: Nx
  INTEGER :: Nx_nopad
  INTEGER :: Nx_tot
  ! no. of padding points for WENO5
  ! padding = wenonum + 1
  INTEGER :: padding

  CONTAINS

  !========================================================================

  SUBROUTINE s_initialize_mpi

    ! mpi_rank = 0 to mpi_size-1
    CALL MPI_INIT( mpi_err )
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD,mpi_size,mpi_err )
    CALL MPI_COMM_RANK( MPI_COMM_WORLD,mpi_rank,mpi_err )

  END SUBROUTINE s_initialize_mpi

  !========================================================================

  SUBROUTINE s_mpi_indices

    ! assume nonperiodic BCs
    INTEGER :: divhigh
    INTEGER :: mpi_mn_nopad, mpi_mx_nopad
    INTEGER :: padminus, padplus
    REAL(KIND(0.D0)) :: tmp

    ! note CEILING(3.7)=4
    tmp = DBLE( Nx_tot )/DBLE( mpi_size )
    divhigh = CEILING( tmp )

    ! min. index for no padding in terms of entire domain
    mpi_mn_nopad = mpi_rank*divhigh + 1
    ! max. index for no padding in terms of entire domain
    mpi_mx_nopad = ( mpi_rank+1 )*divhigh
    ! no. of padding (left/right)
    padding = wenonum + 1
    padminus = padding
    padplus = padding
    IF ( mpi_rank==0 ) padminus = 0
    IF ( mpi_rank==mpi_size-1 ) THEN
       mpi_mx_nopad = Nx_tot
       padplus = 0
    END IF
    ! no. of cells (excluding overlaps) for each processor
    Nx_nopad = mpi_mx_nopad - mpi_mn_nopad + 1

    ! range including overlaps in terms of entire domain
    mpi_mn = mpi_mn_nopad - padminus
    mpi_mx = mpi_mx_nopad + padplus

    ! no. of cells (including overlaps) for each processor
    Nx = mpi_mx - mpi_mn + 1
    ! min. index for no padding for each processor
    i_nopad_mn = 1 + padminus
    ! max. index for no padding for each processor
    i_nopad_mx = Nx - padplus

    ! used for data transfer
    iup = mpi_rank + 1
    idown = mpi_rank - 1
    IF ( mpi_rank==mpi_size-1 ) iup = MPI_PROC_NULL
    IF ( mpi_rank==0 ) idown = MPI_PROC_NULL

  END SUBROUTINE s_mpi_indices

  !========================================================================

  SUBROUTINE s_finalize_mpi

    CALL MPI_FINALIZE( mpi_err )

  END SUBROUTINE s_finalize_mpi

  !========================================================================

END MODULE mpi_setup
