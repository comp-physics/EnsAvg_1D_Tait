!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Module for MPI communication
!
!  Last update: December 12, 2008
!  Author: Keita Ando
!  Department of Mechanical Engineering
!  Division of Engineering and Applied Science
!  California Institute of Technology, Pasadena CA 91125
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mpi_transfer

  USE mpi_setup
  USE m_globalvar
  IMPLICIT NONE

  CONTAINS

  !========================================================================

  SUBROUTINE s_mpi_transfer( qgvn )

    TYPE(coordinate), DIMENSION(Nv), INTENT(INOUT) :: qgvn
    INTEGER :: i
    INTEGER :: iv
    INTEGER :: arg
    REAL(KIND(0.D0)), DIMENSION(Nv*padding) :: qsendu
    REAL(KIND(0.D0)), DIMENSION(Nv*padding) :: qsendd
    REAL(KIND(0.D0)), DIMENSION(Nv*padding) :: qrecvd
    REAL(KIND(0.D0)), DIMENSION(Nv*padding) :: qrecvu

    DO iv = 1,Nv
       DO i = 1,padding
          arg = padding*( iv-1 ) + i
          qsendu(arg) = qgvn(iv)%f(Nx-2*padding+i)
       END DO
    END DO
    DO iv = 1,Nv
       DO i = 1,padding
          arg = padding*( iv-1 ) + i
          qsendd(arg) = qgvn(iv)%f(i+padding)
       END DO
    END DO

    ! send to the right and recieve from the left
    CALL MPI_SENDRECV( qsendu(1),Nv*padding,MPI_DOUBLE_PRECISION,iup,0, &
                       qrecvd(1),Nv*padding,MPI_DOUBLE_PRECISION,idown,0, &
                       MPI_COMM_WORLD,istatus,mpi_err )
    ! send to the left and recieve from the right
    CALL MPI_SENDRECV( qsendd(1),Nv*padding,MPI_DOUBLE_PRECISION,idown,0, &
                       qrecvu(1),Nv*padding,MPI_DOUBLE_PRECISION,iup,0, &
                       MPI_COMM_WORLD,istatus,mpi_err )

    ! pad
    IF ( mpi_rank/=0 ) THEN
       DO iv = 1,Nv
          DO i = 1,padding
             arg = padding*( iv-1 ) + i
             qgvn(iv)%f(i) = qrecvd(arg)
          END DO
       END DO
    END IF
    IF ( mpi_rank/=mpi_size-1 ) THEN
       DO iv = 1,Nv
          DO i = 1,padding
             arg = padding*( iv-1 ) + i
             qgvn(iv)%f(Nx-padding+i) = qrecvu(arg)
          END DO
       END DO
    END IF

  END SUBROUTINE s_mpi_transfer

  !========================================================================

END MODULE mpi_transfer
