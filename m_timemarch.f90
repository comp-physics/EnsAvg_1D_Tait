!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Module solving ODEs (semi-discrete form of PDE), given RHS
!
!  Last update: January 21, 2010
!  Author: Keita Ando
!  Department of Mechanical Engineering
!  Division of Engineering and Applied Science
!  California Institute of Technology, Pasadena CA 91125
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE m_timemarch

  USE mpi_setup
  USE m_globalvar
  USE m_rhs
  USE m_misc
  IMPLICIT NONE

  TYPE(timelevel), DIMENSION(:), ALLOCATABLE :: rk
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: rk_uwall

  CONTAINS

  !========================================================================

  SUBROUTINE s_start_timevar

    INTEGER :: iv
    INTEGER :: ir
    INTEGER :: irk

    ! allocate storage
    IF ( timescheme=='euler' ) THEN
       ALLOCATE( rk(Nv) )
       DO iv = 1,Nv
          ALLOCATE( rk(iv)%l(1)%f(1:Nx) )
       END DO
    ELSE IF ( timescheme=='rk3tvd' ) THEN
       ALLOCATE( rk(Nv) )
       DO irk = 1,2
          DO iv = 1,Nv
             ALLOCATE( rk(iv)%l(irk)%f(1:Nx) )
          END DO
       END DO
       ALLOCATE( rk_uwall(2) )
    END IF

  END SUBROUTINE s_start_timevar

  !========================================================================

  SUBROUTINE s_euler

    INTEGER :: iv

    CALL s_rhs( rk%l(1) )
    DO iv = 1,Nv
       rk(iv)%l(1)%f = rk(iv)%l(1)%f + dt*rhs(iv)%f
    END DO
    uwall = uwall + dt*duwall
    CALL s_check_sign( rk%l(1) )

  END SUBROUTINE s_euler

  !========================================================================

  SUBROUTINE s_rk3tvd

    INTEGER :: iv

    ! step 1
    rk_uwall(1) = uwall
    CALL s_rhs( rk%l(1) )
    DO iv = 1,Nv
       rk(iv)%l(2)%f = rk(iv)%l(1)%f + dt*rhs(iv)%f
    END DO
    rk_uwall(2) = rk_uwall(1) + dt*duwall
    CALL s_check_sign( rk%l(2) )

    ! step 2
    uwall = rk_uwall(2)
    CALL s_rhs( rk%l(2) )
    DO iv = 1,Nv
       rk(iv)%l(2)%f = 0.75D0*rk(iv)%l(1)%f &
                     + 0.25D0*( rk(iv)%l(2)%f+dt*rhs(iv)%f )
    END DO
    rk_uwall(2) = 0.75D0*rk_uwall(1) + 0.25D0*( rk_uwall(2)+dt*duwall )
    CALL s_check_sign( rk%l(2) )

    ! step 3
    uwall = rk_uwall(2)
    CALL s_rhs( rk%l(2) )
    DO iv = 1,Nv
       rk(iv)%l(1)%f = third &
                     * ( rk(iv)%l(1)%f+2.D0*(rk(iv)%l(2)%f+dt*rhs(iv)%f) )
    END DO
    uwall = ( rk_uwall(1)+2.D0*(rk_uwall(2)+dt*duwall) )*third
    CALL s_check_sign( rk%l(1) )

  END SUBROUTINE s_rk3tvd

  !========================================================================

  SUBROUTINE s_stop_timevar

    DEALLOCATE( rk )
    IF ( timescheme/='euler' ) DEALLOCATE( rk_uwall )

  END SUBROUTINE s_stop_timevar

  !========================================================================

END MODULE m_timemarch
